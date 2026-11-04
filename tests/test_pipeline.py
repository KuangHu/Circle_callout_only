import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from circlestream.common import GB, REPO, digest, fingerprint, load_config, lock, read_json, set_status, validate_sample, write_json
from circlestream.controller import admitted, cleanup, download_file, initialize, reconcile_archive, recover, stage_local, status, submit_one, verify_bundle
from circlestream.fastq import extract, materialize
from circlestream.rules import Context, archive, cf_lines
from scripts.cresil_wgls_checked import execute, strict_call
from types import SimpleNamespace


def fq(path, names=("a", "b"), length=100):
    text = "".join(f"@{name}\n{'A'*length}\n+\n{'I'*length}\n" for name in names)
    if path.suffix == ".gz":
        with gzip.open(path,"wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.ref=self.root/"ref.fa"
        self.ref.write_text(">chr1\n"+"ACGT"*100+"\n")
        self.raw=self.root/"input.fastq.gz"
        fq(self.raw)
        self.c=read_json(REPO/"config.example.json")
        for key in ("state_dir","scratch_dir","archive_dir"):
            self.c[key]=str(self.root/key)
        self.c.update(min_free_gb=0.000001,scratch_budget_gb=1,max_active_samples=2)
        self.c["slurm"].update(account="test",partition="test",qos="test")
        self.s={"sample_id":"test1","accession":"fixture","platform":"nanopore","library":"wgs",
            "reference":str(self.ref),"peak_gb":0.01,"files":[{"name":"reads.fastq.gz","url":"https://example.invalid/reads.fastq.gz",
            "size_bytes":self.raw.stat().st_size,"md5":digest(self.raw,"md5")}]}
        self.manifest=self.root/"samples.json"
        write_json(self.manifest,[self.s])
        initialize(self.c,self.manifest)
        self.state=Path(self.c["state_dir"])/"test1"
        self.scratch=Path(self.c["scratch_dir"])/"test1"

    def tearDown(self):
        self.temp.cleanup()

    def staged(self):
        stage_local(self.c,"test1",[self.raw])
        (self.scratch/"work").mkdir()
        return Context(read_json(self.state/"spec.json"))

    def archived(self):
        x=self.staged()
        ev=x.work/"evidence"
        ev.mkdir()
        (ev/"support.bam").write_bytes(b"mock-bam-for-control-plane-test")
        shutil.copy2(self.ref,ev/"reference.fa")
        write_json(self.state/"calls.json",{"candidate_count":0})
        write_json(self.state/"evidence.json",{"read_ids":0})
        archive(x)
        set_status(self.state,"ARCHIVED")
        return x

    def test_manifest_traversal_rejected(self):
        self.s["sample_id"]="../escape"
        with self.assertRaises(ValueError): validate_sample(self.s)

    def test_manifest_missing_checksum_rejected(self):
        self.s["files"][0]["md5"]=""
        with self.assertRaises(ValueError): validate_sample(self.s)

    def test_wrong_library_rejected(self):
        self.s["library"]="rca"
        with self.assertRaises(ValueError): validate_sample(self.s)

    def test_duplicate_ids_rejected(self):
        write_json(self.manifest,[self.s,self.s])
        with self.assertRaises(ValueError): initialize(self.c,self.manifest)

    def test_init_idempotent_and_mutation_rejected(self):
        initialize(self.c,self.manifest)
        self.s["peak_gb"]=0.02
        write_json(self.manifest,[self.s])
        with self.assertRaises(ValueError): initialize(self.c,self.manifest)

    def test_nested_storage_rejected(self):
        self.c["archive_dir"]=str(Path(self.c["scratch_dir"])/"archive")
        path=self.root/"c.json";write_json(path,self.c)
        with self.assertRaises(ValueError): load_config(path)

    def test_lock_excludes_second_owner(self):
        p=self.root/"lock"
        with lock(p):
            with self.assertRaises(RuntimeError):
                with lock(p): pass
        self.assertFalse(p.exists())

    def test_stage_keeps_source_and_checks_hash(self):
        self.staged()
        self.assertTrue(self.raw.is_file())
        self.assertEqual(status(self.state),"READY")
        self.assertEqual(digest(self.raw),digest(self.scratch/"raw"/"reads.fastq.gz"))

    def test_reservations_include_failed_samples(self):
        self.staged();set_status(self.state,"COMPUTE_FAILED")
        self.c["scratch_budget_gb"]=0.015
        self.assertFalse(admitted(self.c,self.s))

    def test_download_checksum_failure_never_finalizes(self):
        raw=self.root/"download";raw.mkdir()
        item=dict(self.s["files"][0]); item["md5"]="0"*32
        def fake_run(cmd):
            shutil.copyfile(self.raw,cmd[cmd.index("--output")+1])
            return subprocess.CompletedProcess(cmd,0)
        with patch("circlestream.controller.subprocess.run",side_effect=fake_run):
            with self.assertRaises(RuntimeError): download_file(item,raw,1)
        self.assertFalse((raw/item["name"]).exists())

    def test_download_valid_partial_finalizes_without_network(self):
        raw=self.root/"download";raw.mkdir()
        item=self.s["files"][0]
        shutil.copyfile(self.raw,raw/(item["name"]+".part"))
        with patch("circlestream.controller.subprocess.run") as run:
            download_file(item,raw,1)
            run.assert_not_called()
        self.assertEqual(digest(raw/item["name"]),digest(self.raw))

    def test_submission_failure_is_held_not_resubmitted(self):
        self.staged()
        with patch("circlestream.controller.subprocess.run",return_value=subprocess.CompletedProcess([],1,"","error")):
            with self.assertRaises(RuntimeError): submit_one(self.c)
        self.assertEqual(status(self.state),"SUBMITTING")
        self.assertFalse(submit_one(self.c))

    def test_fast_worker_status_not_overwritten(self):
        self.staged()
        def fake_run(*args,**kw):
            set_status(self.state,"ARCHIVED")
            return subprocess.CompletedProcess([],0,"123\n","")
        with patch("circlestream.controller.subprocess.run",side_effect=fake_run):
            self.assertTrue(submit_one(self.c))
        self.assertEqual(status(self.state),"ARCHIVED")
        self.assertEqual(read_json(self.state/"job.json")["job_id"],"123")

    def test_active_job_cannot_be_retried(self):
        self.staged();set_status(self.state,"RUNNING")
        write_json(self.state/"job.json",{"job_id":"123"})
        with patch("circlestream.controller.subprocess.run",return_value=subprocess.CompletedProcess([],0,"123|RUNNING|\n","")):
            with self.assertRaises(ValueError): recover(self.c,"test1")

    def test_missing_output_blocks_archive(self):
        x=self.staged();(x.work/"evidence").mkdir()
        with self.assertRaises(ValueError): archive(x)

    def test_archive_cleanup_and_idempotence(self):
        self.archived(); verify_bundle(self.state)
        cleanup(self.c,"test1");cleanup(self.c,"test1")
        self.assertFalse(self.scratch.exists())
        self.assertTrue(self.raw.exists())
        self.assertEqual(status(self.state),"CLEANED")
        self.assertTrue(Path(read_json(self.state/"archive.json")["path"]).exists())

    def test_corrupted_archive_blocks_cleanup(self):
        self.archived()
        with Path(read_json(self.state/"archive.json")["path"]).open("ab") as handle: handle.write(b"corruption")
        with self.assertRaises(ValueError): cleanup(self.c,"test1")
        self.assertTrue(self.scratch.exists())

    def test_wrong_owner_blocks_cleanup(self):
        self.archived();write_json(self.scratch/".owner.json",{})
        with self.assertRaises(ValueError): cleanup(self.c,"test1")

    def test_running_worker_blocks_cleanup(self):
        self.archived()
        with lock(self.state/".worker-lock"):
            with self.assertRaises(ValueError): cleanup(self.c,"test1")

    def test_qc_pair_order_and_short_reads(self):
        a=self.root/"r1.fq";b=self.root/"r2.fq"
        fq(a,("a/1",));fq(b,("b/2",))
        with self.assertRaises(ValueError): materialize([a,b],self.root/"badpair","illumina")
        fq(a,("a/1",),50);fq(b,("a/2",),50)
        with self.assertRaises(ValueError): materialize([a,b],self.root/"short","illumina")

    def test_qc_truncation(self):
        p=self.root/"truncated.fq";p.write_text("@a\nACGT\n+\nII\n")
        with self.assertRaises(ValueError): materialize([p],self.root/"bad","nanopore")

    def test_duplicate_read_ids_rejected(self):
        p=self.root/"duplicate.fq";fq(p,("a","a"))
        with self.assertRaises(ValueError): materialize([p],self.root/"duplicate","nanopore")

    def test_cresil_known_no_candidate_exit(self):
        def call(args): raise SystemExit("[ABORT] no potential merge region was detected\n")
        execute(SimpleNamespace(main=call),SimpleNamespace(triminput=str(self.root/"trim.txt")))
        self.assertTrue((self.root/"no_candidates.json").is_file())
        self.assertEqual(len((self.root/"eccDNA_final.txt").read_text().splitlines()),1)

    def test_cresil_unknown_or_unmapped_exit_fails(self):
        for message in ("[ABORT] zero trimmed region\n","unexpected failure"):
            def call(args): raise SystemExit(message)
            with self.assertRaises(SystemExit):
                execute(SimpleNamespace(main=call),SimpleNamespace(triminput=str(self.root/"trim.txt")))
        self.assertFalse((self.root/"eccDNA_final.txt").exists())

    def test_cresil_subprocess_failure_propagates(self):
        with patch("scripts.cresil_wgls_checked.subprocess.run",side_effect=subprocess.CalledProcessError(1,"aligner")) as run:
            with self.assertRaises(subprocess.CalledProcessError): strict_call("aligner | sorter",shell=True)
            self.assertEqual(run.call_args.args[0][:4],["bash","-euo","pipefail","-c"])

    def test_archive_reconciliation(self):
        self.archived()
        reconcile_archive(self.c,"test1")
        self.assertEqual(status(self.state),"ARCHIVED")

    def test_archive_wrong_sample_spec_rejected(self):
        self.archived()
        spec=read_json(self.state/"spec.json");spec["sample"]["accession"]="different"
        write_json(self.state/"spec.json",spec)
        with self.assertRaises(ValueError): verify_bundle(self.state)

    def test_extract_keeps_full_quality_and_mates(self):
        a=self.root/"r1.fq";b=self.root/"r2.fq"
        fq(a,("a/1","b/1"));fq(b,("a/2","b/2"))
        out=self.root/"ev";out.mkdir()
        self.assertEqual(extract([a,b],{"a"},out,"illumina"),[1,1])
        with gzip.open(out/"support_reads2.fastq.gz","rt") as handle:
            self.assertIn("@a/2",handle.read())
        with self.assertRaises(ValueError): extract([a,b],{"missing"},out,"illumina")

    def test_pinned_circle_finder_split_no_delete(self):
        mapping=cf_lines("mapping"); calling=cf_lines("call")
        self.assertTrue(mapping[0].startswith("bwa mem"))
        self.assertIn("microDNA-JT.txt",calling[-1])
        self.assertFalse(any(l.startswith("rm ") for l in mapping+calling))
        self.assertFalse(any(l.startswith("bwa mem") for l in calling))


if __name__=="__main__": unittest.main()
