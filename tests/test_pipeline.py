import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from circlestream.common import GB, REPO, S3_PART_SIZES, digest, s3_etag_matches, source_ok, fingerprint, load_config, lock, read_json, set_status, validate_sample, write_json
from circlestream.controller import admitted, cleanup, cleanup_archived, download_file, download_one, initialize, reconcile_archive, reconcile_jobs, recover, stage_local, status, submit_one, verify_bundle, validate_slurm, worker, submit_fetch, fetch_specs, reconcile_fetch, fetch_records
from circlestream.fastq import extract, materialize
from circlestream.rules import Context, archive, cf_lines, convert_sra, copy_caller_evidence
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
        write_json(ev/"reference.json",{"path":str(self.ref),"sha256":digest(self.ref)})
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

    def test_qc_drop_pair_policy(self):
        a=self.root/"r1.fq";b=self.root/"r2.fq"
        # pair "s" has a short R2 only; the whole pair must go, other pairs keep full length.
        a.write_text("@k/1\n"+"A"*150+"\n+\n"+"I"*150+"\n@s/1\n"+"A"*150+"\n+\n"+"I"*150+"\n")
        b.write_text("@k/2\n"+"C"*120+"\n+\n"+"I"*120+"\n@s/2\n"+"C"*40+"\n+\n"+"I"*40+"\n")
        with self.assertRaises(ValueError): materialize([a,b],self.root/"strict","illumina")
        stats=materialize([a,b],self.root/"drop","illumina","drop_pair")
        self.assertEqual((stats["input_pairs"],stats["dropped_pairs"],stats["retained_pairs"]),(2,1,1))
        self.assertEqual((stats["reads"],stats["bases"],stats["min_length"]),(2,270,120))
        self.assertEqual(stats["short_read_policy"],"drop_pair")
        self.assertEqual((self.root/"drop"/"reads2.fastq").read_text().count("\n"),4)

    def test_qc_drop_pair_does_not_mask_structural_errors(self):
        a=self.root/"r1.fq";b=self.root/"r2.fq"
        fq(a,("a/1",),50);fq(b,("b/2",),50)
        with self.assertRaises(ValueError): materialize([a,b],self.root/"order","illumina","drop_pair")
        a.write_text("@a/1\nACGT\n+\nII\n");fq(b,("a/2",),4)
        with self.assertRaises(ValueError): materialize([a,b],self.root/"qual","illumina","drop_pair")
        fq(a,("a/1","b/1"),100);fq(b,("a/2",),100)
        with self.assertRaises(ValueError): materialize([a,b],self.root/"count","illumina","drop_pair")
        fq(a,("a/1",),50);fq(b,("a/2",),50)
        with self.assertRaisesRegex(ValueError,"No reads remain"): materialize([a,b],self.root/"empty","illumina","drop_pair")
        with self.assertRaises(ValueError): materialize([a,b],self.root/"bad",'illumina',"truncate")

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

    def fake_identify_wgls(self, line):
        import importlib.util
        path=self.root/"site"/"cresil"/"cli"/"identify_wgls.py"
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(
            # As with pandas, the KeyError is raised outside identify_wgls.py itself.
            "ns = {}\nexec(compile(\"class Frame(dict):\\n    def __getitem__(self, key): raise KeyError(f\\\"None of [Index({key!r}, dtype='object')] are in the [columns]\\\")\\n\", '/site/pandas/core/frame.py', 'exec'), ns)\nFrame = ns['Frame']\n"
            "def main(args):\n    df_rem_bg_aln_reads_merge = Frame()\n    " + line + "\n")
        spec=importlib.util.spec_from_file_location("fake_identify_wgls",path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module

    def test_cresil_empty_region_keyerror_is_no_candidate(self):
        bodies = ("df_rem_bg_aln_reads_merge = df_rem_bg_aln_reads_merge[['chrom','start','end']]",
                  "df_filtered_breaks_bdg = Frame(); breakpointmerge = 2\n    "
                  "df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name'] >= breakpointmerge]")
        for body in bodies:
            for name in ("no_candidates.json","eccDNA_final.txt"): (self.root/name).unlink(missing_ok=True)
            execute(self.fake_identify_wgls(body),SimpleNamespace(triminput=str(self.root/"trim.txt")))
            self.assertIn("empty-frame",(self.root/"no_candidates.json").read_text())
            self.assertEqual(len((self.root/"eccDNA_final.txt").read_text().splitlines()),1)

    def test_cresil_other_keyerror_fails(self):
        for line in ("other = df_rem_bg_aln_reads_merge[['chrom','start','end']]",
                     "raise KeyError('name')  # df_rem_bg_aln_reads_merge = df_rem_bg_aln_reads_merge[['chrom','start','end']]"):
            module=self.fake_identify_wgls(line)
            with self.assertRaises(KeyError):
                execute(module,SimpleNamespace(triminput=str(self.root/"trim.txt")))
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

    def test_held_samples_release_admission_slots_but_keep_disk_reservation(self):
        self.staged()
        self.c["max_active_samples"] = 1
        self.assertFalse(admitted(self.c, self.s))
        for held in ("COMPUTE_FAILED", "DOWNLOAD_FAILED", "CLEANUP_FAILED"):
            set_status(self.state, held)
            self.assertTrue(admitted(self.c, self.s))
            self.c["scratch_budget_gb"] = 0.015
            self.assertFalse(admitted(self.c, self.s))
            self.c["scratch_budget_gb"] = 1

    def test_killed_job_releases_slot_without_unlocking_or_retrying(self):
        self.staged()
        set_status(self.state, "RUNNING")
        write_json(self.state / "job.json", {"job_id": "123"})
        write_json(self.manifest, [self.s, dict(self.s, sample_id="test2")])
        initialize(self.c, self.manifest)
        stage_local(self.c, "test2", [self.raw])
        self.c["max_jobs"] = 1
        calls = []
        def fake_run(cmd, **kwargs):
            calls.append(cmd[0])
            output = {"sacct": "123|OUT_OF_MEMORY|\n", "squeue": "", "sbatch": "456\n"}[cmd[0]]
            return subprocess.CompletedProcess(cmd, 0, output, "")
        with lock(self.state / ".worker-lock"), patch("circlestream.controller.subprocess.run", side_effect=fake_run):
            self.assertTrue(submit_one(self.c))
            self.assertTrue((self.state / ".worker-lock").exists())
        self.assertEqual(calls, ["sacct", "squeue", "sbatch"])
        self.assertEqual(status(self.state), "COMPUTE_FAILED")
        self.assertTrue(self.scratch.is_dir())
        self.assertEqual(read_json(self.state / "job.json")["job_id"], "123")
        self.assertEqual(read_json(self.state.parent / "test2" / "job.json")["job_id"], "456")

    def test_scheduler_unknown_running_and_requeued_remain_active(self):
        set_status(self.state, "RUNNING")
        write_json(self.state / "job.json", {"job_id": "123"})
        for accounting, queue in (("", ""), ("123|RUNNING|\n", "123\n"),
                                   ("123|PREEMPTED|\n", "123\n"), ("123.batch|FAILED|\n", "")):
            results = [subprocess.CompletedProcess([], 0, accounting, ""),
                       subprocess.CompletedProcess([], 0, queue, "")]
            with patch("circlestream.controller.subprocess.run", side_effect=results):
                reconcile_jobs(self.c)
            self.assertEqual(status(self.state), "RUNNING")
        with patch("circlestream.controller.subprocess.run", side_effect=subprocess.TimeoutExpired("sacct", 30)):
            reconcile_jobs(self.c)
        self.assertEqual(status(self.state), "RUNNING")

    def test_reconciliation_does_not_overwrite_worker_final_status(self):
        set_status(self.state, "RUNNING")
        write_json(self.state / "job.json", {"job_id": "123"})
        def fake_run(cmd, **kwargs):
            set_status(self.state, "ARCHIVED")
            return subprocess.CompletedProcess(cmd, 0, "123|COMPLETED|\n" if cmd[0] == "sacct" else "", "")
        with patch("circlestream.controller.subprocess.run", side_effect=fake_run):
            reconcile_jobs(self.c)
        self.assertEqual(status(self.state), "ARCHIVED")

    def test_one_bad_cleanup_does_not_stop_next_sample(self):
        self.archived()
        write_json(self.manifest, [self.s, dict(self.s, sample_id="test2")])
        initialize(self.c, self.manifest)
        stage_local(self.c, "test2", [self.raw])
        state2 = self.state.parent / "test2"
        scratch2 = self.scratch.parent / "test2"
        shutil.copytree(self.scratch / "work", scratch2 / "work")
        write_json(state2 / "calls.json", {"candidate_count": 0})
        write_json(state2 / "evidence.json", {"read_ids": 0})
        archive(Context(read_json(state2 / "spec.json")))
        set_status(state2, "ARCHIVED")
        # An actual damaged archive, not a mocked cleanup failure.
        bundle = Path(read_json(self.state / "archive.json")["path"])
        good = bundle.read_bytes()
        bundle.write_bytes(good + b"corruption")
        cleanup_archived(self.c)
        self.assertEqual(status(self.state), "CLEANUP_FAILED")
        self.assertTrue(self.scratch.exists())
        self.assertEqual(status(state2), "CLEANED")
        self.assertFalse(scratch2.exists())
        bundle.write_bytes(good)
        reconcile_archive(self.c, "test1")
        cleanup_archived(self.c)
        self.assertEqual(status(self.state), "CLEANED")

    def test_job_setup_commands_are_in_script_before_worker(self):
        self.staged()
        self.c["slurm"]["setup_commands"] = ["source /etc/profile", "module load miniforge3"]
        with patch("circlestream.controller.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "123\n", "")):
            submit_one(self.c)
        script = (self.state / "job.sh").read_text()
        self.assertLess(script.index("module load miniforge3"), script.index("exec "))
        self.c["slurm"]["setup_commands"] = "module load miniforge3"
        with self.assertRaises(ValueError): validate_slurm(self.c)

    def test_illumina_retention_excludes_whole_library_tables(self):
        out = self.root / "caller"; out.mkdir()
        ev = self.root / "evidence"; ev.mkdir()
        keep = ["sample-ref.microDNA-JT.txt", "sample-ref.split.txt", "sample-ref.split.id-freq.txt",
                "sample-ref.split_freq2.oneline.txt", "sample-ref.disc.txt"]
        omit = ["sample-ref.concordant.txt", "sample-ref.concordant.id-freq.txt",
                "sample-ref.concordant_freqGr3.txt", "sample-ref.sorted.bam", "unrelated.txt"]
        for name in keep + omit: (out / name).write_text("records\n")
        copy_caller_evidence(out, ev, "illumina", set())
        self.assertTrue(all((ev / name).is_file() for name in keep))
        self.assertFalse(any((ev / name).exists() for name in omit))
        retained = read_json(ev / "caller-retention.json")
        self.assertEqual({p["path"] for p in retained["copied_files"]}, set(keep))
        self.assertGreater(retained["caller_bytes_not_copied_verbatim"], 0)

    def test_ont_retention_keeps_graphs_and_only_support_trim_rows(self):
        out = self.root / "caller"; out.mkdir()
        ev = self.root / "evidence"; ev.mkdir()
        keep = ["eccDNA_final.txt", "cresil_run/subGraphs.summary.txt",
                "cresil_run/assemGraph/ecc1/final_reads.fa", "cresil_run/assemGraph/ecc1/reference_regions.fa"]
        omit = ["cresil_run/tmp/merge_genomecov.txt", "cresil_run/tmp/ref_aln.regions.bed.gz",
                "cresil_run/tmp/ref_aln.bam", "unrecognized.txt"]
        for name in keep + omit:
            p = out / name; p.parent.mkdir(parents=True, exist_ok=True); p.write_text("data\n")
        (out / "trim.txt").write_text("readid\tq_start\tq_end\na\t0\t100\nb\t0\t100\na\t100\t200\n")
        copy_caller_evidence(out, ev, "nanopore", {"a"})
        self.assertTrue(all((ev / "cresil" / name).is_file() for name in keep))
        self.assertFalse(any((ev / "cresil" / name).exists() for name in omit + ["trim.txt"]))
        self.assertEqual((ev / "cresil" / "trim.support.txt").read_text(), "readid\tq_start\tq_end\na\t0\t100\na\t100\t200\n")
        self.assertEqual(read_json(ev / "caller-retention.json")["support_trim_rows"], 2)

    def test_ont_zero_support_retains_trim_header(self):
        out = self.root / "caller"; out.mkdir()
        ev = self.root / "evidence"; ev.mkdir()
        (out / "trim.txt").write_text("readid\tq_start\na\t0\n")
        copy_caller_evidence(out, ev, "nanopore", set())
        self.assertEqual((ev / "cresil" / "trim.support.txt").read_text(), "readid\tq_start\n")

    def test_two_sample_download_submit_archive_cleanup_cycle(self):
        # Exercise real control/archiving code with fake HTTPS transport and Slurm.
        # Fake caller products do NOT test biological correctness.
        write_json(self.manifest, [self.s, dict(self.s, sample_id="test2")])
        initialize(self.c, self.manifest)
        self.c["max_active_samples"] = 1
        self.c["max_jobs"] = 1
        def fake_run(cmd, **kwargs):
            if cmd[0] == "curl":
                shutil.copyfile(self.raw, cmd[cmd.index("--output") + 1])
                return subprocess.CompletedProcess(cmd, 0)
            self.assertEqual(cmd[0], "sbatch")
            return subprocess.CompletedProcess(cmd, 0, "123\n", "")
        with patch("circlestream.controller.socket.getfqdn", return_value=self.c["download_host"]), \
             patch("circlestream.controller.subprocess.run", side_effect=fake_run):
            for sid in ("test1", "test2"):
                self.assertTrue(download_one(self.c))
                self.assertFalse(download_one(self.c))
                self.assertTrue(submit_one(self.c))
                state = self.state.parent / sid
                x = Context(read_json(state / "spec.json"))
                ev = x.work / "evidence"; ev.mkdir(parents=True)
                (ev / "support.bam").write_bytes(b"fake-caller-product")
                write_json(ev / "reference.json", {"path": str(self.ref), "sha256": digest(self.ref)})
                write_json(state / "calls.json", {"candidate_count": 0})
                write_json(state / "evidence.json", {"read_ids": 0})
                archive(x)
                set_status(state, "ARCHIVED")
                cleanup_archived(self.c)
                self.assertEqual(status(state), "CLEANED")
                self.assertFalse(x.scratch.exists())
                verify_bundle(state)
            self.assertFalse(download_one(self.c))
            self.assertFalse(submit_one(self.c))


    # ---- NCBI SRA on AWS source + in-job download ----
    def sra_item(self, path):
        data = path.read_bytes(); part = S3_PART_SIZES[0]
        sums = b"".join(hashlib.md5(data[i:i+part]).digest() for i in range(0, len(data), part))
        n = -(-len(data) // part)
        return {"name": "SRR000001.sra", "url": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001",
                "size_bytes": len(data), "s3_etag": f"{hashlib.md5(sums).hexdigest()}-{n}"}

    def test_s3_multipart_etag_verifies_and_detects_corruption(self):
        obj = self.root / "SRR000001.sra"
        obj.write_bytes(bytes(range(256)) * (S3_PART_SIZES[0] // 256 * 2 + 7))
        item = self.sra_item(obj)
        self.assertTrue(item["s3_etag"].endswith("-3"))
        self.assertTrue(source_ok(obj, item))
        data = bytearray(obj.read_bytes()); data[S3_PART_SIZES[0] + 5] ^= 1; obj.write_bytes(bytes(data))
        self.assertFalse(source_ok(obj, item))
        self.assertFalse(s3_etag_matches(obj, "not-an-etag"))

    def test_sra_manifest_entries_validated(self):
        obj = self.root / "SRR000001.sra"; obj.write_bytes(b"x" * 1000)
        for platform in ("illumina", "nanopore"):
            validate_sample(dict(self.s, platform=platform, files=[self.sra_item(obj)]))
        bad = self.sra_item(obj); del bad["s3_etag"]
        with self.assertRaises(ValueError): validate_sample(dict(self.s, files=[bad]))
        with self.assertRaises(ValueError):
            validate_sample(dict(self.s, platform="illumina", files=[self.sra_item(obj), dict(self.sra_item(obj), name="SRR000002.sra")]))

    def slurm_mode(self, n=1):
        c = dict(self.c, download_mode="slurm", auto_cleanup=True, fetch={"max_jobs": 1, "batch_gb": 1, "time": "01:00:00", "memory_gb": 1})
        for key in ("state_dir", "scratch_dir", "archive_dir"):
            c[key] = str(self.root / ("sl_" + key))
        write_json(self.manifest, [dict(self.s, sample_id=f"t{i}") for i in range(n)])
        initialize(c, self.manifest)
        return c

    def fake_sbatch(self, log, jid="555"):
        def run(cmd, **kwargs):
            if cmd[0] == "sbatch":
                log.append(cmd); return subprocess.CompletedProcess(cmd, 0, jid + "\n", "")
            if cmd[0] == "curl":
                shutil.copyfile(self.raw, cmd[cmd.index("--output") + 1]); return subprocess.CompletedProcess(cmd, 0)
            raise AssertionError(cmd)
        return run

    def test_slurm_mode_fetch_job_is_one_cpu_and_groups_samples(self):
        c = self.slurm_mode(3)
        log = []
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch(log)), \
             patch("circlestream.controller.reconcile_fetch"):
            with self.assertRaises(RuntimeError): download_one(c)
            self.assertTrue(submit_fetch(c))
            self.assertFalse(submit_fetch(c))  # fetch.max_jobs=1
        cmd = log[0]
        self.assertEqual(cmd[cmd.index("--cpus-per-task") + 1], "1")
        rec = read_json(next((Path(c["state_dir"]) / ".fetch").glob("fetch-*.json")))
        self.assertEqual(rec["samples"], ["t0", "t1"])  # max_active_samples=2
        self.assertEqual([status(Path(c["state_dir"]) / f"t{i}") for i in range(3)], ["DOWNLOADING", "DOWNLOADING", "NEW"])

    def test_slurm_mode_fetch_downloads_then_process_job_submitted(self):
        c = self.slurm_mode(1)
        state = Path(c["state_dir"]) / "t0"
        log = []
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch(log)), \
             patch("circlestream.controller.reconcile_fetch"), patch("circlestream.controller.reconcile_jobs"):
            submit_fetch(c)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "555"}):
                self.assertEqual(fetch_specs([state / "spec.json"]), 0)
            self.assertEqual(status(state), "READY")
            self.assertTrue(read_json(state / "download.json")["in_job"])
            self.assertTrue(submit_one(c))
        self.assertEqual(log[-1][log[-1].index("--cpus-per-task") + 1], str(c["slurm"]["cpus"]))

    def test_slurm_mode_ended_fetch_job_holds_unfinished_samples(self):
        c = self.slurm_mode(1)
        state = Path(c["state_dir"]) / "t0"
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"):
            submit_fetch(c)
        def sched(cmd, **kw):
            return subprocess.CompletedProcess(cmd, 0, "555|TIMEOUT\n" if cmd[0] == "sacct" else "", "")
        with patch("circlestream.controller.subprocess.run", side_effect=sched):
            reconcile_fetch(c)
        self.assertEqual(status(state), "DOWNLOAD_FAILED")
        self.assertEqual(fetch_records(c), [])
        recover(c, "t0")
        self.assertEqual(status(state), "NEW")  # next fetch job resumes into the same owned scratch

    def test_slurm_mode_admission_uses_source_bytes_without_tree_walk(self):
        c = self.slurm_mode(1)
        c["scratch_budget_gb"] = self.raw.stat().st_size / GB / 2
        with patch("circlestream.controller.tree_size", side_effect=AssertionError("Lustre walk")):
            self.assertFalse(admitted(c, self.s))
            c["scratch_budget_gb"] = 1
            self.assertTrue(admitted(c, self.s))

    def test_slurm_mode_worker_uses_local_disk_and_cleans_up(self):
        c = self.slurm_mode(1)
        state = Path(c["state_dir"]) / "t0"; scratch = Path(c["scratch_dir"]) / "t0"
        local = self.root / "nodetmp"; local.mkdir()
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"), patch("circlestream.controller.reconcile_jobs"):
            submit_fetch(c)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "555"}):
                fetch_specs([state / "spec.json"])
            submit_one(c)
        seen = {}
        def fake_snakemake(cmd, **kw):
            x = Context(read_json(state / "spec.json"))
            seen["work"] = x.work
            self.assertTrue(str(x.work).startswith(str(local / "circle-stream-777")))
            ev = x.work / "evidence"; ev.mkdir(parents=True)
            (ev / "support.bam").write_bytes(b"fake"); write_json(ev / "reference.json", {"path": str(self.ref), "sha256": digest(self.ref)})
            write_json(state / "calls.json", {"candidate_count": 0}); write_json(state / "evidence.json", {"read_ids": 0})
            archive(x)
            return subprocess.CompletedProcess(cmd, 0)
        with patch.dict("os.environ", {"SLURM_JOB_ID": "777", "CS_LOCAL_ROOT": str(local)}), \
             patch("circlestream.controller.subprocess.run", side_effect=fake_snakemake):
            worker(state / "spec.json", "archive")
        self.assertEqual(status(state), "CLEANED")
        self.assertFalse(scratch.exists())
        self.assertFalse((local / "circle-stream-777").exists())
        verify_bundle(state)

    def test_slurm_mode_failed_worker_still_removes_local_disk(self):
        c = self.slurm_mode(1)
        state = Path(c["state_dir"]) / "t0"
        local = self.root / "nodetmp"; local.mkdir()
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"), patch("circlestream.controller.reconcile_jobs"):
            submit_fetch(c)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "555"}):
                fetch_specs([state / "spec.json"])
            submit_one(c)
        with patch.dict("os.environ", {"SLURM_JOB_ID": "778", "CS_LOCAL_ROOT": str(local)}), \
             patch("circlestream.controller.subprocess.run", side_effect=subprocess.CalledProcessError(1, "snakemake")):
            with self.assertRaises(subprocess.CalledProcessError): worker(state / "spec.json", "archive")
        self.assertEqual(status(state), "COMPUTE_FAILED")
        self.assertFalse((local / "circle-stream-778").exists())
        self.assertTrue((Path(c["scratch_dir"]) / "t0" / "raw" / "reads.fastq.gz").is_file())  # source kept for recovery

    def test_slurm_mode_auto_cleanup_false_keeps_staged_source(self):
        c = self.slurm_mode(1); c["auto_cleanup"] = False
        # config snapshot must match: re-init a separate run root with the flag off
        for key in ("state_dir", "scratch_dir", "archive_dir"):
            c[key] = c[key] + "_keep"
        initialize(c, self.manifest)
        state = Path(c["state_dir"]) / "t0"; scratch = Path(c["scratch_dir"]) / "t0"
        local = self.root / "nodetmp"; local.mkdir()
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"), patch("circlestream.controller.reconcile_jobs"):
            submit_fetch(c)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "555"}):
                fetch_specs([state / "spec.json"])
            submit_one(c)
        def fake_snakemake(cmd, **kw):
            x = Context(read_json(state / "spec.json"))
            ev = x.work / "evidence"; ev.mkdir(parents=True)
            (ev / "support.bam").write_bytes(b"fake")
            write_json(ev / "reference.json", {"path": str(self.ref), "sha256": digest(self.ref)})
            write_json(state / "calls.json", {"candidate_count": 0}); write_json(state / "evidence.json", {"read_ids": 0})
            archive(x)
            return subprocess.CompletedProcess(cmd, 0)
        with patch.dict("os.environ", {"SLURM_JOB_ID": "779", "CS_LOCAL_ROOT": str(local)}), \
             patch("circlestream.controller.subprocess.run", side_effect=fake_snakemake):
            worker(state / "spec.json", "archive")
        self.assertEqual(status(state), "ARCHIVED")
        self.assertTrue((scratch / "raw" / "reads.fastq.gz").is_file())
        self.assertFalse((local / "circle-stream-779").exists())  # node-local work is always removed

    def test_slurm_mode_fetch_groups_are_balanced_and_capped(self):
        c = self.slurm_mode(6)
        c.update(max_active_samples=6, fetch=dict(c["fetch"], max_jobs=3, batch_gb=100))
        records = []
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"):
            while submit_fetch(c):
                pass
        sizes = [len(read_json(r)["samples"]) for r in fetch_records(c)]
        self.assertEqual(sizes, [2, 2, 2])  # 6 admissible over 3 free slots, not 6 in one job
        c2 = self.slurm_mode(5)
        c2.update(max_active_samples=5, fetch=dict(c2["fetch"], max_jobs=1, batch_gb=100, max_samples=3))
        for key in ("state_dir", "scratch_dir", "archive_dir"):
            c2[key] = c2[key] + "_cap"
        initialize(c2, self.manifest)
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"):
            submit_fetch(c2)
        self.assertEqual(len(read_json(fetch_records(c2)[0])["samples"]), 3)

    def test_reference_store_is_deduplicated_read_only_and_checked(self):
        c = dict(self.c, reference_store=str(self.root / "refstore"))
        for key in ("state_dir", "scratch_dir", "archive_dir"):
            c[key] = str(self.root / ("rs_" + key))
        write_json(self.manifest, [dict(self.s, sample_id="r1"), dict(self.s, sample_id="r2")])
        initialize(c, self.manifest)
        sha = digest(self.ref)
        stored = Path(c["reference_store"]) / f"{sha}.fa"
        self.assertEqual(sorted(p.name for p in stored.parent.iterdir()), [stored.name])
        self.assertEqual(digest(stored), sha)
        self.assertFalse(stored.stat().st_mode & 0o222)
        self.assertEqual(read_json(Path(c["state_dir"]) / "r1" / "spec.json")["reference_store_path"], str(stored))
        os.chmod(stored, 0o644); stored.write_text(">bad\nA\n")
        write_json(self.manifest, [dict(self.s, sample_id="r3")])
        with self.assertRaises(ValueError): initialize(c, self.manifest)

    def test_reference_fasta_checks(self):
        from scripts.prepare_references import check_fasta
        norm, hmap, n, total = check_fasta(b">c1 len=4 x\nACGT\n>c2\nNNAC\n")
        self.assertEqual((n, total, hmap), (2, 8, None))
        self.assertEqual(norm, b">c1 len=4 x\nACGT\n>c2\nNNAC\n")  # safe IDs kept verbatim
        norm, hmap, _, _ = check_fasta(b">ENA|X1|X1.1 desc\nACGT\n")
        self.assertEqual(norm, b">ENA_X1_X1.1\nACGT\n"); self.assertEqual(hmap, {"ENA_X1_X1.1": "ENA|X1|X1.1 desc"})
        for bad in (b">a\nACGT\n>a\nACGT\n", b">a\nAC1T\n", b">a\n>b\nACGT\n", b"ACGT\n", b">a\n"):
            with self.assertRaises(ValueError): check_fasta(bad)

    def sra_context(self, platform):
        obj = self.root / "SRR000001.sra"; obj.write_bytes(b"sra")
        x = SimpleNamespace(s={"platform": platform}, cpus=4, state=self.root / "st", scratch=self.root / "sc")
        x.state.mkdir(exist_ok=True); (x.scratch / "raw").mkdir(parents=True, exist_ok=True)
        x.raw_files = [obj]
        fqd = x.scratch / "raw" / "fastq"
        x.inputs = [fqd / "SRR000001_1.fastq", fqd / "SRR000001_2.fastq"] if platform == "illumina" else [fqd / "SRR000001.fastq"]
        return x

    def fake_fqd(self, names):
        def run(cmd, **kwargs):
            out = Path(cmd[cmd.index("-O") + 1])
            for n in names:
                fq(out / n)
            return subprocess.CompletedProcess(cmd, 0, "", "")
        return run

    def test_convert_sra_illumina_drops_orphans_and_records(self):
        x = self.sra_context("illumina")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.subprocess.run", side_effect=self.fake_fqd(["SRR000001_1.fastq", "SRR000001_2.fastq", "SRR000001.fastq"])):
            info = convert_sra(x)
        self.assertEqual(info["dropped_unpaired_files"], ["SRR000001.fastq"])
        self.assertTrue(all(p.is_file() for p in x.inputs))
        self.assertFalse((x.scratch / "raw" / "fastq" / "SRR000001.fastq").exists())
        with patch("circlestream.rules.subprocess.run") as run:
            convert_sra(x); run.assert_not_called()  # completed conversion is reused

    def test_convert_sra_layout_mismatch_fails(self):
        x = self.sra_context("illumina")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.subprocess.run", side_effect=self.fake_fqd(["SRR000001.fastq"])):
            with self.assertRaises(ValueError): convert_sra(x)
        x = self.sra_context("nanopore")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.subprocess.run", side_effect=self.fake_fqd(["SRR000001_1.fastq", "SRR000001_2.fastq"])):
            with self.assertRaises(ValueError): convert_sra(x)


if __name__=="__main__": unittest.main()
