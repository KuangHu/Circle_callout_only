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
from circlestream.controller import admitted, release_failed_raw, only_held_left, cleanup, cleanup_archived, download_file, download_one, initialize, reconcile_archive, reconcile_jobs, recover, stage_local, status, submit_one, verify_bundle, validate_slurm, worker, submit_fetch, fetch_specs, reconcile_fetch, fetch_records
from circlestream.fastq import extract, materialize
from circlestream.rules import Context, archive, cf_lines, convert_sra, copy_caller_evidence, library_depth
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
        keep = ["eccDNA_final.txt", "cresil_run/subGraphs.summary.txt", "cresil_run/assemGraph/ecc1/reference_regions.fa",
                "ont_junctions.tsv", "ont_junction_events.tsv", "ont_junction_ids.txt", "ont_junctions.support.paf.gz",
                "cresil_trace.tsv", "ont_library.json", "contig_depth.tsv"]
        # final_reads.fa: read sequence for every graph (whole replicons in deep samples); IDs are kept elsewhere.
        omit = ["cresil_run/tmp/merge_genomecov.txt", "cresil_run/tmp/ref_aln.regions.bed.gz",
                "cresil_run/tmp/ref_aln.bam", "unrecognized.txt", "cresil_run/assemGraph/ecc1/final_reads.fa", "junctions.paf"]
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

    def test_slurm_mode_release_failed_raw_frees_budget_and_recover_refetches(self):
        c = self.slurm_mode(2); c["fetch"]["max_samples"] = 1; c["release_failed_raw"] = True
        state, scratch = Path(c["state_dir"]) / "t0", Path(c["scratch_dir"]) / "t0"
        with patch("circlestream.controller.subprocess.run", side_effect=self.fake_sbatch([])), \
             patch("circlestream.controller.reconcile_fetch"):
            submit_fetch(c)
            with patch.dict("os.environ", {"SLURM_JOB_ID": "555"}):
                fetch_specs([state / "spec.json"])
        write_json(state / "job.json", {"job_id": "777"})
        set_status(state, "COMPUTE_FAILED", error="scratch reservation exceeded")
        c["scratch_budget_gb"] = 1.5 * self.raw.stat().st_size / GB
        t1 = read_json(Path(c["state_dir"]) / "t1" / "spec.json")["sample"]
        self.assertFalse(admitted(c, t1))           # held raw bytes still count against the budget
        self.assertTrue(only_held_left(dict(c)) is False)
        queue = lambda out: (lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, out, ""))  # noqa: E731
        with patch("circlestream.controller.subprocess.run", side_effect=queue("777\n")):
            release_failed_raw(c)
        self.assertTrue(scratch.exists())           # job still in the queue: untouched
        with patch("circlestream.controller.subprocess.run", side_effect=queue("")):
            release_failed_raw(c)
        self.assertFalse(scratch.exists())
        info = read_json(state / "status.json")
        self.assertEqual((info["status"], info["raw_released"], info["error"]), ("COMPUTE_FAILED", True, "scratch reservation exceeded"))
        self.assertTrue(admitted(c, t1))
        recover(c, "t0")
        self.assertEqual((status(state), read_json(state / "status.json")["raw_refetch"]), ("NEW", True))

    def test_release_failed_raw_is_opt_in(self):
        c = self.slurm_mode(1)
        state, scratch = Path(c["state_dir"]) / "t0", Path(c["scratch_dir"]) / "t0"
        scratch.mkdir(); set_status(state, "COMPUTE_FAILED")
        with patch("circlestream.controller.subprocess.run", side_effect=AssertionError("no queue query")):
            release_failed_raw(c)
        self.assertTrue(scratch.exists())

    def test_watch_stops_when_only_held_samples_remain(self):
        c = self.slurm_mode(2)
        set_status(Path(c["state_dir"]) / "t0", "CLEANED")
        self.assertFalse(only_held_left(c))         # t1 is NEW
        set_status(Path(c["state_dir"]) / "t1", "COMPUTE_FAILED")
        self.assertTrue(only_held_left(c))
        set_status(Path(c["state_dir"]) / "t1", "CLEANED")
        self.assertFalse(only_held_left(c))         # all cleaned: the normal exit, not a held stop

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
        x = SimpleNamespace(s={"platform": platform}, c={}, cpus=4, state=self.root / "st", scratch=self.root / "sc")
        x.state.mkdir(exist_ok=True); (x.scratch / "raw").mkdir(parents=True, exist_ok=True)
        x.raw_files = [obj]
        fqd = x.scratch / "raw" / "fastq"
        x.inputs = [fqd / "SRR000001_1.fastq", fqd / "SRR000001_2.fastq"] if platform == "illumina" else [fqd / "SRR000001.fastq"]
        return x

    def fake_fqd(self, names):
        def run(cmd, timeout_s, cwd=None):
            assert cwd == Path(cmd[-1]).parent  # refseqs resolve from the working directory
            out = Path(cmd[cmd.index("-O") + 1])
            for n in names:
                fq(out / n)
            return 0, "", "", False
        return run

    def test_convert_sra_illumina_drops_orphans_and_records(self):
        x = self.sra_context("illumina")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.run_limited", side_effect=self.fake_fqd(["SRR000001_1.fastq", "SRR000001_2.fastq", "SRR000001.fastq"])):
            info = convert_sra(x)
        self.assertEqual(info["dropped_unpaired_files"], ["SRR000001.fastq"])
        self.assertTrue(all(p.is_file() for p in x.inputs))
        self.assertFalse((x.scratch / "raw" / "fastq" / "SRR000001.fastq").exists())
        with patch("circlestream.rules.run_limited") as run:
            convert_sra(x); run.assert_not_called()  # completed conversion is reused

    def test_convert_sra_layout_mismatch_fails(self):
        x = self.sra_context("illumina")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.run_limited", side_effect=self.fake_fqd(["SRR000001.fastq"])):
            with self.assertRaises(ValueError): convert_sra(x)
        x = self.sra_context("nanopore")
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.run_limited", side_effect=self.fake_fqd(["SRR000001_1.fastq", "SRR000001_2.fastq"])):
            with self.assertRaises(ValueError): convert_sra(x)

    def test_convert_sra_timeout_and_network_isolation(self):
        x = self.sra_context("illumina"); x.c = {"sra_conversion_timeout_minutes": 2}
        seen = {}
        def hang(cmd, timeout_s, cwd=None):
            seen.update(cmd=cmd, timeout=timeout_s)
            return -9, "", "killed", True
        with patch("circlestream.rules.sra_tool", return_value="fasterq-dump"), \
             patch("circlestream.rules.network_isolation", return_value=["unshare", "-rn"]), \
             patch("circlestream.rules.run_limited", side_effect=hang):
            with self.assertRaisesRegex(RuntimeError, "sra_conversion_timeout_minutes"):
                convert_sra(x)
        self.assertEqual(seen["timeout"], 120)
        self.assertEqual(seen["cmd"][:3], ["unshare", "-rn", "fasterq-dump"])

    def test_run_limited_kills_process_group(self):
        from circlestream.sra import run_limited
        rc, _, _, timed_out = run_limited(["bash", "-c", "sleep 30 & sleep 30; wait"], 1)
        self.assertTrue(timed_out)

    def test_refseq_deps_parsing_and_fetch(self):
        from circlestream import sra
        raw = self.root / "raw"; raw.mkdir(); obj = raw / "ERR1.sra"; obj.write_bytes(b"x")
        calls = []
        outputs = iter(["BX571856.1,gi|49240382|emb|BX571856.1|,false,remote::https://sra-download.ncbi.nlm.nih.gov/traces/refseq/BX571856.1\n",
                        "BX571856.1,gi|49240382|emb|BX571856.1|,false,remote:" + str(raw / "BX571856.1") + "\n"])
        def fake_run(cmd, **kw):
            calls.append(cmd)
            if cmd[0].endswith("align-info"):
                self.assertEqual(kw.get("cwd"), raw)
            if cmd[0].endswith("align-info"):
                return subprocess.CompletedProcess(cmd, 0, next(outputs), "")
            if cmd[0] == "curl":
                Path(cmd[cmd.index("--output") + 1]).write_bytes(b"refseq")
                return subprocess.CompletedProcess(cmd, 0, "", "")
            return subprocess.CompletedProcess(cmd, 0, "Table 'BX571856.1' is consistent", "")  # vdb-validate
        with patch("circlestream.sra.sra_tool", side_effect=lambda n: n), patch("circlestream.sra.subprocess.run", side_effect=fake_run):
            recs = sra.fetch_refseqs(obj, 2)
        self.assertEqual(recs[0]["seqid"], "BX571856.1"); self.assertEqual(recs[0]["size_bytes"], 6)
        self.assertTrue((raw / "BX571856.1").is_file()); self.assertFalse((raw / "BX571856.1.part").exists())
        self.assertIn("https://sra-download.ncbi.nlm.nih.gov/traces/refseq/BX571856.1", [c[-1] for c in calls if c[0] == "curl"])
        # Still remote after download -> hold the sample instead of letting conversion hang.
        outputs = iter(["A1.1,a,false,remote::https://sra-download.ncbi.nlm.nih.gov/traces/refseq/A1.1\n"] * 2)
        with patch("circlestream.sra.sra_tool", side_effect=lambda n: n), patch("circlestream.sra.subprocess.run", side_effect=fake_run):
            with self.assertRaisesRegex(RuntimeError, "unresolved"):
                sra.fetch_refseqs(obj, 1)
        # Unaligned runs: align-info prints only a warning, no dependencies.
        with patch("circlestream.sra.sra_tool", side_effect=lambda n: n), \
             patch("circlestream.sra.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "", "warn: not a database")):
            self.assertEqual(sra.fetch_refseqs(obj, 1), [])

    def test_unexpected_refseq_location_rejected(self):
        from circlestream import sra
        raw = self.root / "raw"; raw.mkdir(); obj = raw / "ERR1.sra"; obj.write_bytes(b"x")
        with patch("circlestream.sra.sra_tool", side_effect=lambda n: n), \
             patch("circlestream.sra.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "X1.1,a,false,remote::http://evil/X1.1\n", "")):
            with self.assertRaisesRegex(ValueError, "unexpected location"):
                sra.fetch_refseqs(obj, 1)

    def test_fetch_inputs_records_refseqs_for_sra(self):
        from circlestream.controller import fetch_inputs
        state = self.root / "st"; scratch = self.root / "sc"; state.mkdir(); (scratch / "raw").mkdir(parents=True)
        sample = {"files": [{"name": "ERR1.sra", "url": "https://x/ERR1", "size_bytes": 1, "s3_etag": "0" * 32}]}
        with patch("circlestream.controller.download_file"), \
             patch("circlestream.controller.fetch_refseqs", return_value=[{"seqid": "BX571856.1"}]) as fr:
            fetch_inputs(self.c, sample, state, scratch)
        fr.assert_called_once()
        self.assertEqual(read_json(state / "download.json")["refseqs"], [{"seqid": "BX571856.1"}])

    def grep_binary(self):
        for cand in [Path(os.environ.get("CONDA_PREFIX", "/nonexistent")) / "bin" / "grep", shutil.which("grep")]:
            if cand and Path(cand).is_file():
                out = subprocess.run([str(cand), "--version"], capture_output=True, text=True).stdout
                if "GNU grep" in out:
                    return str(cand)
        self.skipTest("GNU grep not available")

    def wordgrep(self, pats, table, env=None):
        p = self.root / "pats.txt"; t = self.root / "table.txt"
        p.write_bytes(pats); t.write_bytes(table)
        ours = subprocess.run([os.sys.executable, str(REPO / "scripts" / "wordgrep.py"), str(p), str(t)], capture_output=True,
                              env=dict(os.environ, **(env or {})))
        ref = subprocess.run([self.grep_binary(), "-w", "-F", "-f", str(p), str(t)], capture_output=True, env=dict(os.environ, LC_ALL="C"))
        return ours, ref

    def test_wordgrep_matches_gnu_grep_on_randomized_tables(self):
        import random
        rng = random.Random(7)
        for trial in range(40):
            ids = [f"ERR{rng.randint(1, 3)}.{rng.randint(1, 300)}" for _ in range(60)] + [str(rng.randint(1, 500)) for _ in range(5)] \
                + [f"r{rng.randint(1, 50)}-x:{rng.randint(1, 9)}" for _ in range(5)] + ["NC_1.2", "101M", "+", "first"]
            rows = []
            for _ in range(400):
                rid = rng.choice(ids[:-4] + [f"ERR1.{rng.randint(1, 3000)}", f"ERR2.{rng.randint(1, 30)}.5", "a_b.c"])
                cigar = rng.choice(["101M", "50M51S", "30S71M", "60H41M"])
                rows.append(f"{rng.choice(['NC_1.2', 'chr1', 'contig_7.1'])}\t{rng.randint(0, 600)}\t{rng.randint(0, 600)}\t{rid}\t"
                            f"{rng.choice('12')}\t{rng.randint(0, 60)}\t{rng.choice('+-')}\t{cigar}" + rng.choice(["", "\tfirst", "\tsecond"]))
            pats = rng.sample(ids, rng.randint(0, 30))
            table = ("\n".join(rows) + rng.choice(["\n", ""])).encode()
            ours, ref = self.wordgrep(("\n".join(pats) + ("\n" if pats else "")).encode(), table)
            self.assertEqual(ours.returncode, 0, ours.stderr)
            self.assertEqual(ours.stdout, ref.stdout, f"trial {trial}")

    def test_wordgrep_edge_cases(self):
        cases = [
            (b"", b"a\tb\n"),                          # no patterns -> nothing
            (b"\n", b"a\tb\n\n"),                      # one empty pattern -> grep's own behaviour
            (b"a\n\nb\n", b"x\ty\na\n"),               # empty pattern among others
            (b"12\n", b"c\t12\t5\tERR1.3\n"),          # numeric ID also matching a coordinate column
            (b"ERR1.1\n", b"c\t1\t2\tERR1.12\nc\t1\t2\tERR1.1\n"),  # prefix must not match
            (b"ERR1.1\n", b"c\t1\t2\tERR1.1.5\n"),     # '.' is a word boundary
            (b"a b\n", b"x\ta b\tz\n"),                # pattern with a space
            (b"ERR1.1", b"c\t1\t2\tERR1.1"),           # no trailing newlines
        ]
        for pats, table in cases:
            ours, ref = self.wordgrep(pats, table)
            self.assertIn(ours.returncode, (0,), (pats, table, ours.stderr))
            self.assertEqual(ours.stdout, ref.stdout, (pats, table))

    def test_wordgrep_verify_mode(self):
        ours, ref = self.wordgrep(b"ERR1.1\n", b"c\t1\t2\tERR1.1\n", env={"CS_WORDGREP_VERIFY": "1"})
        self.assertEqual(ours.returncode, 0); self.assertIn(b"IDENTICAL", ours.stderr); self.assertEqual(ours.stdout, ref.stdout)

    def test_circle_finder_adapter_replaces_every_word_grep(self):
        from circlestream.rules import WORD_GREP
        lines = cf_lines("call")
        greps = [l for l in lines if "grep" in l]
        self.assertEqual(len(greps), 5)
        self.assertTrue(all(WORD_GREP.fullmatch(l) for l in greps))

    def junctions(self, paf_rows, contigs={"chr1": 30000, "pl": 3000}, extra=(), seqs=None, extra_paf=()):
        """paf_rows: (read, qlen, [(qs, qe, strand, contig, ts, te[, mapq[, tp]])]). Returns (rows, events, ids)."""
        import csv as _csv, random as _random
        paf = self.root / "j.paf"; ref = self.root / "j.fa"; out = self.root / "j.tsv"; ev = self.root / "j.events"; ids = self.root / "j.ids"
        rng = _random.Random(7)
        seqs = dict(seqs or {})
        for k, v in contigs.items():
            seqs.setdefault(k, "".join(rng.choice("ACGT") for _ in range(v)))
        ref.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        lines = []
        for read, qlen, segs in paf_rows:
            for seg in segs:
                qs, qe, strand, contig, ts, te = seg[:6]
                mapq = seg[6] if len(seg) > 6 else 60
                tp = seg[7] if len(seg) > 7 else "P"
                lines.append("\t".join(map(str, [read, qlen, qs, qe, strand, contig, len(seqs[contig]), ts, te, te - ts, te - ts,
                                                 mapq, f"tp:A:{tp}", f"AS:i:{te - ts}"])))
        paf.write_text("\n".join(lines + list(extra_paf)) + "\n")
        r = subprocess.run([os.sys.executable, str(REPO / "scripts" / "ont_junctions.py"), str(paf), str(ref), str(out),
                            "--events", str(ev), "--support-ids", str(ids), *extra], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        return (list(_csv.DictReader(out.open(), delimiter="\t")), list(_csv.DictReader(ev.open(), delimiter="\t")),
                ids.read_text().split())

    def long_background(self, n=200, qlen=20000):
        # Library reads long enough that a tandem copy of a 2 kb unit would reach outside sequence.
        return [(f"bg{i}", qlen, [(0, qlen, "+", "chr1", 5000, 5000 + qlen)]) for i in range(n)]

    def circle_reads(self, n, prefix="c", jitter=0, both_strands=True):
        # Linearized 2 kb circle [10000,12000): each read <= one unit, crossing the junction once.
        out = []
        for i in range(n):
            d = (i * 7) % (jitter + 1) if jitter else 0
            cut = 10500 + 37 * i
            x, y = 12000 - d, 10000 + d
            la = x - cut; lb = 1900 - la
            segs = [(0, la, "+", "chr1", cut, x), (la, 1900, "+", "chr1", y, y + lb)]
            if both_strands and i % 2:  # the same kind of molecule read from the other strand
                segs = [(0, lb, "-", "chr1", y, y + lb), (lb, 1900, "-", "chr1", cut, x)]
            out.append((f"{prefix}{i}", 1900, segs))
        return out

    def test_ont_junction_circle_reads_are_junction_supported_and_compatible(self):
        rows, events, ids = self.junctions(self.long_background() + self.circle_reads(10))
        self.assertEqual(len(rows), 1); r = rows[0]
        self.assertEqual((r["contig"], r["start"], r["end"], r["location"]), ("chr1", "10000", "12000", "internal"))
        self.assertEqual((r["junction_reads"], r["independent_molecules"], r["anchored_molecules"]), ("10", "10", "10"))
        self.assertEqual((r["label"], r["circle_compatible"], r["tandem_flank_molecules"]), ("junction_supported", "yes", "0"))
        self.assertEqual(r["support_ge3"], "yes")
        self.assertEqual(sorted(ids), sorted(f"c{i}" for i in range(10)))
        self.assertEqual({e["junction_id"] for e in events}, {r["junction_id"]})
        self.assertTrue(all(e["raw_start"] == "10000" and e["raw_end"] == "12000" for e in events))

    def test_ont_junction_reads_into_unique_flank_are_tandem_supported(self):
        reads = self.long_background()
        for i in range(6):  # genome ...L [10000,12000) [10000,12000) R...: reads run from L through the junction
            reads.append((f"t{i}", 5000, [(0, 3000 + 37 * i, "+", "chr1", 9000 - 37 * i, 12000), (3000 + 37 * i, 5000, "+", "chr1", 10000, 12000 - 37 * i)]))
        rows, _, _ = self.junctions(reads)
        r = rows[0]
        self.assertEqual((r["label"], r["circle_compatible"], r["tandem_flank_molecules"]), ("tandem_supported", "no", "6"))

    def test_ont_junction_cigar_insertion_through_tandem_counts_as_tandem(self):
        # Long reads through ...L [10000,12000)[10000,12000) R... aligned as ONE alignment with a 2 kb insertion.
        reads = self.long_background() + self.circle_reads(4)
        paf_extra = []
        for i in range(3):
            ts, te = 8000 - 37 * i, 14000 + 37 * i
            ql = te - ts + 2000
            paf_extra.append("\t".join(map(str, [f"ins{i}", ql, 0, ql, "+", "chr1", 30000, ts, te, te - ts, ql, 60,
                                                  "tp:A:P", f"AS:i:{te - ts}", f"cg:Z:{12000 - ts}M2000I{te - 12000}M"])))
        rows, events, ids = self.junctions(reads, extra_paf=paf_extra)
        r = rows[0]
        self.assertEqual((r["insertion_tandem_molecules"], r["tandem_flank_molecules"], r["label"], r["circle_compatible"]),
                         ("3", "3", "tandem_supported", "no"))
        self.assertEqual((r["split_molecules"], r["independent_molecules"]), ("4", "7"))
        self.assertTrue({"ins0", "ins1", "ins2"} <= set(ids))
        self.assertEqual(sum(e["evidence"] == "insertion" for e in events), 3)

    def test_ont_junction_insertion_is_credited_to_one_nearby_junction_only(self):
        # A second, weaker junction 70 bp away (2 split reads) also fits the 2 kb insertions; they
        # must count once, for the junction with more split-read events.
        reads = self.long_background() + self.circle_reads(4)
        for i in range(2):
            cut = 10600 + 37 * i
            reads.append((f"w{i}", 1800, [(0, 12000 - cut, "+", "chr1", cut, 12000), (12000 - cut, 1800, "+", "chr1", 10070, 10070 + 1800 - (12000 - cut))]))
        paf_extra = []
        for i in range(3):
            ts, te = 8000 - 37 * i, 14000 + 37 * i
            ql = te - ts + 2000
            paf_extra.append("\t".join(map(str, [f"ins{i}", ql, 0, ql, "+", "chr1", 30000, ts, te, te - ts, ql, 60,
                                                  "tp:A:P", f"AS:i:{te - ts}", f"cg:Z:{12000 - ts}M2000I{te - 12000}M"])))
        rows, events, _ = self.junctions(reads, extra_paf=paf_extra)
        by = {(r["start"], r["end"]): r for r in rows}
        self.assertEqual(set(by), {("10000", "12000"), ("10070", "12000")})
        self.assertEqual(by[("10000", "12000")]["insertion_tandem_molecules"], "3")
        self.assertEqual((by[("10070", "12000")]["insertion_tandem_molecules"], by[("10070", "12000")]["independent_molecules"]), ("0", "2"))
        self.assertEqual(sum(e["evidence"] == "insertion" for e in events), 3)

    def test_ont_junction_single_flank_read_is_not_tandem_and_overlap_proves_nothing(self):
        flank = [("x0", 5000, [(0, 3000, "+", "chr1", 9000, 12000), (3000, 5000, "+", "chr1", 10000, 12000)])]
        twice = [("m0", 3800, [(0, 1500, "+", "chr1", 10500, 12000), (1500, 3500, "+", "chr1", 10000, 12000),
                               (3500, 3800, "+", "chr1", 10000, 10300)])]  # passes the unit twice, never leaves it
        rows, _, _ = self.junctions(self.long_background() + self.circle_reads(20) + flank + twice)
        r = rows[0]
        self.assertEqual(r["tandem_flank_molecules"], "1")
        self.assertEqual(r["label"], "junction_supported")           # one flank read: not tandem_supported
        self.assertEqual(r["circle_compatible"], "no")                 # ...but it does exclude "compatible"
        self.assertGreaterEqual(int(r["multipass_molecules"]), 1)     # recorded, not used for topology

    def test_ont_junction_repeat_or_short_anchors_are_ambiguous(self):
        reads = self.long_background()
        for i in range(5):  # anchor a has an equally good secondary alignment elsewhere
            reads.append((f"r{i}", 1900, [(0, 1000, "+", "chr1", 11000 - 37 * i, 12000), (1000, 1900, "+", "chr1", 10000, 10900 - 37 * i),
                                          (0, 1000, "+", "chr1", 20000, 21000 + 37 * i, 0, "S")]))
        rows, _, _ = self.junctions(reads)
        self.assertEqual((rows[0]["label"], rows[0]["repeat_ambiguous_molecules"]), ("ambiguous", "5"))
        short = [(f"s{i}", 600, [(0, 250, "+", "chr1", 11750 - 37 * i, 12000), (250, 600, "+", "chr1", 10000, 10350 - 37 * i)]) for i in range(5)]
        rows, _, _ = self.junctions(self.long_background() + short)
        self.assertEqual(rows[0]["label"], "ambiguous")

    def test_ont_junction_duplex_pair_counts_as_one_molecule(self):
        mol = self.circle_reads(1, both_strands=False)[0]
        (qa0, qa1, _, _, ta0, ta1), (qb0, qb1, _, _, tb0, tb1) = mol[2]
        comp = ("c0_complement", 1900, [(0, qb1 - qb0, "-", "chr1", tb0, tb1), (qb1 - qb0, 1900, "-", "chr1", ta0, ta1)])
        rows, _, _ = self.junctions(self.long_background() + [mol, comp] + self.circle_reads(3, prefix="o")[1:])
        r = rows[0]
        self.assertEqual((r["junction_reads"], r["independent_molecules"]), ("4", "3"))

    def test_ont_junction_microhomology_equivalent_breakpoints_merge(self):
        import random as _random
        rng = _random.Random(3)
        s = [rng.choice("ACGT") for _ in range(30000)]
        s[11880:12000] = s[9880:10000]  # 120 bp identical before Y=10000 and X=12000
        seq = "".join(s)
        reads = self.long_background()
        for i, shift in enumerate([0] * 4 + [100] * 4):  # same junction reported 100 bp apart (> tolerance)
            x, y = 12000 - shift, 10000 - shift
            reads.append((f"h{i}", 1900, [(0, 900, "+", "chr1", x - 900 + 37 * i, x), (900, 1900, "+", "chr1", y, y + 1000 - 37 * i)]))
        rows, events, _ = self.junctions(reads, seqs={"chr1": seq})
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["raw_clusters"].count(":"), 2)          # both raw clusters kept
        self.assertGreaterEqual(int(r["microhomology_bp"]), 100)
        self.assertEqual(r["independent_molecules"], "8")
        self.assertEqual({e["raw_start"] for e in events if e["read"].startswith("h")}, {"10000", "9900"})

    def test_ont_junction_endpoint_jitter_stays_one_junction(self):
        rows, _, _ = self.junctions(self.long_background() + self.circle_reads(30, jitter=40))
        self.assertEqual(len(rows), 1); self.assertEqual(rows[0]["junction_reads"], "30")

    def test_ont_junction_whole_contig_with_assembler_overlap_is_separate_class(self):
        # Circular 3 kb plasmid assembled as 3100 bp: the last 100 bp repeat the first 100 bp.
        reads = self.long_background()
        for i in range(8):
            o = (i * 13) % 101
            reads.append((f"p{i}", 2000, [(0, 1000, "+", "pl", 2000 + o, 3000 + o), (1000, 2000, "+", "pl", o, 1000 + o)]))
        rows, _, _ = self.junctions(reads, contigs={"chr1": 30000, "pl": 3100})
        self.assertEqual([(r["contig"], r["start"], r["end"], r["location"]) for r in rows], [("pl", "0", "3100", "whole_contig")])

    def test_ont_junction_contig_edge_flagged_not_dropped(self):
        reads = [(f"e{i}", 1500, [(0, 700, "+", "chr1", 29300 - 37 * i, 30000 - 40), (700, 1500, "+", "chr1", 25000, 25800 - 37 * i)]) for i in range(4)]
        rows, _, _ = self.junctions(self.long_background() + reads)
        self.assertEqual((rows[0]["location"], rows[0]["label"]), ("contig_edge", "junction_supported"))

    def test_ont_junction_short_library_is_not_called_compatible(self):
        reads = [(f"bg{i}", 1500, [(0, 1500, "+", "chr1", 5000, 6500)]) for i in range(200)]
        reads += [(f"c{i}", 1500, [(0, 800, "+", "chr1", 11200 - 9 * i, 12000), (800, 1500, "+", "chr1", 10000, 10700 - 9 * i)]) for i in range(5)]
        rows, _, _ = self.junctions(reads)
        self.assertEqual((rows[0]["label"], rows[0]["circle_compatible"]), ("junction_supported", "na"))

    def test_ont_junction_ignores_indels_and_other_contigs(self):
        reads = [("d", 3000, [(0, 1000, "+", "chr1", 1000, 2000), (1000, 3000, "+", "chr1", 2500, 4500)]),  # deletion: jump forward
                 ("o", 3000, [(0, 1000, "+", "chr1", 1000, 2000), (1000, 3000, "+", "pl", 0, 2000)]),       # translocation
                 ("s", 3000, [(0, 1000, "+", "chr1", 1000, 2000), (1000, 3000, "-", "chr1", 500, 1500)]),   # inversion
                 ("i", 3000, [(0, 1000, "+", "chr1", 1000, 2000), (1000, 3000, "+", "chr1", 1960, 3960)])]  # small overlap
        rows, _, ids = self.junctions(reads)
        self.assertEqual(rows, []); self.assertEqual(ids, [])

    def test_library_depth_counts_unmapped_reads_from_input(self):
        fai = self.root / "ref.fa.fai"; fai.write_text("chr\t200000\t5\t60\t61\npl\t10000\t9\t60\t61\n")
        paf = self.root / "lib.paf"
        paf.write_text("r1\t1000\t0\t1000\t+\tchr\t200000\t0\t1000\t1000\t1000\t60\ttp:A:P\n"
                       "r2\t3000\t0\t3000\t+\tpl\t10000\t0\t3000\t3000\t3000\t60\ttp:A:P\n")
        fa = self.root / "reads.fa"; fa.write_text(">r1\n" + "A" * 1000 + "\n>r2\n" + "A" * 3000 + "\n>u1\n" + "A" * 2000 + "\n>u2\n" + "A" * 2000 + "\n")
        j = library_depth(paf, fai, self.root / "depth.tsv", fa)
        self.assertEqual((j["reads_input"], j["reads_in_paf"], j["reads_mapped_primary"]), (4, 2, 2))
        self.assertEqual((j["mapped_fraction"], j["aligned_base_fraction"], j["read_n50"]), (0.5, 0.5, 2000))
        self.assertIn("pl\t10000\t3000\t0.30", (self.root / "depth.tsv").read_text())

    def ont_rule_context(self, cap=2):
        x = SimpleNamespace(s={"platform": "nanopore", "reference": str(self.ref)}, spec={"reference_sha256": "0" * 64},
                            c={"caller": {"cresil_mode": "linkage", "min_region_size": 200, "breakpoint_depth": 3,
                                          "average_depth": 5, "cresil_graph_read_cap": cap}},
                            cpus=4, ref=self.ref, state=self.root / "st", work=self.root / "w")
        x.out, x.log = x.work / "caller", x.state / "logs"
        for d in (x.out, x.log):
            d.mkdir(parents=True, exist_ok=True)
        (self.root / "ref.fa.fai").write_text("chr1\t400\t6\t400\t401\n")
        return x

    def test_ont_call_keeps_junction_result_when_cresil_fails(self):
        from circlestream.rules import call
        x = self.ont_rule_context()
        (x.out / "eccDNA_final.txt").write_text("stale")
        def run(args, **kw):
            if any("cresil_wgls_checked.py" in str(a) for a in args):
                (x.log / "call.1.1.stderr.log").write_text("worker OOM\nAttributeError: 'DataFrame' object has no attribute 'ovl_5end'\n")
                raise subprocess.CalledProcessError(1, args)
            raise AssertionError(f"nothing may run after a CReSIL failure: {args}")
        x.run = run
        with patch("circlestream.rules.ont_junctions", return_value={"junctions": 7}):
            call(x)
        j = read_json(x.state / "calls.json")
        self.assertEqual(j["ont_junctions"], {"junctions": 7})
        self.assertEqual((j["cresil"]["status"], j["candidate_count"]), ("failed", None))
        self.assertIn("ovl_5end", j["cresil"]["error_tail"][-1])
        self.assertFalse((x.out / "eccDNA_final.txt").exists())

    def ont_evidence_fixture(self, x, bam_records=None):
        reads = ["j1", "g1", "g2", "g3", "g4", "other"]
        x.inputs = [self.root / "reads.fastq"]; fq(x.inputs[0], reads)
        (x.out / "eccDNA_final.txt").write_text("id\tmerge_region\tmerge_len\tnum_region\tctc\tnumreads\n"
                                                "ec1\tchr1:1-300_+\t300\t1\tTrue\t5\nec2\tchr1:1-50_+\t50\t1\tFalse\t1\n")
        for gid, members in (("ec1", ["j1", "g1", "g2", "g3", "g4"]), ("ec2", ["other"])):
            d = x.out / "cresil_run" / "assemGraph" / gid; d.mkdir(parents=True)
            (d / "final_reads.fa").write_text("".join(f">{gid}_{r}_0_100\nACGT\n" for r in members))
        (x.out / "ont_junction_ids.txt").write_text("j1\n")
        (x.out / "trim.txt").write_text("readid\tq_start\nj1\t0\nother\t0\n")
        def shell(code, **kw):
            (x.work / "evidence" / "support.bam").write_bytes(b"bam")
        def run(args, stdout=None, **kw):
            if args[:3] == ["samtools", "view", "-c"]:
                with gzip.open(x.work / "evidence" / "support_reads1.fastq.gz", "rt") as h:
                    n = sum(1 for i, _ in enumerate(h) if i % 4 == 0)
                Path(stdout).write_text(f"{n if bam_records is None else bam_records}\n")
        x.shell, x.run = shell, run

    def test_ont_evidence_caps_cyclic_graph_reads_and_keeps_bam_not_fastq(self):
        from circlestream.rules import evidence
        x = self.ont_rule_context(cap=2); self.ont_evidence_fixture(x)
        evidence(x)
        ev = x.work / "evidence"
        ids = (ev / "read_ids.txt").read_text().split()
        self.assertIn("j1", ids)                     # junction reads always kept (and ranked first in the graph)
        self.assertEqual(len([r for r in ids if r.startswith("g")]), 1)  # cap 2 = j1 + one graph read
        self.assertNotIn("other", ids)               # non-cyclic graph: IDs only
        with gzip.open(ev / "cresil" / "graph_reads.tsv.gz", "rt") as h:
            self.assertEqual(len(h.read().splitlines()), 1 + 6)  # every graph read ID retained
        self.assertEqual(read_json(ev / "cresil" / "graph_read_sequences.json")["graphs"]["ec1"], {"reads": 5, "sequences_kept": 2})
        self.assertFalse((ev / "support_reads1.fastq.gz").exists())
        self.assertTrue((ev / "support.bam").is_file())
        self.assertIn("samtools fastq -F 0x900", read_json(x.state / "evidence.json")["support_reads"])

    def test_ont_evidence_refuses_to_drop_fastq_if_bam_is_incomplete(self):
        from circlestream.rules import evidence
        x = self.ont_rule_context(cap=2); self.ont_evidence_fixture(x, bam_records=1)
        with self.assertRaisesRegex(ValueError, "primary records"):
            evidence(x)

    def test_ont_evidence_without_cresil_output_keeps_junction_reads(self):
        from circlestream.rules import evidence, copy_caller_evidence
        x = self.ont_rule_context(); self.ont_evidence_fixture(x)
        shutil.rmtree(x.out / "cresil_run"); (x.out / "eccDNA_final.txt").unlink(); (x.out / "trim.txt").unlink()
        evidence(x)
        self.assertEqual((x.work / "evidence" / "read_ids.txt").read_text().split(), ["j1"])
        self.assertEqual(read_json(x.work / "evidence" / "caller-retention.json")["support_trim_rows"], 0)

    def test_config_validates_graph_cap_and_release_flag(self):
        for bad in ({"caller": dict(self.c["caller"], cresil_graph_read_cap=-1)}, {"release_failed_raw": "yes"},
                    {"release_failed_raw": True, "download_mode": "dtn"}):
            p = self.root / "cfg.json"; write_json(p, dict(self.c, **bad))
            with self.assertRaises(ValueError):
                load_config(p)

    def test_cresil_trace_requires_spanning_entry_not_sub_unit_overlap(self):
        from scripts.cresil_trace import final_spans, span_matches
        p = self.root / "eccDNA_final.txt"
        p.write_text("id\tmerge_region\tmerge_len\tnum_region\tctc\n"
                     "ec1\tpl:7-996_-\t990\t1\tTrue\n"                                   # cyclic sub-unit of a 3066 bp plasmid
                     "ec2\tchr:2882024-2885212_+,chr:2885213-2900183_+\t18160\t2\tFalse\n"  # one unit split in two regions
                     "ec3\tpl2:3001-5000_+,pl2:1-3000_+\t5000\t2\tTrue\n")                 # closure written as two arcs
        spans = final_spans(p, 100)
        self.assertIn(("chr", 2882023, 2900183, "ec2", False), spans)
        self.assertIn(("pl2", 0, 5000, "ec3", True), spans)
        hit = lambda c, s, e: [x[3] for x in spans if span_matches(c, s, e, 0, x)]  # noqa: E731
        self.assertEqual(hit("pl", 0, 3066), [])             # overlap with a sub-unit is not agreement
        self.assertEqual(hit("chr", 2882000, 2900186), ["ec2"])
        self.assertEqual(hit("chr", 2850683, 2900186), [])   # same end, different start
        self.assertEqual(hit("pl2", 0, 5000), ["ec3"])

    def test_config_rejects_unknown_ont_junction_key(self):
        c = dict(self.c, caller=dict(self.c["caller"], ont_junction={"min_reads": 3, "bogus": 1}))
        p = self.root / "cfg.json"; write_json(p, c)
        with self.assertRaisesRegex(ValueError, "ont_junction"):
            load_config(p)


if __name__=="__main__": unittest.main()
