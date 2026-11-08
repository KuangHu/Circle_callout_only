import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time

from .common import GB, REPO, digest, fingerprint, is_sra, read_json, source_ok, state_paths, work_root, write_json
from .fastq import extract, materialize
from .sra import network_isolation, refseq_deps, run_limited, sra_tool, unresolved

CF_SHA256 = "a7df874bb976cd8761840f91190df43ac3dc949b7e97796441b86b0bb699acd2"
CRESIL_COMMIT = "78356affa6adca6a41ea297c3d7d93f5ecb5ae6a"
# Upstream read-ID filters, e.g. grep -w -Ff $6-$7\.split.id-freq2.txt $6-$7\.split.txt > $6-$7\.split_freq2.txt
WORD_GREP = re.compile(r"grep -w -Ff (\S+) (\S+) > (\S+)")


class Context:
    def __init__(self, spec):
        self.spec = spec
        self.c, self.s = spec["config"], spec["sample"]
        self.state, self.scratch = state_paths(self.c, self.s["sample_id"])
        self.root = work_root(self.c, self.s["sample_id"])
        self.work = self.root / "work"
        self.raw_files = [self.scratch / "raw" / f["name"] for f in self.s["files"]]
        self.sra = any(is_sra(f) for f in self.s["files"])
        if self.sra:
            # FASTQ converted from the SRA object during qc (fasterq-dump --split-3 naming).
            run = self.raw_files[0].stem
            fq_dir = self.root / "fastq"
            self.inputs = ([fq_dir / f"{run}_1.fastq", fq_dir / f"{run}_2.fastq"]
                           if self.s["platform"] == "illumina" else [fq_dir / f"{run}.fastq"])
        else:
            self.inputs = self.raw_files
        self.prefix = self.c["tools"][self.s["platform"] + "_prefix"]
        self.cpus = self.c["slurm"]["cpus"]
        self.ref = self.work / "reference" / "genome.fa"
        self.out = self.work / "caller"
        self.log = self.state / "logs"
        self.log.mkdir(exist_ok=True)
        self.counter = 0

    def run(self, args, *, cwd=None, stdout=None):
        self.counter += 1
        cmd = self.prefix + list(map(str, args))
        name = f"{os.environ.get('CS_RULE','step')}.{time.time_ns()}.{self.counter}"
        write_json(self.log / (name + ".command.json"), {"argv": cmd, "cwd": str(cwd or self.work)})
        with (self.log / (name + ".stderr.log")).open("w") as err:
            if stdout is not None:
                with Path(stdout).open("w") as out:
                    subprocess.run(cmd, cwd=cwd or self.work, stdout=out, stderr=err, check=True)
            else:
                with (self.log / (name + ".stdout.log")).open("w") as out:
                    subprocess.run(cmd, cwd=cwd or self.work, stdout=out, stderr=err, check=True)

    def shell(self, code, *, cwd=None, args=()):
        self.run(["bash", "-euo", "pipefail", "-c", code, "circle-stream", *args], cwd=cwd)


def preflight(x):
    if read_json(x.scratch / ".owner.json") != {"spec_sha256": fingerprint(x.spec)}:
        raise ValueError("Owned scratch/spec mismatch")
    if digest(x.s["reference"]) != x.spec["reference_sha256"]:
        raise ValueError("Reference changed since init")
    # Upstream CReSIL constructs shell commands without escaping; constrain its paths.
    for p in (x.work, x.state, Path(x.s["reference"])):
        if not re.fullmatch(r"[A-Za-z0-9_./-]+", str(p)):
            raise ValueError(f"Caller paths must not contain spaces/shell metacharacters: {p}")
    for path, item in zip(x.raw_files, x.s["files"]):
        if not source_ok(path, item):
            raise ValueError(f"Input checksum mismatch: {path}")
    if x.sra:
        # Internal per-column MD5s of the SRA object, in addition to the S3 transfer checksum.
        result = subprocess.run([sra_tool("vdb-validate"), str(x.raw_files[0])], capture_output=True, text=True,
                                cwd=x.raw_files[0].parent)
        (x.state / "vdb-validate.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode != 0 or "is consistent" not in result.stdout + result.stderr:
            raise ValueError(f"vdb-validate failed for {x.raw_files[0]}")
        # Aligned-read runs: every reference sequence must already sit beside the object
        # (fetched by the download step), so conversion never touches the network.
        deps = refseq_deps(x.raw_files[0])
        write_json(x.state / "sra-refseqs.json", [{"seqid": i, "location": loc} for i, loc in deps])
        if unresolved(deps):
            raise ValueError(f"SRA object needs reference sequences that were not downloaded: {unresolved(deps)}")
    binaries = ["samtools", "bedtools", "bash"]
    binaries += ["bwa", "samblaster", "awk", "grep", "sed", "sort"] if x.s["platform"] == "illumina" else ["cresil", "minimap2", "mosdepth"]
    x.shell("command -v " + " ".join(binaries))
    x.run(["samtools", "--version"], stdout=x.state / "samtools-version.txt")
    env_code = "import sys,pathlib,json; root=pathlib.Path(sys.prefix)/'conda-meta'; rows=[json.loads(p.read_text()) for p in sorted(root.glob('*.json'))]; print(json.dumps({'python':sys.version,'prefix':sys.prefix,'packages':[{k:r.get(k) for k in ('name','version','build','subdir')} for r in rows]},indent=2))"
    x.run(["python", "-c", env_code], stdout=x.state / "environment-packages.json")
    if x.s["platform"] == "nanopore":
        verify_code = "import cresil,pathlib,hashlib; root=pathlib.Path(cresil.__file__).parent; expected={'identify_wgls.py':'093c6c92babf4decc59d3c888ecbc5234fa19a2a','trim.py':'cb50b2f257ad0365660cbb15fdbf5e8ffcafb618'}; actual={n:hashlib.sha1(b'blob '+str(len((root/'cli'/n).read_bytes())).encode()+b'\\0'+(root/'cli'/n).read_bytes()).hexdigest() for n in expected}; print(actual); assert actual==expected, 'CReSIL source differs from pinned revision'"
        x.run(["python", "-c", verify_code], stdout=x.state / "cresil-source-check.txt")
        x.run(["cresil", "identify_wgls", "--help"], stdout=x.state / "cresil-help.txt")
        x.run(["python", "-c", "import cresil,inspect; print(inspect.getfile(cresil))"], stdout=x.state / "cresil-install.txt")
        x.run(["python", "-m", "pip", "freeze"], stdout=x.state / "python-packages.txt")
    else:
        source = (REPO / "third_party" / "Circle_finder" / "original.sh").read_text()
        if hashlib.sha256(source.encode()).hexdigest() != CF_SHA256:
            raise ValueError("Vendored Circle_finder source differs from reviewed version")
    write_json(x.state / "provenance.json", {
        "spec": x.spec, "reference_sha256": digest(x.s["reference"]),
        "circle_finder_commit": "3eb333db2ea6277dde36cbf640be9afeb710c717",
        "cresil_expected_commit": CRESIL_COMMIT,
        "pipeline_sources": {str(p.relative_to(REPO)): digest(p) for p in sorted(REPO.rglob("*.py")) if "__pycache__" not in p.parts},
        "started": time.time(), "slurm_job_id": os.environ.get("SLURM_JOB_ID")
    })


def convert_sra(x):
    """fasterq-dump into raw/fastq; random-access temp files go to node-local disk."""
    out = x.inputs[0].parent
    if all(p.is_file() for p in x.inputs) and (out / ".complete").is_file():
        return read_json(out / ".complete")
    if out.exists():
        shutil.rmtree(out)  # partial conversion from an interrupted job; source object is retained
    out.mkdir(parents=True)
    local = os.environ.get("TMPDIR") if os.environ.get("TMPDIR", "").startswith("/tmp") else "/tmp"
    start = time.time()
    limit = 60 * x.c.get("sra_conversion_timeout_minutes", 60)
    isolate = network_isolation()  # refseqs are local (checked in preflight); fail fast instead of hanging
    with tempfile.TemporaryDirectory(prefix="cs-fqd-", dir=local) as tmp:
        rc, so, se, timed_out = run_limited(isolate + [sra_tool("fasterq-dump"), "--split-3", "--skip-technical",
                                                       "-e", str(x.cpus), "-t", tmp, "-O", str(out), str(x.raw_files[0])], limit,
                                             cwd=x.raw_files[0].parent)  # refseqs resolve from the working directory
    (x.state / "fasterq-dump.log").write_text(so + se, encoding="utf-8")
    if timed_out:
        raise RuntimeError(f"fasterq-dump exceeded sra_conversion_timeout_minutes ({limit // 60}); see fasterq-dump.log")
    if rc != 0:
        raise RuntimeError(f"fasterq-dump failed ({rc}); see fasterq-dump.log")
    produced = sorted(p.name for p in out.glob("*.fastq"))
    expected = sorted(p.name for p in x.inputs)
    orphans = [n for n in produced if n not in expected]
    if [n for n in expected if n not in produced]:
        raise ValueError(f"SRA layout does not match platform {x.s['platform']}: produced {produced}, expected {expected}")
    if x.s["platform"] == "nanopore" and orphans:
        raise ValueError(f"Unexpected extra reads from ONT SRA object: {orphans}")
    # Illumina --split-3 puts mate-less reads in <run>.fastq; record and drop them (pairs only).
    info = {"seconds": time.time() - start, "tool": "fasterq-dump --split-3 --skip-technical",
            "network_isolated": bool(isolate), "timeout_seconds": limit,
            "outputs": {n: (out / n).stat().st_size for n in produced}, "dropped_unpaired_files": orphans}
    for n in orphans:
        (out / n).unlink()
    write_json(out / ".complete", info)
    return info


def qc(x):
    conversion = convert_sra(x) if x.sra else None
    stats = materialize(x.inputs, x.work / "prepared", x.s["platform"],
                        x.c["caller"].get("short_read_policy", "reject"))
    if conversion:
        stats["sra_conversion"] = conversion
    write_json(x.state / "qc.json", stats)


def reference(x):
    x.ref.parent.mkdir(exist_ok=True)
    shutil.copyfile(x.s["reference"], x.ref)
    if digest(x.ref) != x.spec["reference_sha256"]:
        raise ValueError("Reference changed while preparing indexes")
    seen = set()
    with x.ref.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in seen or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                    raise ValueError("Reference contig IDs must be unique and use only letters, digits, _, ., -")
                seen.add(name)
    if not seen:
        raise ValueError("Reference has no FASTA headers")
    x.run(["samtools", "faidx", x.ref])
    if x.s["platform"] == "illumina":
        x.run(["bwa", "index", x.ref])
    else:
        x.run(["minimap2", "-x", "map-ont", "-t", x.cpus, "-d", str(x.ref)+".mmi", x.ref])
        x.run(["samtools", "faidx", x.work / "prepared" / "reads.fa"])


def cf_lines(stage):
    source = (REPO / "third_party" / "Circle_finder" / "original.sh").read_text()
    if hashlib.sha256(source.encode()).hexdigest() != CF_SHA256:
        raise ValueError("Unexpected Circle_finder source")
    current = 0
    selected = []
    for line in source.splitlines():
        match = re.match(r"#\s*step\s*(\d+)", line, re.I)
        if match:
            current = int(match.group(1))
        if not line.strip() or line.startswith("#"):
            continue
        if line == "rm *hg38.sam *hg38.bam":
            continue
        if (stage == "mapping" and current <= 2) or (stage == "call" and current >= 3):
            selected.append(line)
    return selected


def circle_finder(x, stage):
    # Preserve scientific shell expressions, change execution/error/cleanup handling only.
    preamble = 'grep() { local rc=0; command grep "$@" || rc=$?; if [ "$rc" -gt 1 ]; then return "$rc"; fi; return 0; }; export LC_ALL=C; '
    args = [str(max(1, x.cpus-2)), str(x.ref), str(x.work/"prepared"/"reads1.fastq"),
            str(x.work/"prepared"/"reads2.fastq"), str(x.c["caller"]["circle_finder_min_non_overlap"]), "sample", "ref"]
    for index, code in enumerate(cf_lines(stage), 1):
        if code == "samtools view -@ $1 -bS $6-$7\\.unmap.sam > $6-$7\\.unmap.bam" and (x.out / "sample-ref.unmap.sam").stat().st_size == 0:
            # samblaster writes no header when nothing is unmapped/clipped; unmap.bam is unused downstream.
            write_json(x.log / f"cf-{stage}-{index}.json", {"command": code, "skipped": "empty unmap.sam"})
            continue
        record = {"command": code}
        m = WORD_GREP.fullmatch(code)
        if m:
            # Same selection as upstream `grep -w -Ff`, without GNU grep's slow path (scripts/wordgrep.py).
            record["executed"] = code = (f"{shlex.quote(sys.executable)} {shlex.quote(str(REPO / 'scripts' / 'wordgrep.py'))} "
                                         f"{m.group(1)} {m.group(2)} > {m.group(3)}")
        started = time.monotonic()
        x.shell(preamble + code, cwd=x.out, args=args)
        record["seconds"] = time.monotonic() - started
        write_json(x.log / f"cf-{stage}-{index}.json", record)


def mapping(x):
    if x.out.exists():
        # Never mix a retried caller with partial output from an earlier attempt.
        x.out.rename(x.work / ("caller.failed." + str(time.time_ns())))
    x.out.mkdir()
    if x.s["platform"] == "illumina":
        circle_finder(x, "mapping")
        x.run(["samtools", "quickcheck", "-v", x.out / "sample-ref.sorted.bam"])
    else:
        x.run(["cresil", "trim", "-t", max(4, x.cpus-2), "-fq", x.work/"prepared"/"reads.fa",
               "-r", str(x.ref)+".mmi", "-o", x.out])
        with (x.out / "trim.txt").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"readid", "q_start", "q_end", "ref", "r_start", "r_end"} <= set(reader.fieldnames or []):
                raise ValueError("CReSIL trim schema mismatch")


def call(x):
    if x.s["platform"] == "illumina":
        circle_finder(x, "call")
        result = x.out / "sample-ref.microDNA-JT.txt"
        rows = []
        with result.open() as handle:
            for line in handle:
                a = line.rstrip().split("\t")
                if len(a) != 4 or not (0 <= int(a[1]) < int(a[2])) or int(a[3]) <= 0:
                    raise ValueError("Circle_finder result schema mismatch")
                rows.append(a)
        write_json(x.state / "calls.json", {"candidate_count": len(rows), "classification": "eccDNA_or_tandem_duplication", "caller": "Circle_finder"})
    else:
        # The upstream caller appends some intermediates; remove neither nor reuse these on retry.
        nested = x.out / "cresil_run"
        if nested.exists():
            nested.rename(x.work / ("cresil_run.failed." + str(time.time_ns())))
        (x.out / "eccDNA_final.txt").unlink(missing_ok=True)
        p = x.c["caller"]
        x.run(["python", REPO/"scripts"/"cresil_wgls_checked.py", "-t", max(4,x.cpus-2), "-m", p["cresil_mode"],
               "-r", str(x.ref)+".mmi", "-fa", x.ref, "-fai", str(x.ref)+".fai",
               "-fq", x.work/"prepared"/"reads.fa", "-trim", x.out/"trim.txt",
               "-minrsize", p["min_region_size"], "-break", p["breakpoint_depth"], "-depth", p["average_depth"]])
        with (x.out / "eccDNA_final.txt").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"id", "ctc", "merge_region", "numreads"} <= set(reader.fieldnames or []):
                raise ValueError("CReSIL result schema mismatch")
            rows = list(reader)
        if any(r["ctc"] not in ("True", "False") for r in rows):
            raise ValueError("Unexpected CReSIL cyclic flag")
        write_json(x.state / "calls.json", {"candidate_count": sum(r["ctc"] == "True" for r in rows),
            "all_graph_count": len(rows), "classification": "cyclic_graph_candidate_not_topology_proof", "caller": "CReSIL"})


def copy_caller_evidence(out, ev, platform, wanted):
    """Explicit allowlist: never archive whole-library/depth intermediates."""
    selected = []
    omitted_bytes = 0
    for p in out.rglob("*"):
        if p.is_symlink():
            raise ValueError("Unexpected symlink in caller output")
        if not p.is_file():
            continue
        rel = p.relative_to(out)
        if platform == "illumina":
            keep = len(rel.parts) == 1 and p.suffix == ".txt" and (
                p.name == "sample-ref.microDNA-JT.txt" or
                p.name.startswith(("sample-ref.split.", "sample-ref.split_", "sample-ref.disc.")))
            dest = ev / rel
        else:
            keep = rel.as_posix() in {
                "eccDNA_final.txt", "no_candidates.json", "cresil_run/subGraphs.summary.txt"}
            keep = keep or (len(rel.parts) == 4 and rel.parts[:2] == ("cresil_run", "assemGraph")
                            and p.name in {"final_reads.fa", "reference_regions.fa"})
            dest = ev / "cresil" / rel
        if keep:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest)
            selected.append({"path": rel.as_posix(), "bytes": p.stat().st_size})
        else:
            omitted_bytes += p.stat().st_size
    trim_rows = 0
    if platform == "nanopore":
        dest = ev / "cresil" / "trim.support.txt"
        dest.parent.mkdir(parents=True, exist_ok=True)
        with (out / "trim.txt").open() as source, dest.open("w", newline="") as target:
            reader = csv.DictReader(source, delimiter="\t")
            if "readid" not in (reader.fieldnames or []):
                raise ValueError("CReSIL trim schema mismatch")
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["readid"] in wanted:
                    writer.writerow(row)
                    trim_rows += 1
    write_json(ev / "caller-retention.json", {"copied_files": selected,
        "caller_bytes_not_copied_verbatim": omitted_bytes, "support_trim_rows": trim_rows,
        "policy": "split/discordant tables and final calls; CReSIL graph sequences/tables and support-only trim rows"})


def evidence(x):
    ev = x.work / "evidence"
    if ev.exists():
        ev.rename(x.work / ("evidence.failed." + str(time.time_ns())))
    ev.mkdir()
    wanted = set()
    if x.s["platform"] == "illumina":
        # Conservative union: retain all split/discordant reads, not only final positives.
        for kind in ("split", "disc"):
            sam = ev / f"{kind}.sam"
            x.run(["samtools", "view", x.out/f"sample-ref.{kind}.bam"], stdout=sam)
            with sam.open() as handle:
                for line in handle:
                    name = line.split("\t",1)[0]
                    # samblaster's splitter file appends _1/_2 (read of pair) to IDs.
                    if kind == "split" and name.endswith(("_1", "_2")):
                        name = name[:-2]
                    wanted.add(name)
            sam.unlink()
        (ev / "read_ids.txt").write_text("".join(n+"\n" for n in sorted(wanted)))
        x.run(["samtools", "view", "-b", "-N", ev/"read_ids.txt", "-o", ev/"support.bam", x.out/"sample-ref.sorted.bam"])
    else:
        for p in x.out.glob("cresil_run/assemGraph/*/final_reads.fa"):
            with p.open() as handle:
                for line in handle:
                    if line.startswith(">"):
                        # Upstream format: graphID_originalReadID_queryStart_queryEnd.
                        header = line[1:].split()[0]
                        wanted.add(header.split("_",1)[1].rsplit("_",2)[0])
        (ev / "read_ids.txt").write_text("".join(n+"\n" for n in sorted(wanted)))
    copy_caller_evidence(x.out, ev, x.s["platform"], wanted)
    counts = extract(x.inputs, wanted, ev, x.s["platform"])
    if x.s["platform"] == "nanopore":
        # Preserve original read alignments with CIGAR/SA/secondary records for downstream review.
        q = shlex.quote
        if wanted:
            x.shell(f"minimap2 -ax map-ont --secondary=yes -t {max(1,x.cpus-2)} {q(str(x.ref))} {q(str(ev/'support_reads1.fastq.gz'))} | samtools sort -@ 1 -o {q(str(ev/'support.bam'))}")
        else:
            (ev / "empty.sam").write_text("@HD\tVN:1.6\tSO:coordinate\n" + "".join(
                f"@SQ\tSN:{a[0]}\tLN:{a[1]}\n" for a in (l.split("\t") for l in Path(str(x.ref)+".fai").read_text().splitlines())))
            x.run(["samtools", "view", "-b", "-o", ev/"support.bam", ev/"empty.sam"])
            (ev / "empty.sam").unlink()
    x.run(["samtools", "quickcheck", "-v", ev/"support.bam"])
    x.run(["samtools", "index", ev/"support.bam"])
    # Reference by identity, not by copy: one genome copy per sample dominated archive size.
    # support.bam can be viewed against the original file, verified by its SHA-256.
    shutil.copy2(str(x.ref)+".fai", ev / "reference.fa.fai")
    write_json(ev / "reference.json", {"path": x.s["reference"], "sha256": x.spec["reference_sha256"],
                                       "store_path": x.spec.get("reference_store_path"),
                                       "note": "Archives do not embed the genome; keep the reference store (or this file)"})
    # Freeze the actual executable workflow, not just a citation to a changing repository.
    code_dest = ev / "pipeline_source"
    for p in list(REPO.glob("*.py")) + list((REPO/"circlestream").glob("*.py")) + list((REPO/"scripts").glob("*.py")) + [REPO/"workflow"/"Snakefile", REPO/"third_party"/"Circle_finder"/"original.sh", REPO/"third_party"/"Circle_finder"/"LICENSE"]:
        dest = code_dest / p.relative_to(REPO)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    write_json(x.state / "evidence.json", {"read_ids": len(wanted), "records_per_file": counts,
        "bytes_before_archive": sum(p.stat().st_size for p in ev.rglob("*") if p.is_file()),
        "scope": "all split/discordant read IDs" if x.s["platform"]=="illumina" else "reads in all reported CReSIL graphs, including acyclic",
        "limitation": "Original full-library alternatives require re-download; no physical-topology validation is claimed"})


def archive(x):
    ev = x.work / "evidence"
    # Include all completed per-rule measurements and command logs available before archiving.
    paths = {"evidence/"+str(p.relative_to(ev)).replace(os.sep,"/"): p for p in ev.rglob("*") if p.is_file()}
    for p in x.state.rglob("*"):
        if p.is_file() and not any(part.startswith(".") for part in p.relative_to(x.state).parts) and not p.name.startswith(("archive.", "slurm-")):
            paths["metadata/"+str(p.relative_to(x.state)).replace(os.sep,"/")] = p
    required = [x.state/"calls.json", x.state/"evidence.json", ev/"support.bam", ev/"reference.json"]
    if any(not p.is_file() for p in required):
        raise ValueError("Missing required evidence")
    if any(p.is_symlink() for p in paths.values()):
        raise ValueError("Archive refuses symlinks")
    total = sum(p.stat().st_size for p in paths.values())
    root = Path(x.c["archive_dir"])
    if shutil.disk_usage(root).free < total + x.c["min_free_gb"] * GB:
        raise ValueError("Insufficient archive capacity; scratch will be retained")
    manifest = {name: digest(p) for name,p in paths.items()}
    final = root / (x.s["sample_id"] + ".tar.gz")
    if final.exists():
        raise ValueError("Archive already exists; verify/reconcile before overwriting")
    temp = root / (x.s["sample_id"] + ".partial.tar.gz")
    with tarfile.open(temp, "w:gz", compresslevel=6) as bundle:
        payload = json.dumps(manifest, sort_keys=True).encode()
        info = tarfile.TarInfo("MANIFEST.json")
        info.size = len(payload)
        bundle.addfile(info, io.BytesIO(payload))
        for name, p in paths.items():
            bundle.add(p, arcname=name, recursive=False)
    with temp.open("r+b") as handle:
        os.fsync(handle.fileno())
    receipt = {"path": str(temp), "sha256": digest(temp), "bytes": temp.stat().st_size}
    write_json(x.state / "archive.json", receipt)
    from .controller import verify_bundle
    verify_bundle(x.state)
    os.replace(temp, final)
    receipt["path"] = str(final)
    write_json(x.state / "archive.json", receipt)
    verify_bundle(x.state)


RULES = {f.__name__: f for f in (preflight,qc,reference,mapping,call,evidence,archive)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--rule", choices=RULES, required=True)
    a = ap.parse_args()
    os.environ["CS_RULE"] = a.rule
    x = Context(read_json(a.spec))
    os.chdir(x.work)
    RULES[a.rule](x)
