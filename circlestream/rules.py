import argparse
import collections
import csv
import gzip
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


CRESIL_FAILED = "cresil_failed.json"


def stderr_tail(x, rule, lines=3):
    log = sorted(x.log.glob(f"{rule}.*.stderr.log"), key=lambda q: q.stat().st_mtime)
    return log[-1].read_text(errors="replace").strip().splitlines()[-lines:] if log else []


def cresil_failed(x, junctions, failure):
    """Archive the junction-counter result with the CReSIL failure recorded (no partial CReSIL output)."""
    for name in ("eccDNA_final.txt", "cresil_trace.tsv"):
        (x.out / name).unlink(missing_ok=True)
    write_json(x.state / "calls.json", {"caller": "ont_junctions+CReSIL",
        "classification": "evidence categories, not topology proof",
        "ont_junctions": junctions,
        "cresil": {"mode": x.c["caller"]["cresil_mode"], "status": "failed", **failure},
        "cresil_trace": {}, "candidate_count": None})


def mapping(x):
    if x.out.exists():
        # Never mix a retried caller with partial output from an earlier attempt.
        x.out.rename(x.work / ("caller.failed." + str(time.time_ns())))
    x.out.mkdir()
    if x.s["platform"] == "illumina":
        circle_finder(x, "mapping")
        x.run(["samtools", "quickcheck", "-v", x.out / "sample-ref.sorted.bam"])
    else:
        # trim is a CReSIL-only prerequisite. If it fails (e.g. OOM-killed), record the failure and
        # let call() still run the independent junction scan instead of failing the sample.
        try:
            x.run(["cresil", "trim", "-t", max(4, x.cpus-2), "-fq", x.work/"prepared"/"reads.fa",
                   "-r", str(x.ref)+".mmi", "-o", x.out])
        except subprocess.CalledProcessError as error:
            (x.out / "trim.txt").unlink(missing_ok=True)
            write_json(x.out / CRESIL_FAILED, {"stage": "trim", "returncode": error.returncode,
                                               "error_tail": stderr_tail(x, "mapping")})
            return
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
        junctions = ont_junctions(x)  # genome-wide scan first: independent of CReSIL's outcome
        # CReSIL is the secondary, comparison caller. If it fails (e.g. a worker OOM-killed in a
        # deep library, seen in 3 of 1,440 samples at 24 GB), the junction-counter result is
        # still archived and the failure is recorded instead of discarding the whole sample.
        if (x.out / CRESIL_FAILED).is_file():  # trim already failed in mapping()
            cresil_failed(x, junctions, json.loads((x.out / CRESIL_FAILED).read_text()))
            return
        try:
            x.run(["python", REPO/"scripts"/"cresil_wgls_checked.py", "-t", max(4,x.cpus-2), "-m", p["cresil_mode"],
                   "-r", str(x.ref)+".mmi", "-fa", x.ref, "-fai", str(x.ref)+".fai",
                   "-fq", x.work/"prepared"/"reads.fa", "-trim", x.out/"trim.txt",
                   "-minrsize", p["min_region_size"], "-break", p["breakpoint_depth"], "-depth", p["average_depth"]])
        except subprocess.CalledProcessError as error:
            cresil_failed(x, junctions, {"stage": "identify", "returncode": error.returncode,
                                         "error_tail": stderr_tail(x, "call")})
            return
        with (x.out / "eccDNA_final.txt").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"id", "ctc", "merge_region", "numreads"} <= set(reader.fieldnames or []):
                raise ValueError("CReSIL result schema mismatch")
            rows = list(reader)
        if any(r["ctc"] not in ("True", "False") for r in rows):
            raise ValueError("Unexpected CReSIL cyclic flag")
        x.run([sys.executable, REPO / "scripts" / "cresil_trace.py", x.out / "trim.txt", str(x.ref) + ".fai",
               x.out / "ont_junctions.tsv", x.out, x.out / "cresil_trace.tsv",
               "--minrsize", p["min_region_size"], "--depth", p["average_depth"]])
        with (x.out / "cresil_trace.tsv").open() as handle:
            trace = collections.Counter((r["location"], r["excluded_at"]) for r in csv.DictReader(handle, delimiter="\t"))
        write_json(x.state / "calls.json", {"caller": "ont_junctions+CReSIL",
            "classification": "evidence categories, not topology proof",
            "ont_junctions": junctions,
            "cresil": {"mode": p["cresil_mode"], "status": "ok", "cyclic_graphs": sum(r["ctc"] == "True" for r in rows), "all_graphs": len(rows)},
            "cresil_trace": {f"{k[0]}/{k[1]}": v for k, v in sorted(trace.items())},
            # Backwards-compatible field: CReSIL cyclic graphs.
            "candidate_count": sum(r["ctc"] == "True" for r in rows)})


def library_depth(paf, fai, out_tsv, reads_fa=None):
    """Read-length and per-contig depth summary from the whole-library PAF (primary chains).

    Unmapped reads never appear in a PAF, so input read counts, N50 and the mapped fractions
    come from reads_fa (all prepared reads) when given.

    Per-contig depth against the BioSample's assembly is the data-side check that reads and
    assembly come from the same culture: a plasmid contig without reads (or extra coverage
    elsewhere) flags a different isolate / DNA prep."""
    lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fai.read_text().splitlines() if l.strip()}
    aligned = collections.Counter()
    qlen, mapped = {}, set()
    with paf.open() as handle:
        for line in handle:
            f = line.split("\t", 13)
            qlen[f[0]] = int(f[1])
            if "tp:A:P" in line and int(f[11]) >= 1:
                aligned[f[5]] += int(f[8]) - int(f[7]); mapped.add(f[0])
    total = sum(lengths.values())
    depth = {c: aligned[c] / L for c, L in lengths.items()}
    by_len = sorted(lengths, key=lambda c: -lengths[c])
    med = sorted(depth[c] for c in by_len if lengths[c] >= 100000) or sorted(depth.values())
    median = med[len(med) // 2] if med else 0.0
    with out_tsv.open("w") as h:
        h.write("contig\tlength\taligned_bases\tmean_depth\tdepth_vs_large_contig_median\n")
        for c in by_len:
            h.write(f"{c}\t{lengths[c]}\t{aligned[c]}\t{depth[c]:.2f}\t{(depth[c] / median if median else 0):.3f}\n")
    if reads_fa is not None:
        ls, cur = [], -1
        with reads_fa.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    if cur >= 0:
                        ls.append(cur)
                    cur = 0
                else:
                    cur += len(line.rstrip("\n"))
        if cur >= 0:
            ls.append(cur)
        ls.sort(reverse=True)
    else:
        ls = sorted(qlen.values(), reverse=True)
    half, run, n50 = sum(ls) / 2, 0, 0
    for v in ls:
        run += v
        if run >= half:
            n50 = v; break
    return {"reads_input": len(ls), "bases_input": sum(ls), "reads_in_paf": len(qlen), "reads_mapped_primary": len(mapped),
            "mapped_fraction": round(len(mapped) / max(1, len(ls)), 4),
            "aligned_base_fraction": round(sum(aligned.values()) / max(1, sum(ls)), 4), "read_n50": n50,
            "mean_read_length": round(sum(ls) / max(1, len(ls)), 1), "assembly_bp": total,
            "mean_depth_assembly": round(sum(aligned.values()) / max(1, total), 2),
            "median_depth_large_contigs": round(median, 2),
            "contigs_below_0.1x_median_ge5kb": sum(1 for c in lengths if lengths[c] >= 5000 and median and depth[c] < 0.1 * median)}


def ont_junctions(x):
    """Direct long-read junction counts (scripts/ont_junctions.py), independent of CReSIL."""
    q = shlex.quote
    paf = x.out / "junctions.paf"
    # Secondary alignments (score >= 0.9 x primary) only flag repeat / multi-mapping anchors.
    x.shell(f"minimap2 -c --secondary=yes -N 5 -p 0.9 -x map-ont -t {max(1, x.cpus - 1)} {q(str(x.ref))} "
            f"{q(str(x.work / 'prepared' / 'reads.fa'))} > {q(str(paf))}")
    flags = []
    for key, value in x.c["caller"].get("ont_junction", {}).items():
        flags += ["--" + key.replace("_", "-"), str(value)]
    x.run([sys.executable, REPO / "scripts" / "ont_junctions.py", paf, x.ref, x.out / "ont_junctions.tsv",
           "--events", x.out / "ont_junction_events.tsv", "--support-ids", x.out / "ont_junction_ids.txt", *flags])
    write_json(x.out / "ont_library.json", library_depth(paf, Path(str(x.ref) + ".fai"), x.out / "contig_depth.tsv",
                                                             x.work / "prepared" / "reads.fa"))
    with (x.out / "ont_junction_ids.txt").open() as handle:
        ids = {l.strip() for l in handle if l.strip()}
    # Junction reads' alignments (incl. secondaries) are the classification basis; keep them.
    with paf.open() as source, gzip.open(x.out / "ont_junctions.support.paf.gz", "wt") as target:
        for line in source:
            if line.split("\t", 1)[0] in ids:
                target.write(line)
    paf.unlink()  # whole-library alignment: never archived, regenerated from the raw run if needed
    with (x.out / "ont_junctions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    context = junction_context(x)
    by = collections.Counter(f"{r['location']}/{r['label']}" for r in rows)
    thresholds = {f">={t}": collections.Counter(f"{r['location']}/{r['label']}" for r in rows
                                                if int(r["independent_molecules"]) >= t) for t in (1, 2, 3, 5, 10)}
    return {"junctions": len(rows), "by_location_label": dict(sorted(by.items())),
            "by_min_independent_molecules": {k: dict(sorted(v.items())) for k, v in thresholds.items()},
            "primary_threshold": "independent_molecules >= 3",
            "circle_compatible_yes": sum(r["circle_compatible"] == "yes" for r in rows),
            "junction_context": context,
            "note": "evidence categories, not topology proof; whole_contig = replicon closures (technical control)"}


CONTEXT_FLANK = 1000
JUNCTION_CONTEXT = ("junction_context.tsv.gz", "junction_context.fa.gz", "junction_context.ge3.fa.gz")


def junction_context(x):
    """[upstream][circle][downstream] of every junction from the sample's own assembly
    (scripts/junction_context.py, CONTEXT_FLANK bp flanks, genome copies of the circle sequence).
    It only reads ont_junctions.tsv and the reference. A failure is recorded, never fatal: the
    archive keeps ont_junctions.tsv and the reference SHA-256, so the context can be regenerated."""
    try:
        x.run([sys.executable, REPO / "scripts" / "junction_context.py", "--junctions", x.out / "ont_junctions.tsv",
               "--reference", x.ref, "--sample", x.s["sample_id"], "--out", x.out / "junction_context",
               "--flank", CONTEXT_FLANK])
    except subprocess.CalledProcessError as error:
        for name in JUNCTION_CONTEXT:
            (x.out / name).unlink(missing_ok=True)
        return {"status": "failed", "returncode": error.returncode,
                "error_tail": stderr_tail(x, os.environ.get("CS_RULE", "call"))}
    with gzip.open(x.out / JUNCTION_CONTEXT[0], "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {"status": "ok", "files": list(JUNCTION_CONTEXT), "flank_bp": CONTEXT_FLANK, "junctions": len(rows),
            "multi_copy": sum(r["genome_copies"] not in ("", "0", "1") for r in rows),
            "called_locus_not_found": sum(r["called_locus_in_copies"] == "no" for r in rows)}


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
                "eccDNA_final.txt", "no_candidates.json", "cresil_run/subGraphs.summary.txt",
                "ont_junctions.tsv", "ont_junction_events.tsv", "ont_junction_ids.txt",
                "ont_junctions.support.paf.gz", "cresil_trace.tsv", "ont_library.json", "contig_depth.tsv",
                *JUNCTION_CONTEXT}
            # final_reads.fa holds read sequence for EVERY graph (whole replicons in deep samples:
            # gigabytes); keep graph -> read ID tables instead (graph_reads.tsv.gz, below).
            keep = keep or (len(rel.parts) == 4 and rel.parts[:2] == ("cresil_run", "assemGraph")
                            and p.name == "reference_regions.fa")
            dest = ev / "cresil" / rel
        if keep:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest)
            selected.append({"path": rel.as_posix(), "bytes": p.stat().st_size})
        else:
            omitted_bytes += p.stat().st_size
    trim_rows = 0
    if platform == "nanopore" and (out / "trim.txt").is_file():  # absent when CReSIL failed early
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
        cyclic = set()
        if (x.out / "eccDNA_final.txt").is_file():
            with (x.out / "eccDNA_final.txt").open() as handle:
                cyclic = {r["id"] for r in csv.DictReader(handle, delimiter="\t") if r["ctc"] == "True"}
        junction_ids = x.out / "ont_junction_ids.txt"
        if junction_ids.is_file():  # reads behind every junction cluster, whatever its label: all kept
            wanted.update(l.strip() for l in junction_ids.open() if l.strip())
        # A cyclic graph over a high-copy replicon can hold most of the library (99% of the
        # support reads in the largest 1 TB-run archives). All graph read IDs go to
        # graph_reads.tsv.gz; sequences are kept for at most `cap` reads per cyclic graph,
        # chosen deterministically (junction reads first, then by a hash of the read ID).
        cap = int(x.c["caller"].get("cresil_graph_read_cap", 50))
        members = collections.defaultdict(set)
        (ev / "cresil").mkdir(parents=True, exist_ok=True)
        with gzip.open(ev / "cresil" / "graph_reads.tsv.gz", "wt") as table:
            table.write("graph_id\tcyclic\tread_id\tq_start\tq_end\n")
            for p in sorted(x.out.glob("cresil_run/assemGraph/*/final_reads.fa")):
                with p.open() as handle:
                    for line in handle:
                        if line.startswith(">"):
                            # Upstream format: graphID_originalReadID_queryStart_queryEnd.
                            header = line[1:].split()[0]
                            gid, rest = header.split("_", 1)
                            rid, qs, qe = rest.rsplit("_", 2)
                            table.write(f"{gid}\t{gid in cyclic}\t{rid}\t{qs}\t{qe}\n")
                            if gid in cyclic:
                                members[gid].add(rid)
        junction_reads = set(wanted)
        graph_kept = {}
        for gid in sorted(members):
            ranked = sorted(members[gid], key=lambda r: (r not in junction_reads, hashlib.sha256(r.encode()).hexdigest()))
            keep = ranked[:cap]
            wanted.update(keep)
            graph_kept[gid] = {"reads": len(members[gid]), "sequences_kept": len(keep)}
        write_json(ev / "cresil" / "graph_read_sequences.json", {"cap_per_cyclic_graph": cap, "graphs": graph_kept})
        (ev / "read_ids.txt").write_text("".join(n+"\n" for n in sorted(wanted)))
    copy_caller_evidence(x.out, ev, x.s["platform"], wanted)
    counts = extract(x.inputs, wanted, ev, x.s["platform"])
    if x.s["platform"] == "nanopore":
        # Preserve original read alignments with CIGAR/SA/secondary records for downstream review.
        q = shlex.quote
        if wanted:
            x.shell(f"minimap2 -ax map-ont --secondary=yes -t {max(1,x.cpus-2)} {q(str(x.ref))} {q(str(ev/'support_reads1.fastq.gz'))} | samtools sort -@ 1 -o {q(str(ev/'support.bam'))}")
            # Every read has exactly one primary (or unmapped) record carrying its full sequence
            # and qualities, so the FASTQ is redundant (about half the archive). Verify before
            # removing: `samtools fastq -F 0x900 support.bam` recovers it.
            with gzip.open(ev / "support_reads1.fastq.gz", "rt") as handle:
                fastq_reads = sum(1 for i, _ in enumerate(handle) if i % 4 == 0)
            count = x.work / "support.primary.count"
            x.run(["samtools", "view", "-c", "-F", "0x900", ev / "support.bam"], stdout=count)
            primary = int(count.read_text().strip()); count.unlink()
            if primary != fastq_reads:
                raise ValueError(f"support.bam primary records ({primary}) != support reads ({fastq_reads})")
            (ev / "support_reads1.fastq.gz").unlink()
        else:
            (ev / "support_reads1.fastq.gz").unlink(missing_ok=True)
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
        "scope": "all split/discordant read IDs" if x.s["platform"]=="illumina" else "all junction-counter reads (every junction, any label) + up to cresil_graph_read_cap reads per CReSIL cyclic graph (cresil/graph_read_sequences.json); read IDs of all graphs in cresil/graph_reads.tsv.gz",
        **({} if x.s["platform"] == "illumina" else {"support_reads": "support.bam only: samtools fastq -F 0x900 support.bam recovers every support read with qualities"}),
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
