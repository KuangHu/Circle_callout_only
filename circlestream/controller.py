"""Separate commands for DTN download and Slurm submission; no remote SSH calls."""
import argparse
import csv
import getpass
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import tarfile
import time

from .common import (GB, REPO, digest, fingerprint, is_sra, load_config, lock, read_json,
                     local_mode, safe_id, set_status, source_ok, state_paths, tree_size, validate_sample, work_root,
                     write_json, workflow_fingerprint)
from .sra import fetch_refseqs, stage_refseqs

ACTIVE_JOBS = {"SUBMITTING", "QUEUED", "RUNNING"}
HELD = {"DOWNLOAD_FAILED", "COMPUTE_FAILED", "CLEANUP_FAILED"}
TERMINAL_JOBS = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                 "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}


def entries(c):
    root = Path(c["state_dir"])
    return sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".") and (p / "spec.json").is_file()) if root.exists() else []


def status(path):
    return read_json(path / "status.json")["status"]


def initialize(c, manifest):
    samples = read_json(manifest)
    ids = [s["sample_id"] for s in samples]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs")
    for sample in samples:
        validate_sample(sample)
    for key in ("state_dir", "scratch_dir", "archive_dir"):
        Path(c[key]).mkdir(parents=True, exist_ok=True)
    with lock(Path(c["state_dir"]) / ".init-lock"):
        config_snapshot = Path(c["state_dir"]) / "config.json"
        if config_snapshot.exists() and read_json(config_snapshot) != c:
            raise ValueError("Config changed: use a new state/scratch/archive set for a new benchmark profile")
        write_json(config_snapshot, c)
        for sample in samples:
            state, scratch = state_paths(c, sample["sample_id"])
            spec = {"config": c, "sample": sample, "repo": str(REPO), "python": sys.executable,
                    "workflow_sha256": workflow_fingerprint(), "reference_sha256": digest(sample["reference"])}
            if c.get("reference_store"):
                spec["reference_store_path"] = store_reference(c, sample["reference"], spec["reference_sha256"])
            if state.exists():
                if read_json(state / "spec.json") != spec:
                    raise ValueError(f"Changed sample or software location: {state.name}; use a new ID/run root")
                continue
            if scratch.exists():
                raise ValueError(f"Unowned scratch directory already exists: {scratch}")
            pending = state.with_name(".init-"+state.name)
            pending.mkdir()
            write_json(pending / "spec.json", spec)
            set_status(pending, "NEW")
            os.replace(pending,state)
    print(f"Initialized {len(samples)} samples; no downloads or jobs started")


def store_reference(c, path, sha):
    """Permanent, deduplicated copy named by checksum; archives refer to it instead of embedding it."""
    store = Path(c["reference_store"])
    store.mkdir(parents=True, exist_ok=True)
    dest = store / f"{sha}.fa"
    if dest.exists():
        if digest(dest) != sha:
            raise ValueError(f"Reference store entry is corrupted: {dest}")
        return str(dest)
    temp = store / f".{sha}.{os.getpid()}.tmp"
    shutil.copyfile(path, temp)
    if digest(temp) != sha:
        temp.unlink()
        raise ValueError(f"Reference changed while copying into store: {path}")
    os.chmod(temp, 0o444)
    os.replace(temp, dest)
    return str(dest)


def check_snapshot(c):
    if read_json(Path(c["state_dir"]) / "config.json") != c:
        raise ValueError("Configuration differs from init snapshot; create a separate benchmark run")


def dependency_bytes(state):
    """Bytes of refseq dependencies a fetch job staged for this sample (0 until recorded)."""
    d = state / "download.json"
    if not d.is_file():
        return 0
    return sum(r.get("size_bytes") or 0 for r in read_json(d).get("refseqs", []))


def raw_released(state):
    return bool(read_json(state / "status.json").get("raw_released"))


def admitted(c, sample):
    # Failed samples whose raw objects were released (release_failed_raw) hold no disk.
    active = [p for p in entries(c) if status(p) not in ("NEW", "CLEANED") and not raw_released(p)]
    reserved = sum(read_json(p / "spec.json")["sample"]["peak_gb"] for p in active)
    # Held samples retain their full disk reservation, but do not monopolize slots.
    live = [p for p in active if status(p) not in HELD]
    if local_mode(c):
        # Shared scratch only holds source objects; node-local disk holds the work.
        raw = sum(f["size_bytes"] for p in active for f in read_json(p / "spec.json")["sample"]["files"])
        # Aligned SRA runs also stage refseq dependencies beside the object; their sizes are only
        # known once fetched (download.json), so charge them from then on.
        raw += sum(dependency_bytes(p) for p in active)
        need = sum(f["size_bytes"] for f in sample["files"])
        if len(live) >= c["max_active_samples"] or (raw + need) / GB > c["scratch_budget_gb"]:
            return False
        return shutil.disk_usage(c["scratch_dir"]).free / GB >= need / GB + c["min_free_gb"]
    if len(live) >= c["max_active_samples"] or reserved + sample["peak_gb"] > c["scratch_budget_gb"]:
        return False
    # Account for remaining reserved expansion, not just currently free bytes.
    used = tree_size(c["scratch_dir"]) / GB
    needed = max(0, reserved - used) + sample["peak_gb"] + c["min_free_gb"]
    return shutil.disk_usage(c["scratch_dir"]).free / GB >= needed


def download_file(item, raw, retries):
    final = raw / item["name"]
    partial = raw / (item["name"] + ".part")
    def valid(path):
        return source_ok(path, item)
    if valid(final):
        return
    if final.exists():
        raise ValueError(f"Previously finalized input is corrupted: {final}")
    for attempt in range(retries):
        if partial.exists() and valid(partial):
            os.replace(partial, final)
            return
        if partial.exists() and partial.stat().st_size >= item["size_bytes"]:
            partial.unlink()  # only this managed, checksum-invalid download
        result = subprocess.run([
            "curl", "--fail", "--location", "--proto", "=https", "--proto-redir", "=https",
            "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "120",
            "--max-filesize", str(item["size_bytes"]), "--continue-at", "-",
            "--output", str(partial), item["url"]
        ])
        if result.returncode == 0 and valid(partial):
            os.replace(partial, final)
            return
        if result.returncode == 33 and partial.exists():
            partial.unlink()  # server does not support range resume
        if attempt + 1 < retries:
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"Download/checksum failed after {retries} attempts: {item['name']}")


def require_dtn(c):
    actual = socket.getfqdn().split(".")[0]
    if actual != c["download_host"].split(".")[0]:
        raise RuntimeError(f"Run download on configured DTN {c['download_host']}, not {socket.getfqdn()}")


def claim_scratch(state, scratch, spec):
    """Create owned scratch for an admitted NEW sample (caller holds the admission lock).
    An existing directory is reused only if this exact spec owns it (download resume)."""
    owner = {"spec_sha256": fingerprint(spec)}
    if scratch.exists():
        if read_json(scratch / ".owner.json") != owner:
            raise ValueError(f"Unowned scratch directory already exists: {scratch}")
    else:
        scratch.mkdir()
        write_json(scratch / ".owner.json", owner)
    (scratch / "raw").mkdir(exist_ok=True)


def fetch_inputs(c, sample, state, scratch, **extra):
    """Download (or resume) source files and SRA refseq dependencies into owned scratch and
    record timing. Every download path (fetch job, DTN, DTN resume) goes through here."""
    start = time.time()
    refseqs = []
    for item in sample["files"]:
        download_file(item, scratch / "raw", c["download_retries"])
        if is_sra(item):
            # Aligned-read runs need their reference sequences beside the object (see sra.py).
            refseqs += fetch_refseqs(scratch / "raw" / item["name"], c["download_retries"])
    write_json(state / "download.json", {
        "refseqs": refseqs,
        "seconds": time.time() - start,
        "compressed_bytes": sum(f["size_bytes"] for f in sample["files"]),
        "host": socket.getfqdn(), "in_job": bool(os.environ.get("SLURM_JOB_ID")),
        "files": sample["files"], **extra
    })


def download_one(c, sid=None):
    if local_mode(c):
        raise RuntimeError("download_mode is 'slurm': 'submit --watch' launches 1-CPU fetch jobs; no DTN process")
    require_dtn(c)
    with lock(Path(c["state_dir"]) / ".download-lock"):
        for state in entries(c):
            if sid and state.name != sid:
                continue
            if status(state) != "NEW":
                continue
            spec = read_json(state / "spec.json")
            sample = spec["sample"]
            if not admitted(c, sample):
                continue
            _, scratch = state_paths(c, sample["sample_id"])
            claim_scratch(state, scratch, spec)
            set_status(state, "DOWNLOADING")
            try:
                fetch_inputs(c, sample, state, scratch)
                set_status(state, "READY")
            except Exception as error:
                set_status(state, "DOWNLOAD_FAILED", error=str(error))
                raise
            print(f"Downloaded {sample['sample_id']}", flush=True)
            return True
    return False


def fetch_records(c):
    d = Path(c["state_dir"]) / ".fetch"
    return sorted(d.glob("fetch-*.json")) if d.exists() else []


def reconcile_fetch(c):
    """Fetch jobs that ended leave no sample in DOWNLOADING (held, not retried)."""
    records = fetch_records(c)
    if not records:
        return
    ids = {read_json(r)["job_id"]: r for r in records}
    try:
        acct = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ",".join(ids), "--format=JobIDRaw,State%40"],
                              check=True, capture_output=True, text=True, timeout=30)
        queue = subprocess.run(["squeue", "--noheader", "--user", getpass.getuser(), "--format=%A"],
                               check=True, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        print(f"Scheduler reconciliation unavailable; retaining fetch slots: {error}", file=sys.stderr, flush=True)
        return
    terminal = {}
    for line in acct.stdout.splitlines():
        f = line.split("|")
        if len(f) >= 2 and f[1].strip():
            terminal[f[0].strip()] = f[1].split()[0].rstrip("+")
    live = set(queue.stdout.split())
    for jid, record in ids.items():
        if terminal.get(jid) not in TERMINAL_JOBS or jid in live:
            continue
        for sid in read_json(record)["samples"]:
            state = Path(c["state_dir"]) / sid
            if status(state) == "DOWNLOADING":
                set_status(state, "DOWNLOAD_FAILED", fetch_job=jid, error=f"Fetch job ended ({terminal[jid]}) before this sample finished")
                print(f"Held {sid}: fetch job {jid} {terminal[jid]}", flush=True)
        os.replace(record, record.with_suffix(".done"))


def submit_fetch(c):
    """Launch one 1-CPU fetch job for a group of admitted NEW samples (slurm mode)."""
    validate_slurm(c)
    f, s = c["fetch"], c["slurm"]
    with lock(Path(c["state_dir"]) / ".submit-lock"):
        reconcile_fetch(c)
        free_slots = f["max_jobs"] - len(fetch_records(c))
        if free_slots <= 0:
            return False
        # Spread admissible samples over all free fetch slots (each fetch job downloads
        # serially), so several fetch jobs feed processing concurrently.
        states = entries(c)
        live = sum(1 for p in states if status(p) not in ("NEW", "CLEANED") and status(p) not in HELD)
        share = -(-max(0, c["max_active_samples"] - live) // free_slots)
        cap = max(1, min(f.get("max_samples", 10), share))
        chosen, size = [], 0
        for state in states:
            if status(state) != "NEW":
                continue
            spec = read_json(state / "spec.json")
            need = sum(x["size_bytes"] for x in spec["sample"]["files"])
            if len(chosen) >= cap or (chosen and size + need > f["batch_gb"] * GB):
                break
            if not admitted(c, spec["sample"]):
                break
            claim_scratch(state, state_paths(c, state.name)[1], spec)
            set_status(state, "DOWNLOADING", fetch_job="pending")
            chosen.append(state)
            size += need
        if not chosen:
            return False
        d = Path(c["state_dir"]) / ".fetch"
        d.mkdir(exist_ok=True)
        tag = f"fetch-{time.time_ns()}"
        python = read_json(chosen[0] / "spec.json")["python"]
        cmd = [python, str(REPO / "stream.py"), "fetch"] + [x for st in chosen for x in ("--spec", str(st / "spec.json"))]
        script = d / f"{tag}.sh"
        script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(s.get("setup_commands", [])) +
                          "\nexec " + shlex.join(cmd) + "\n", encoding="utf-8", newline="\n")
        result = subprocess.run([
            "sbatch", "--parsable", "--account", s["account"], "--partition", s["partition"], "--qos", s["qos"],
            "--cpus-per-task", "1", "--mem", f"{f['memory_gb']}G", "--time", f["time"],
            "--job-name", "csf_" + fingerprint(tag)[:12], "--output", str(d / f"{tag}-%j.out"), str(script)
        ], text=True, capture_output=True)
        jobid = result.stdout.strip().split(";")[0]
        if result.returncode != 0 or not jobid.isdigit():
            for st in chosen:
                set_status(st, "DOWNLOAD_FAILED", error=f"Fetch submission failed: {result.stderr.strip()[:300]}")
            raise RuntimeError(f"Fetch job submission failed: {result.stderr.strip()}")
        write_json(d / f"{tag}.json", {"job_id": jobid, "samples": [st.name for st in chosen],
                                       "bytes": size, "submitted": time.time()})
        for st in chosen:
            if status(st) == "DOWNLOADING":
                set_status(st, "DOWNLOADING", fetch_job=jobid)
        print(f"Fetch job {jobid}: {len(chosen)} samples, {size/1e9:.2f} GB", flush=True)
        return True


def fetch_specs(spec_paths):
    """Runs inside a 1-CPU fetch job: download each sample in turn; one failure holds only that sample."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("fetch must run inside a Slurm allocation")
    failed = 0
    for path in spec_paths:
        spec = read_json(path)
        c, sample = spec["config"], spec["sample"]
        state, scratch = state_paths(c, sample["sample_id"])
        if status(state) != "DOWNLOADING":
            print(f"Skipping {state.name}: {status(state)}", flush=True)
            continue
        try:
            if read_json(scratch / ".owner.json") != {"spec_sha256": fingerprint(spec)}:
                raise ValueError("Scratch ownership mismatch")
            fetch_inputs(c, sample, state, scratch)
            set_status(state, "READY")
            print(f"Fetched {state.name}", flush=True)
        except Exception as error:
            failed += 1
            set_status(state, "DOWNLOAD_FAILED", error=str(error))
            print(f"Held {state.name}: {error}", file=sys.stderr, flush=True)
    return failed


def validate_slurm(c):
    for key in ("account", "partition", "qos"):
        value = c["slurm"][key]
        if not value or "CHANGE_ME" in value or any(x.isspace() for x in value):
            raise ValueError(f"Set Slurm {key}")
    setup = c["slurm"].get("setup_commands", [])
    if not isinstance(setup, list) or any(not isinstance(line, str) or "\n" in line or "\r" in line for line in setup):
        raise ValueError("slurm.setup_commands must be a list of single-line shell commands")


def reconcile_jobs(c):
    """Called under the submit lock; release slots only on confirmed scheduler exit.

    Never retry jobs or remove worker locks. An unknown submission/accounting state
    remains active so that an accounting outage cannot create duplicate work.
    """
    jobs = {}
    for state in entries(c):
        if status(state) in ACTIVE_JOBS and (state / "job.json").is_file():
            jid = read_json(state / "job.json")["job_id"]
            if not isinstance(jid, str) or not jid.isdigit():
                raise ValueError(f"Invalid recorded job ID: {state.name}")
            jobs[jid] = state
    if not jobs:
        return
    try:
        result = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ",".join(jobs),
                                 "--format=JobIDRaw,State%40"],
                                check=True, capture_output=True, text=True, timeout=30)
        terminal = {}
        for line in result.stdout.splitlines():
            fields = line.split("|")
            if len(fields) >= 2 and fields[1].strip():
                terminal[fields[0].strip()] = fields[1].split()[0].rstrip("+")
        # A requeued/preempted job may still be live despite an earlier terminal
        # accounting record. Absence from the queue is required as well.
        queue = subprocess.run(["squeue", "--noheader", "--user", getpass.getuser(), "--format=%A"],
                               check=True, capture_output=True, text=True, timeout=30)
        live = set(queue.stdout.split())
    except (OSError, subprocess.SubprocessError) as error:
        print(f"Scheduler reconciliation unavailable; retaining job slots: {error}", file=sys.stderr, flush=True)
        return
    for jid, state in jobs.items():
        if terminal.get(jid) in TERMINAL_JOBS and jid not in live and status(state) in ACTIVE_JOBS:
            set_status(state, "COMPUTE_FAILED", job_id=jid, scheduler_state=terminal[jid],
                       error="Job ended without worker final status; inspect retained scratch/locks before recovery")
            print(f"Held {state.name}: Slurm {terminal[jid]}; job slot released, scratch retained", flush=True)


def only_held_left(c):
    """True when no sample can progress without operator action (all CLEANED or held)."""
    states = [status(p) for p in entries(c)]
    return bool(states) and all(st == "CLEANED" or st in HELD for st in states) and any(st in HELD for st in states)


def release_failed_raw(c):
    """Opt-in (release_failed_raw, slurm mode). Held samples keep their raw objects on shared
    staging, and those bytes count against scratch_budget_gb: in the 1 TB run, 32 failures held
    43 GB and admission stalled. Once a COMPUTE_FAILED sample's job has left the queue, delete
    its owned staging directory; `recover` then re-fetches (checksum-gated) instead of resuming.

    Returns the number of failed samples whose release is deferred (worker still finishing, job
    still queued, or a concurrent recover holding the sample's .release-lock), so the watcher
    keeps running until they can be released."""
    if not (local_mode(c) and c.get("release_failed_raw")):
        return 0
    failed = [p for p in entries(c) if status(p) == "COMPUTE_FAILED" and not raw_released(p)]
    if not failed:
        return 0
    try:
        queue = subprocess.run(["squeue", "--noheader", "--user", getpass.getuser(), "--format=%A"],
                               check=True, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        print(f"Queue unavailable; not releasing failed raw: {error}", file=sys.stderr, flush=True)
        return len(failed)
    live = set(queue.stdout.split())
    root = Path(c["scratch_dir"]).resolve()
    deferred = 0
    for state in failed:
        try:
            # recover() takes the same per-sample lock, so the status cannot change between
            # these checks and the deletion (no READY + new worker on a deleted directory).
            with lock(state / ".release-lock"):
                if status(state) != "COMPUTE_FAILED" or raw_released(state):
                    continue
                job = state / "job.json"
                if (any((state / name).exists() for name in (".worker-lock", ".cleanup-lock"))
                        or (job.exists() and read_json(job)["job_id"] in live)):
                    deferred += 1  # the worker may still be finishing on its node
                    continue
                scratch = Path(c["scratch_dir"]) / state.name
                if scratch.exists():
                    spec = read_json(state / "spec.json")
                    if (scratch.is_symlink() or scratch.resolve().parent != root
                            or read_json(scratch / ".owner.json") != {"spec_sha256": fingerprint(spec)}
                            or any(p.is_symlink() for p in scratch.rglob("*"))):
                        print(f"Not releasing {state.name}: scratch ownership/symlink check failed", file=sys.stderr, flush=True)
                        continue
                    shutil.rmtree(scratch)
                info = read_json(state / "status.json")
                info.pop("updated", None)
                set_status(state, info.pop("status"), raw_released=True, **info)
                print(f"Released raw staging of failed {state.name}; recover re-fetches it", flush=True)
        except RuntimeError:
            deferred += 1  # .release-lock held by a concurrent recover
    return deferred


def cleanup_archived(c):
    """Isolate cleanup failures so a single sample cannot terminate the watcher."""
    for state in entries(c):
        if status(state) != "ARCHIVED" or any((state / name).exists() for name in (".worker-lock", ".cleanup-lock")):
            continue
        try:
            cleanup(c, state.name)
        except Exception as error:
            # A concurrent explicit cleanup may have acquired its per-sample lock.
            if status(state) == "ARCHIVED" and not (state / ".cleanup-lock").exists():
                set_status(state, "CLEANUP_FAILED", error=str(error))
            print(f"Cleanup held for {state.name}: {error}; continuing other samples", file=sys.stderr, flush=True)


def submit_one(c, through="archive", sid=None):
    validate_slurm(c)
    with lock(Path(c["state_dir"]) / ".submit-lock"):
        reconcile_jobs(c)
        active = [p for p in entries(c) if status(p) in ACTIVE_JOBS]
        if len(active) >= c["max_jobs"]:
            return False
        for state in entries(c):
            if sid and state.name != sid:
                continue
            if status(state) != "READY":
                continue
            spec = read_json(state / "spec.json")
            script = state / "job.sh"
            cmd = [spec["python"], str(REPO / "stream.py"), "worker", "--spec", str(state / "spec.json"), "--through", through]
            setup = "\n".join(c["slurm"].get("setup_commands", []))
            script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + setup + "\nexec " + shlex.join(cmd) + "\n", encoding="utf-8", newline="\n")
            s = c["slurm"]
            name = "cs_" + fingerprint(str(state))[:16]
            set_status(state, "SUBMITTING", job_name=name)
            # A crash/timeout after sbatch can mean the job exists: never blindly resubmit.
            result = subprocess.run([
                "sbatch", "--parsable", "--account", s["account"], "--partition", s["partition"],
                "--qos", s["qos"], "--cpus-per-task", str(s["cpus"]), "--mem", f"{s['memory_gb']}G",
                "--time", s["time"], "--job-name", name, "--output", str(state / "slurm-%j.out"), str(script)
            ], text=True, capture_output=True)
            (state / "submission.stderr").write_text(result.stderr, encoding="utf-8")
            if result.returncode != 0:
                raise RuntimeError(f"Submission uncertain/failed for {state.name}; inspect SUBMITTING state before recovery")
            jobid = result.stdout.strip().split(";")[0]
            if not jobid.isdigit():
                raise RuntimeError(f"Unexpected sbatch response: {result.stdout}")
            write_json(state / "job.json", {"job_id": jobid, "job_name": name, "submitted": time.time(), "through": through})
            # Worker may already have started or completed. Do not overwrite its status.
            print(f"Submitted {state.name}: {jobid}", flush=True)
            return True
    return False


def stage_local(c, sid, inputs):
    """Copy existing cluster fixtures into owned scratch; never delete source inputs."""
    state, scratch = state_paths(c, sid)
    with lock(Path(c["state_dir"]) / ".download-lock"):
        if status(state) != "NEW":
            raise ValueError("stage-local requires NEW sample")
        spec = read_json(state / "spec.json")
        sample = spec["sample"]
        if len(inputs) != len(sample["files"]) or not admitted(c, sample):
            raise ValueError("Wrong input count or insufficient reserved capacity")
        for source,item in zip(inputs,sample["files"]):
            if not source_ok(source, item):
                raise ValueError("Local input does not match manifest")
        scratch.mkdir(exist_ok=False)
        write_json(scratch / ".owner.json", {"spec_sha256": fingerprint(spec)})
        (scratch / "raw").mkdir()
        set_status(state,"DOWNLOADING")
        start=time.time()
        try:
            refseqs = []
            for source,item in zip(inputs,sample["files"]):
                target=scratch / "raw" / item["name"]
                shutil.copy2(source,target)
                if not source_ok(target, item):
                    raise ValueError("Staged input checksum mismatch")
                if is_sra(item):  # aligned runs: refseqs beside the source, else fetched
                    refseqs += stage_refseqs(source, target, c["download_retries"])
            write_json(state/"download.json", {"seconds":time.time()-start,"local_copy":True,"refseqs":refseqs,
                "compressed_bytes":sum(f["size_bytes"] for f in sample["files"])})
            set_status(state,"READY")
        except Exception as error:
            set_status(state,"DOWNLOAD_FAILED",error=str(error))
            raise


def worker(spec_path, through):
    spec = read_json(spec_path)
    c, sample = spec["config"], spec["sample"]
    state, scratch = state_paths(c, sample["sample_id"])
    if spec["workflow_sha256"] != workflow_fingerprint():
        raise ValueError("Workflow changed since init; use a separate run root rather than reuse old rule markers")
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("worker must run inside a Slurm allocation; unit tests do not invoke real callers")
    with lock(state / ".worker-lock"):
        if status(state) not in ("SUBMITTING", "READY"):
            raise RuntimeError(f"Cannot start worker from {status(state)}")
        set_status(state, "RUNNING", job_id=os.environ["SLURM_JOB_ID"], started=time.time(),
                   host=socket.getfqdn())
        root = work_root(c, sample["sample_id"])
        if local_mode(c) and root.exists():
            raise RuntimeError(f"Node-local work directory already exists: {root}")
        work = root / "work"
        work.mkdir(parents=True, exist_ok=True)
        target = {"preflight": "preflight.ok", "qc": "qc.ok", "reference": "reference.ok", "mapping": "mapping.ok",
                  "call": "call.ok", "evidence": "evidence.ok", "archive": "archive.ok"}[through]
        cmd = [spec["python"], "-m", "snakemake", "--snakefile", str(REPO / "workflow" / "Snakefile"),
               "--directory", str(work), "--cores", str(c["slurm"]["cpus"]), "--rerun-incomplete",
               "--config", f"spec={Path(spec_path).resolve()}", "--", target]
        try:
            subprocess.run(cmd, check=True)
            set_status(state, "ARCHIVED" if through == "archive" else "PAUSED", through=through)
        except Exception as error:
            set_status(state, "COMPUTE_FAILED", error=str(error))
            raise
        finally:
            if local_mode(c):
                remove_local(root)
    if local_mode(c) and through == "archive" and c["auto_cleanup"]:
        # Verify the bundle and release the staged source object from inside the job.
        # With auto_cleanup=false the staged copy stays until an explicit 'cleanup'.
        try:
            cleanup(c, sample["sample_id"])
        except Exception as error:  # e.g. the login loop holds .cleanup-lock; it will finish
            print(f"In-job cleanup deferred for {sample['sample_id']}: {error}", file=sys.stderr, flush=True)


def remove_local(root):
    """Delete this job's node-local tree; logs and metrics already live in the state dir."""
    base = Path(os.environ.get("CS_LOCAL_ROOT", "/tmp")).resolve()
    job = root.parent
    if job.parent.resolve() != base or job.name != f"circle-stream-{os.environ.get('SLURM_JOB_ID')}" or job.is_symlink():
        raise ValueError(f"Refusing to remove unexpected local path {job}")
    shutil.rmtree(job, ignore_errors=True)


def verify_bundle(state):
    receipt = read_json(state / "archive.json")
    bundle = Path(receipt["path"])
    if digest(bundle) != receipt["sha256"]:
        raise ValueError("Archive checksum mismatch; refusing cleanup")
    with tarfile.open(bundle, "r:gz") as archive:
        manifest = json.load(archive.extractfile("MANIFEST.json"))
        if json.load(archive.extractfile("metadata/spec.json")) != read_json(state / "spec.json"):
            raise ValueError("Archive belongs to a different sample or run")
        names = [m.name for m in archive.getmembers()]
        if len(names) != len(set(names)) or set(names) != set(manifest) | {"MANIFEST.json"}:
            raise ValueError("Archive member set mismatch")
        import hashlib
        for name, expected in manifest.items():
            member = archive.getmember(name)
            if not member.isfile():
                raise ValueError("Non-file archive member")
            h = hashlib.sha256()
            with archive.extractfile(member) as handle:
                for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                    h.update(block)
            if h.hexdigest() != expected:
                raise ValueError(f"Archive member checksum mismatch: {name}")
    return receipt


def reconcile_archive(c, sid):
    state, scratch = state_paths(c,sid)
    with lock(state / ".cleanup-lock"):
        if (state/".worker-lock").exists():
            raise ValueError("Worker lock present: confirm job completion before archive recovery")
        old = status(state)
        if old not in ("ARCHIVED","CLEANED","CLEANUP_FAILED"):
            if not (state/"job.json").exists():
                raise ValueError("No job ID available; reconcile ambiguous submission first")
            jid=read_json(state/"job.json")["job_id"]
            result=subprocess.run(["sacct","-X","-n","-P","-j",jid,"--format=JobIDRaw,State"],check=True,capture_output=True,text=True)
            rows=[l.split("|") for l in result.stdout.splitlines() if l.split("|")[0]==jid]
            if not rows or rows[0][1].split()[0].rstrip("+") not in {"COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE"}:
                raise ValueError("Job active or unknown")
        receipt=read_json(state/"archive.json")
        final=Path(c["archive_dir"])/(sid+".tar.gz")
        partial=Path(c["archive_dir"])/(sid+".partial.tar.gz")
        if Path(receipt["path"]) not in (final,partial):
            raise ValueError("Unexpected archive path")
        if not Path(receipt["path"]).exists() and final.exists() and digest(final)==receipt["sha256"]:
            receipt["path"]=str(final);write_json(state/"archive.json",receipt)
        verify_bundle(state)
        if Path(receipt["path"])==partial:
            if final.exists():
                raise ValueError("Both partial and final exist; review before overwrite")
            os.replace(partial,final)
            receipt["path"]=str(final);write_json(state/"archive.json",receipt)
        if old!="CLEANED":
            set_status(state,"ARCHIVED",reconciled=True)
        print(f"Archive verified: {sid}")


def cleanup(c, sid):
    state, scratch = state_paths(c, sid)
    with lock(state / ".cleanup-lock"):
        if status(state) == "CLEANED":
            return
        if status(state) != "ARCHIVED":
            raise ValueError("Cleanup requires ARCHIVED status")
        verify_bundle(state)
        spec = read_json(state / "spec.json")
        root = Path(c["scratch_dir"]).resolve()
        if scratch.is_symlink() or scratch.resolve().parent != root or scratch.name != safe_id(sid):
            raise ValueError("Unsafe cleanup target")
        if scratch.exists():
            if read_json(scratch / ".owner.json") != {"spec_sha256": fingerprint(spec)}:
                raise ValueError("Scratch ownership mismatch")
            if (state / ".worker-lock").exists():
                raise ValueError("Worker still holds lock")
            # Reject all symlinks to avoid following caller-created links in cleanup.
            if any(p.is_symlink() for p in scratch.rglob("*")):
                raise ValueError("Symlink in scratch; manual review required")
            shutil.rmtree(scratch)
        set_status(state, "CLEANED")
        print(f"Cleaned {sid}; verified bundle retained", flush=True)


def recover(c, sid):
    """Explicit operator action, not automatic retry of ambiguous Slurm submissions."""
    state, _ = state_paths(c, sid)
    if status(state) == "COMPUTE_FAILED":
        # Serialize with release_failed_raw: it must not delete staging that this recovery is
        # about to hand to a new worker, and recovery must see whether raw was released.
        with lock(state / ".release-lock"):
            return _recover(c, sid)
    return _recover(c, sid)


def _recover(c, sid):
    state, scratch = state_paths(c, sid)
    current = status(state)
    for name in (".worker-lock", ".cleanup-lock"):
        if (state / name).exists():
            raise ValueError(f"Inspect stale lock first: {state / name}")
    if current == "DOWNLOAD_FAILED" and local_mode(c):
        # The next fetch job resumes into the same owned scratch (partial files are checksum-gated).
        set_status(state, "NEW", recovered_from="DOWNLOAD_FAILED")
        return
    if current in ("DOWNLOAD_FAILED", "DOWNLOADING"):
        require_dtn(c)
        if (Path(c["state_dir"]) / ".download-lock").exists():
            raise ValueError("Download lock is still present")
        with lock(Path(c["state_dir"]) / ".download-lock"):
            spec = read_json(state / "spec.json")
            if read_json(scratch / ".owner.json") != {"spec_sha256": fingerprint(spec)}:
                raise ValueError("Ownership mismatch")
            set_status(state, "DOWNLOADING")
            try:
                fetch_inputs(c, spec["sample"], state, scratch, resumed=True)  # incl. SRA refseqs
                set_status(state, "READY")
            except Exception as error:
                set_status(state, "DOWNLOAD_FAILED", error=str(error))
                raise
    elif current == "COMPUTE_FAILED" and raw_released(state):
        # Raw objects were released: download again through a fetch job, then rerun.
        set_status(state, "NEW", recovered_from="COMPUTE_FAILED", raw_refetch=True)
    elif current in ("COMPUTE_FAILED", "PAUSED", "RUNNING", "SUBMITTING"):
        job = state / "job.json"
        if not job.exists():
            raise ValueError("Missing job ID: inspect squeue/sacct using status.json job_name; reconcile manually")
        jid = read_json(job)["job_id"]
        check = subprocess.run(["sacct", "-X", "-n", "-P", "-j", jid, "--format=JobIDRaw,State"],
                               check=True, capture_output=True, text=True)
        rows = [l.split("|") for l in check.stdout.splitlines() if l.split("|")[0] == jid]
        terminal = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"}
        if not rows or rows[0][1].split()[0].rstrip("+") not in terminal:
            raise ValueError("Job is active or accounting unavailable: do not retry")
        # Preserve failed partial results and benchmarks. Rules rerun with isolated caller outputs.
        set_status(state, "READY", recovered_job=jid)
    else:
        raise ValueError(f"Cannot recover state {current}")


def report(c):
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", "platform", "status", "download_seconds", "input_bytes", "rule", "wall_seconds", "cpu_seconds", "max_rss_kb", "sampled_peak_scratch_bytes", "exit_code"])
    for state in entries(c):
        spec = read_json(state / "spec.json")
        d = read_json(state / "download.json") if (state / "download.json").exists() else {}
        metrics = sorted((state / "metrics").glob("*.json")) or [None]
        for m in metrics:
            r = read_json(m) if m else {}
            writer.writerow([state.name, spec["sample"]["platform"], status(state), d.get("seconds", ""),
                d.get("compressed_bytes", ""), r.get("rule", ""), r.get("wall_seconds", ""), r.get("cpu_seconds", ""),
                r.get("max_rss_kb", ""), r.get("sampled_peak_scratch_bytes", ""), r.get("exit_code", "")])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=["init", "download", "stage-local", "submit", "cleanup", "verify", "recover", "status", "report", "worker", "fetch"])
    ap.add_argument("--config", default="config.local.json")
    ap.add_argument("--manifest", default="samples.local.json")
    ap.add_argument("--sample")
    ap.add_argument("--spec", action="append")
    ap.add_argument("--inputs", nargs="+")
    ap.add_argument("--watch", action="store_true")
    ap.add_argument("--through", choices=["preflight", "qc", "reference", "mapping", "call", "evidence", "archive"], default="archive")
    a = ap.parse_args()
    if a.command == "worker":
        worker(a.spec[0], a.through)
        return
    if a.command == "fetch":
        sys.exit(1 if fetch_specs(a.spec) else 0)
    c = load_config(a.config)
    if a.command == "init":
        initialize(c, a.manifest)
        return
    check_snapshot(c)
    if a.sample:
        safe_id(a.sample)
        if not (Path(c["state_dir"])/a.sample/"spec.json").is_file():
            raise ValueError("Unknown sample ID")
    if a.command == "stage-local":
        if not a.sample or not a.inputs:
            ap.error("stage-local requires --sample and --inputs")
        stage_local(c,a.sample,a.inputs)
        return
    if a.command == "status":
        for p in entries(c):
            print(p.name, json.dumps(read_json(p / "status.json")))
        return
    if a.command == "report":
        report(c)
        return
    if a.command in ("cleanup", "recover", "verify"):
        if not a.sample:
            ap.error("--sample required")
        {"cleanup":cleanup,"recover":recover,"verify":reconcile_archive}[a.command](c, a.sample)
        return
    while True:
        deferred = 0
        try:
            if a.command == "submit" and c["auto_cleanup"]:
                cleanup_archived(c)
            if a.command == "submit":
                deferred = release_failed_raw(c)
            if a.command == "download":
                progressed = download_one(c, a.sample)
            else:
                fetched = submit_fetch(c) if local_mode(c) else False
                progressed = submit_one(c, a.through, a.sample) or fetched
        except Exception as error:
            if not a.watch:
                raise
            print(f"Held for review: {error}", file=sys.stderr, flush=True)
            progressed = False
        if not a.watch:
            break
        if all(status(p) == "CLEANED" for p in entries(c)):
            break
        if a.command == "submit" and only_held_left(c) and not deferred:
            held = sum(status(p) in HELD for p in entries(c))
            print(f"Nothing left to run: {held} held sample(s) need review/recover; exiting watch", flush=True)
            break
        if not progressed:
            counts = {}
            for p in entries(c):
                key = status(p)
                counts[key] = counts.get(key, 0) + 1
            print(f"Waiting ({a.command}): {json.dumps(counts, sort_keys=True)}; held samples retain disk reservations", flush=True)
        time.sleep(1 if progressed else c["poll_seconds"])
