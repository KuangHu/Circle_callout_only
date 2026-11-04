#!/usr/bin/env python3
"""Profile a rule, supervise its process group, then publish a completion marker."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from circlestream.common import GB, REPO, read_json, state_paths, tree_size, write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--rule", required=True)
    ap.add_argument("--inner", action="store_true")
    a = ap.parse_args()
    if a.inner:
        from circlestream.rules import main as inner
        sys.argv = [sys.argv[0], "--spec", a.spec, "--rule", a.rule]
        inner()
        return
    if os.name != "posix" or not Path("/usr/bin/time").exists():
        raise RuntimeError("Real rules require Linux and GNU /usr/bin/time")
    spec = read_json(a.spec)
    c, s = spec["config"], spec["sample"]
    state, scratch = state_paths(c, s["sample_id"])
    metrics = state / "metrics"
    metrics.mkdir(exist_ok=True)
    attempt = f"{a.rule}.{time.time_ns()}"
    timefile = metrics / (attempt + ".time.tsv")
    tmp = scratch / "work" / "tmp"
    tmp.mkdir(exist_ok=True)
    env = dict(os.environ, TMPDIR=str(tmp), TEMP=str(tmp), TMP=str(tmp),
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", MPLBACKEND="Agg")
    cmd = ["/usr/bin/time", "-f", "%e\t%U\t%S\t%M", "-o", str(timefile), sys.executable,
           str(REPO / "run_rule.py"), "--inner", "--spec", a.spec, "--rule", a.rule]
    started = time.monotonic()
    peak = tree_size(scratch)
    reason = None
    proc = subprocess.Popen(cmd, env=env, start_new_session=True)
    try:
        while proc.poll() is None:
            used = tree_size(scratch)
            peak = max(peak, used)
            if used > s["peak_gb"] * GB:
                reason = "scratch reservation exceeded; increase peak_gb in a new run or inspect intermediates"
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                break
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
    except BaseException:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGTERM)
        raise
    finally:
        rc = proc.wait()
        peak = max(peak, tree_size(scratch))
        stat = {"rule": a.rule, "wall_seconds": time.monotonic()-started, "exit_code": rc,
                "sampled_peak_scratch_bytes": peak, "reservation_exceeded": bool(reason),
                "allocated_cpus": c["slurm"]["cpus"], "job_id": os.environ.get("SLURM_JOB_ID")}
        if timefile.exists():
            try:
                elapsed,user,system,rss = timefile.read_text().strip().splitlines()[-1].split("\t")
                stat.update(cpu_seconds=float(user)+float(system), max_rss_kb=int(rss))
            except (ValueError, IndexError):
                stat["measurement_incomplete"] = True
        write_json(metrics / (attempt + ".json"), stat)
    if rc or reason:
        raise RuntimeError(reason or f"{a.rule} exited {rc}; inputs retained")
    write_json(scratch / "work" / (a.rule + ".ok"), stat)


if __name__ == "__main__":
    main()
