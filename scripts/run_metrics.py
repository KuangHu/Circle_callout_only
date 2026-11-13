#!/usr/bin/env python3
"""Per-sample engineering metrics for one or more run state dirs (needs sacct).

    run_metrics.py STATE_DIR [STATE_DIR ...] --out metrics.tsv [--fetch-log submit.log ...]

Per sample: status, compute job wall seconds and allocated core-hours (sacct), sum of rule
wall seconds and CPU seconds, peak sampled node-local scratch, peak_gb reserved, downloaded
bytes, archive bytes. Fetch jobs (1 CPU) are summarized from the submit logs. Totals go to
stdout as JSON. Raw call counts are deliberately not an efficiency metric here.
"""
import argparse
import collections
import datetime as dt
import glob
import json
import os
import re
import subprocess
from pathlib import Path


def sacct(ids):
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), 200):
        r = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ",".join(ids[i:i + 200]),
                            "-o", "JobID,State,Submit,Start,End,ElapsedRaw,AllocCPUS"], capture_output=True, text=True, check=True)
        for line in r.stdout.splitlines():
            f = line.split("|")
            out[f[0]] = {"state": f[1].split()[0], "submit": f[2], "start": f[3], "end": f[4],
                         "elapsed": int(f[5] or 0), "cpus": int(f[6] or 0)}
    return out


def ts(x):
    try:
        return dt.datetime.fromisoformat(x).timestamp()
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("state", nargs="+"); ap.add_argument("--out", required=True)
    ap.add_argument("--fetch-log", nargs="*", default=[])
    a = ap.parse_args()
    rows, jobs = [], []
    for sd in a.state:
        for d in sorted(Path(sd).iterdir()):
            if not (d / "status.json").is_file() or not (d / "spec.json").is_file():
                continue
            spec = json.load(open(d / "spec.json"))["sample"]
            job = json.load(open(d / "job.json")) if (d / "job.json").is_file() else {}
            down = json.load(open(d / "download.json")) if (d / "download.json").is_file() else {}
            wall = cpu = peak = 0.0
            for m in glob.glob(str(d / "metrics" / "*.json")):
                j = json.load(open(m))
                wall += j.get("wall_seconds", 0); cpu += j.get("cpu_seconds", 0)
                peak = max(peak, j.get("sampled_peak_scratch_bytes", 0))
            arch = json.load(open(d / "archive.json")) if (d / "archive.json").is_file() else {}
            rows.append({"run": Path(sd).name, "sample": d.name, "status": json.load(open(d / "status.json"))["status"],
                         "job_id": job.get("job_id", ""), "download_bytes": sum(f["size_bytes"] for f in spec["files"]),
                         "refseq_bytes": sum(r.get("size_bytes") or 0 for r in down.get("refseqs", [])),
                         "rule_wall_s": round(wall), "rule_cpu_s": round(cpu), "peak_scratch_gb": round(peak / 1e9, 2),
                         "peak_gb_reserved": spec["peak_gb"],
                         "archive_bytes": arch.get("bytes", arch.get("size_bytes", ""))})
            if job.get("job_id"):
                jobs.append(job["job_id"])
    acct = sacct(jobs) if jobs else {}
    for r in rows:
        j = acct.get(r["job_id"], {})
        r["job_state"] = j.get("state", ""); r["job_wall_s"] = j.get("elapsed", "")
        r["core_h"] = round(j["elapsed"] * j["cpus"] / 3600, 3) if j else ""
        r["queue_wait_s"] = round(ts(j["start"]) - ts(j["submit"])) if j and ts(j["start"]) and ts(j["submit"]) else ""
    fetch = []
    for log in a.fetch_log:
        fetch += re.findall(r"Fetch job (\d+)", open(log).read())
    facct = sacct(fetch) if fetch else {}
    cols = list(rows[0]) if rows else []
    with open(a.out, "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")
    done = [r for r in rows if r["status"] == "CLEANED"]
    num = lambda k, rs: [r[k] for r in rs if isinstance(r[k], (int, float))]  # noqa: E731
    med = lambda v: sorted(v)[len(v) // 2] if v else None  # noqa: E731
    total = {"samples": len(rows), "status": dict(collections.Counter(r["status"] for r in rows)),
             "compute_core_h_all": round(sum(num("core_h", rows)), 2),
             "compute_core_h_completed": round(sum(num("core_h", done)), 2),
             "core_h_per_completed_sample_median": med(num("core_h", done)),
             "job_wall_s_median": med(num("job_wall_s", done)), "job_wall_s_max": max(num("job_wall_s", done), default=None),
             "download_gb": round(sum(r["download_bytes"] for r in rows) / 1e9, 2),
             "peak_scratch_gb_max": max(num("peak_scratch_gb", rows), default=None),
             "peak_over_reserved_max": max([r["peak_scratch_gb"] / r["peak_gb_reserved"] for r in rows if r["peak_gb_reserved"]], default=None),
             "fetch_jobs": len(facct), "fetch_cpus": sorted({j["cpus"] for j in facct.values()}),
             "fetch_core_h": round(sum(j["elapsed"] * j["cpus"] for j in facct.values()) / 3600, 3),
             "queue_wait_s_median": med(num("queue_wait_s", rows))}
    print(json.dumps(total, indent=2))


if __name__ == "__main__":
    main()
