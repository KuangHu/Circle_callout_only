#!/usr/bin/env python3
"""Build manifest entries for runs hosted in NCBI's SRA Open Data bucket on AWS.

Input TSV columns (header required): sample_id, accession, platform, reference
Optional column: peak_gb. Without it, peak_gb is estimated from the SRA object size
using factors measured on LRC (batch1 peaks plus the uncompressed fasterq-dump output).
Any other non-empty column (e.g. Scientific_Name, ont_model, library_construction_protocol)
is carried into the entry as "metadata" and so archived with the sample (metadata/spec.json).

    python scripts/sra_manifest.py runs.tsv > samples.batchN.json

Only HEAD requests are made (size and S3 ETag); nothing is downloaded.
"""
import argparse
import csv
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request

BUCKET = "https://sra-pub-run-odp.s3.amazonaws.com/sra"
# Nanopore: node-local peak / SRA object size over 1,039 completed 1 TB-run samples was median
# 6.2x, p99 18x, max 20.6x; 10 of 1,440 failed at 20x (up to ~23x, in the call step).
# ont1tb3 at 30x (1,447 samples): completed median 6.5x, p99 22.9x, max 28.8x; 5 exceeded it at
# 31-32x (call step). 40x: peak_gb is node-local /tmp (1.8 TB per lr7 node), not shared scratch.
PEAK_FACTOR = {"illumina": 30, "nanopore": 40}
PEAK_MIN_GB = 2
COLUMNS = {"sample_id", "accession", "platform", "reference", "peak_gb"}


def head(url, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "circle-stream-manifest"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return int(r.headers["Content-Length"]), r.headers["ETag"].strip('"')
        except urllib.error.HTTPError as error:
            if error.code in (403, 404):
                raise LookupError(f"not in SRA ODP bucket (HTTP {error.code})")
            last = error
        except (urllib.error.URLError, TimeoutError) as error:
            last = error
        time.sleep(2 ** attempt)
    raise RuntimeError(f"HEAD failed: {last}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    a = ap.parse_args()
    out, missing = [], []
    with open(a.tsv, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc, platform = row["accession"].strip(), row["platform"].strip()
            if not re.fullmatch(r"[SED]RR[0-9]{6,}", acc):
                raise ValueError(f"Not a run accession: {acc!r}")
            if platform not in PEAK_FACTOR:
                raise ValueError(f"platform must be illumina or nanopore: {platform!r}")
            url = f"{BUCKET}/{acc}/{acc}"
            try:
                size, etag = head(url)
            except LookupError as error:
                missing.append(f"{acc}\t{error}")
                continue
            peak = float(row["peak_gb"]) if row.get("peak_gb") else max(PEAK_MIN_GB, math.ceil(size / 1e9 * PEAK_FACTOR[platform] + 1))
            entry = {"sample_id": row["sample_id"].strip(), "accession": acc, "platform": platform,
                     "library": "wgs", "reference": row["reference"].strip(), "peak_gb": peak,
                     "files": [{"name": f"{acc}.sra", "url": url, "size_bytes": size, "s3_etag": etag}]}
            meta = {k: v.strip() for k, v in row.items() if k and k not in COLUMNS and v and v.strip()}
            if meta:
                entry["metadata"] = meta
            out.append(entry)
    json.dump(out, sys.stdout, indent=2)
    sys.stdout.write("\n")
    for line in missing:
        print(f"SKIPPED {line}", file=sys.stderr)
    total = sum(s["files"][0]["size_bytes"] for s in out)
    print(f"{len(out)} runs, {total/1e9:.2f} GB of SRA objects; {len(missing)} skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
