#!/usr/bin/env python3
"""Download, check and permanently store the per-sample references of a selection.

Run as ONE 1-CPU Slurm job for a whole selection (not one job per reference):

    python scripts/prepare_references.py selection.tsv --store /path/to/reference_store --out refs.tsv

For each row (Sample_ID, run, Ref_FTP): download over HTTPS, decompress, check FASTA
structure (headers, unique IDs, IUPAC bases, non-empty records), and store as
<store>/<sha256-of-FASTA>.fa, read-only (shared with circlestream's reference_store).
Sequence IDs are kept as-is when their first word is shell/BED-safe; otherwise the
header is reduced to a sanitized first word and the original header map is saved
next to the store entry. Failures are recorded in the output with the reason; a
failed reference means that sample is left out of the manifest, never paired with
another assembly.
"""
import argparse
import csv
import gzip
import hashlib
import io
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

SAFE_ID = re.compile(r"[A-Za-z0-9._-]+")
IUPAC = set(b"ACGTNRYKMSWBDHVacgtnrykmswbdhv")


def fetch(url, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "circle-stream-references"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as error:  # noqa: BLE001
            if attempt + 1 == retries:
                raise RuntimeError(f"download failed: {error}")
            time.sleep(2 ** attempt)


def check_fasta(data):
    """Return (normalized_bytes, header_map or None, n_records, total_bases)."""
    out, ids, header_map, n, total, seq_len = io.BytesIO(), set(), {}, 0, 0, None
    renamed = False
    for raw in data.splitlines():
        line = raw.rstrip(b"\r")
        if not line:
            continue
        if line.startswith(b">"):
            if seq_len == 0:
                raise ValueError("empty sequence record")
            header = line[1:].decode("ascii", "replace")
            first = header.split()[0] if header.split() else ""
            if not first:
                raise ValueError("empty FASTA header")
            new = first if SAFE_ID.fullmatch(first) else re.sub(r"[^A-Za-z0-9._-]", "_", first)
            if new != first:
                renamed = True
            if new in ids:
                raise ValueError(f"duplicate sequence ID {new}")
            ids.add(new)
            header_map[new] = header
            out.write(b">" + (new.encode() if new != first else line[1:]) + b"\n")
            n += 1
            seq_len = 0
        else:
            if n == 0:
                raise ValueError("sequence before first header")
            if not set(line) <= IUPAC:
                raise ValueError("non-IUPAC characters in sequence")
            seq_len += len(line)
            total += len(line)
            out.write(line + b"\n")
    if n == 0 or seq_len == 0:
        raise ValueError("no sequences or empty last record")
    return out.getvalue(), (header_map if renamed else None), n, total


def store(data, store_dir):
    sha = hashlib.sha256(data).hexdigest()
    dest = store_dir / f"{sha}.fa"
    if dest.exists():
        if hashlib.sha256(dest.read_bytes()).hexdigest() != sha:
            raise ValueError(f"store entry corrupted: {dest}")
        return sha, dest
    tmp = store_dir / f".{sha}.{os.getpid()}.tmp"
    tmp.write_bytes(data)
    os.chmod(tmp, 0o444)
    os.replace(tmp, dest)
    return sha, dest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("selection")
    ap.add_argument("--store", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    store_dir = Path(a.store)
    store_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(a.selection, newline=""), delimiter="\t"))
    cols = ["Sample_ID", "run", "status", "reason", "ref_url", "source_gz_md5", "sha256", "store_path",
            "records", "total_bases", "renamed_headers", "seconds"]
    ok = 0
    with open(a.out, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            url = "https://" + row["Ref_FTP"].strip().removeprefix("ftp://").removeprefix("https://")
            rec = {"Sample_ID": row["Sample_ID"], "run": row["run"], "ref_url": url}
            start = time.time()
            try:
                gz = fetch(url)
                rec["source_gz_md5"] = hashlib.md5(gz).hexdigest()
                data = gzip.decompress(gz) if url.endswith(".gz") else gz
                norm, header_map, n, total = check_fasta(data)
                if total < 100_000:
                    raise ValueError(f"assembly too small for a bacterial genome ({total} bp)")
                sha, dest = store(norm, store_dir)
                if header_map:
                    (store_dir / f"{sha}.headers.tsv").write_text(
                        "".join(f"{k}\t{v}\n" for k, v in header_map.items()), encoding="utf-8")
                rec.update(status="OK", sha256=sha, store_path=str(dest), records=n, total_bases=total,
                           renamed_headers=bool(header_map))
                ok += 1
            except Exception as error:  # noqa: BLE001
                rec.update(status="FAILED", reason=str(error)[:300])
            rec["seconds"] = round(time.time() - start, 2)
            w.writerow(rec)
            handle.flush()
            print(f"{rec['status']}\t{row['Sample_ID']}\t{rec.get('reason','')}", flush=True)
    print(f"{ok}/{len(rows)} references stored in {store_dir}", file=sys.stderr)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
