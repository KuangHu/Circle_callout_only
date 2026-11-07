#!/usr/bin/env python3
"""Stratified, seeded selection of Illumina paired-end WGS runs from master_list.tsv.

Strata: the top-N species (by row count) get --per-species runs each; all remaining
species form one "other" stratum with --other runs. Within each stratum, candidates are
drawn at random, checked, and picked evenly across global size tertiles of their
NCBI SRA (AWS) object sizes. Each row is one BioSample with one run; the reference is
that row's own assembly (Ref_FTP), never a species-level substitute.

Checks per candidate (failures are written to skipped.tsv with the reason):
  - run present in NCBI's SRA Open Data bucket on AWS (HEAD: size + ETag)
  - ENA read_run metadata: instrument_platform=ILLUMINA, library_strategy=WGS,
    library_layout=PAIRED, sample_accession == Sample_ID
  - Ref_FTP present

Outputs in --out: selection.tsv, skipped.tsv, selection.json (parameters + seed).
For engineering trials; this is not a sample for estimating prevalence.
"""
import argparse
import collections
import csv
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sra_manifest import BUCKET, head  # noqa: E402

ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,sample_accession,instrument_platform,instrument_model,library_strategy,library_layout,base_count,read_count"


def ena_meta(runs):
    out = {}
    for i in range(0, len(runs), 200):
        data = urllib.parse.urlencode({"result": "read_run", "fields": FIELDS, "format": "tsv",
                                       "includeAccessionType": "run", "includeAccessions": ",".join(runs[i:i+200])}).encode()
        for attempt in range(4):
            try:
                req = urllib.request.Request(ENA, data=data, headers={"User-Agent": "circle-stream-select"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    rows = list(csv.DictReader(r.read().decode().splitlines(), delimiter="\t"))
                break
            except Exception as error:  # noqa: BLE001
                if attempt == 3:
                    raise RuntimeError(f"ENA metadata query failed: {error}")
                time.sleep(2 ** attempt)
        out.update({row["run_accession"]: row for row in rows})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_list")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261107)
    ap.add_argument("--top-species", type=int, default=10)
    ap.add_argument("--per-species", type=int, default=6)
    ap.add_argument("--other", type=int, default=20)
    ap.add_argument("--pool-factor", type=int, default=4, help="candidates drawn per wanted run")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)

    by_species = collections.defaultdict(list)
    with open(a.master_list, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_species[row["Scientific_Name"]].append(row)
    ranked = sorted(by_species, key=lambda s: (-len(by_species[s]), s))
    strata = {s: (by_species[s], a.per_species) for s in ranked[:a.top_species]}
    strata["other"] = ([r for s in ranked[a.top_species:] for r in by_species[s]], a.other)

    skipped, candidates = [], {}
    for name, (rows, want) in strata.items():
        pool = rng.sample(rows, min(len(rows), want * a.pool_factor))
        checked = []
        for row in pool:
            # Run accession from the R1 file name (path directories are shard prefixes like SRR286/).
            m = re.search(r"/([SED]RR\d+)_1\.fastq\.gz(?:;|$)", row["Read_FTP"])
            run = m.group(1) if m else ""
            if not run or not row["Ref_FTP"].strip():
                skipped.append((row["Sample_ID"], run, name, "missing run or Ref_FTP"))
                continue
            try:
                size, etag = head(f"{BUCKET}/{run}/{run}")
            except LookupError as error:
                skipped.append((row["Sample_ID"], run, name, f"AWS: {error}"))
                continue
            checked.append(dict(row, run=run, stratum=name, size_bytes=size, s3_etag=etag))
        candidates[name] = (checked, want)

    meta = ena_meta([c["run"] for checked, _ in candidates.values() for c in checked])
    sizes = sorted(c["size_bytes"] for checked, _ in candidates.values() for c in checked)
    cut = [sizes[len(sizes) // 3], sizes[2 * len(sizes) // 3]] if sizes else [0, 0]
    tertile = lambda b: 0 if b < cut[0] else (1 if b < cut[1] else 2)  # noqa: E731

    selected = []
    for name, (checked, want) in candidates.items():
        ok = []
        for c in checked:
            m = meta.get(c["run"])
            reason = None
            if not m:
                reason = "ENA: no read_run metadata"
            elif m["instrument_platform"] != "ILLUMINA":
                reason = f"ENA platform {m['instrument_platform']}"
            elif m["library_strategy"] != "WGS":
                reason = f"ENA strategy {m['library_strategy']}"
            elif m["library_layout"] != "PAIRED":
                reason = f"ENA layout {m['library_layout']}"
            elif m["sample_accession"] != c["Sample_ID"]:
                reason = f"ENA sample {m['sample_accession']} != {c['Sample_ID']}"
            if reason:
                skipped.append((c["Sample_ID"], c["run"], name, reason))
                continue
            ok.append(dict(c, instrument_model=m["instrument_model"], base_count=m["base_count"]))
        # Even pick across size tertiles (round-robin), random order within each tertile.
        bins = [[c for c in ok if tertile(c["size_bytes"]) == t] for t in range(3)]
        picked, t = [], 0
        while len(picked) < want and any(bins):
            if bins[t % 3]:
                picked.append(bins[t % 3].pop(0))
            t += 1
        if len(picked) < want:
            print(f"WARNING: stratum {name!r} has only {len(picked)}/{want} eligible runs", file=sys.stderr)
        for c in ok:
            if c not in picked:
                skipped.append((c["Sample_ID"], c["run"], name, "eligible, not needed (stratum full)"))
        selected += picked

    cols = ["Sample_ID", "run", "Scientific_Name", "stratum", "size_bytes", "s3_etag", "instrument_model", "base_count", "Ref_FTP"]
    with open(out / "selection.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows([[c[k] for k in cols] for c in selected])
    with open(out / "skipped.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["Sample_ID", "run", "stratum", "reason"])
        w.writerows(skipped)
    json.dump({"seed": a.seed, "top_species": a.top_species, "per_species": a.per_species, "other": a.other,
               "pool_factor": a.pool_factor, "size_tertile_cuts_bytes": cut, "selected": len(selected),
               "selected_bytes": sum(c["size_bytes"] for c in selected), "skipped": len(skipped),
               "master_list": str(Path(a.master_list).resolve()),
               "purpose": "engineering trial; not a prevalence sample"},
              open(out / "selection.json", "w"), indent=2)
    print(f"selected {len(selected)} runs ({sum(c['size_bytes'] for c in selected)/1e9:.1f} GB); skipped {len(skipped)}", file=sys.stderr)


if __name__ == "__main__":
    main()
