#!/usr/bin/env python3
"""Seeded, stratified selection of BioSamples with BOTH an Illumina paired-end run and
an Oxford Nanopore run, plus the sample's own assembly, for a same-sample comparison.

    python scripts/select_paired.py master_list.tsv ena_ont_runs.tsv --out DIR

ena_ont_runs.tsv is a saved ENA portal snapshot (result=read_run, ONT, WGS, GENOMIC,
bacteria) with run_accession, sample_accession, base_count; it is copied into DIR.
Strata: top-N species by number of paired samples get --per-species each; the rest
form one "other" stratum. Within a stratum: random draw, checks, then an even pick
across global size tertiles of the ONT object. Per sample one ONT run (seeded random
among that sample's runs within the size limits) and the master_list Illumina run.

Checks (failures -> skipped.tsv with reason): both runs in the SRA Open Data bucket
on AWS (HEAD: size + ETag); ENA metadata ILLUMINA/WGS/PAIRED and
OXFORD_NANOPORE/WGS with sample_accession == Sample_ID; SRA object sizes within
--max-gb; ONT runs with >= --min-ont-bases (ENA base_count); Ref_FTP present. Engineering/method pilot, not a prevalence sample.
"""
import argparse
import collections
import csv
import json
import random
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sra_manifest import BUCKET, head  # noqa: E402
from select_runs import ena_meta  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_list"); ap.add_argument("ont_runs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--top-species", type=int, default=10)
    ap.add_argument("--per-species", type=int, default=3)
    ap.add_argument("--other", type=int, default=10)
    ap.add_argument("--pool-factor", type=int, default=4)
    ap.add_argument("--max-gb", type=float, default=3.0, help="max SRA object size per run")
    ap.add_argument("--min-ont-bases", type=float, default=100e6, help="min ONT bases (ENA base_count), ~20x of 5 Mb")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(a.ont_runs, out / "ena_ont_runs.snapshot.tsv")

    ont = collections.defaultdict(list)
    with open(a.ont_runs, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["base_count"] and float(row["base_count"]) >= a.min_ont_bases:
                ont[row["sample_accession"]].append(row["run_accession"])
    by_species = collections.defaultdict(list)
    with open(a.master_list, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Sample_ID"] in ont:
                by_species[row["Scientific_Name"]].append(row)
    ranked = sorted(by_species, key=lambda s: (-len(by_species[s]), s))
    strata = {s: (by_species[s], a.per_species) for s in ranked[:a.top_species]}
    strata["other"] = ([r for s in ranked[a.top_species:] for r in by_species[s]], a.other)

    skipped, candidates, limit = [], {}, a.max_gb * 1e9
    for name, (rows, want) in strata.items():
        checked = []
        for row in rng.sample(rows, min(len(rows), want * a.pool_factor)):
            sid = row["Sample_ID"]
            m = re.search(r"/([SED]RR\d+)_1\.fastq\.gz(?:;|$)", row["Read_FTP"])
            ill = m.group(1) if m else ""
            onts = sorted(ont[sid]); rng.shuffle(onts)
            if not ill or not row["Ref_FTP"].strip():
                skipped.append((sid, ill, "", name, "missing Illumina run or Ref_FTP")); continue
            try:
                ill_size, ill_etag = head(f"{BUCKET}/{ill}/{ill}")
            except LookupError as error:
                skipped.append((sid, ill, "", name, f"Illumina AWS: {error}")); continue
            pick, reasons = None, []
            for run in onts:
                try:
                    size, etag = head(f"{BUCKET}/{run}/{run}")
                except LookupError as error:
                    reasons.append(f"{run} AWS: {error}"); continue
                if size > limit:
                    reasons.append(f"{run} {size/1e9:.1f} GB > {a.max_gb}"); continue
                pick = (run, size, etag); break
            if not pick:
                skipped.append((sid, ill, ";".join(onts), name, "no usable ONT run: " + "; ".join(reasons)[:200])); continue
            if ill_size > limit:
                skipped.append((sid, ill, pick[0], name, f"Illumina {ill_size/1e9:.1f} GB > {a.max_gb}")); continue
            checked.append(dict(row, stratum=name, ill_run=ill, ill_size=ill_size, ill_etag=ill_etag,
                                ont_run=pick[0], ont_size=pick[1], ont_etag=pick[2]))
        candidates[name] = (checked, want)

    meta = ena_meta([c[k] for checked, _ in candidates.values() for c in checked for k in ("ill_run", "ont_run")])
    sizes = sorted(c["ont_size"] for checked, _ in candidates.values() for c in checked)
    cut = [sizes[len(sizes) // 3], sizes[2 * len(sizes) // 3]] if sizes else [0, 0]
    tertile = lambda b: 0 if b < cut[0] else (1 if b < cut[1] else 2)  # noqa: E731

    def problem(c):
        for key, platform, layout in (("ill_run", "ILLUMINA", "PAIRED"), ("ont_run", "OXFORD_NANOPORE", None)):
            m = meta.get(c[key])
            if not m:
                return f"ENA: no metadata for {c[key]}"
            if m["instrument_platform"] != platform or m["library_strategy"] != "WGS":
                return f"ENA {c[key]}: {m['instrument_platform']}/{m['library_strategy']}"
            if layout and m["library_layout"] != layout:
                return f"ENA {c[key]}: layout {m['library_layout']}"
            if m["sample_accession"] != c["Sample_ID"]:
                return f"ENA {c[key]}: sample {m['sample_accession']}"
        return None

    selected = []
    for name, (checked, want) in candidates.items():
        ok = []
        for c in checked:
            reason = problem(c)
            if reason:
                skipped.append((c["Sample_ID"], c["ill_run"], c["ont_run"], name, reason))
            else:
                ok.append(dict(c, ill_model=meta[c["ill_run"]]["instrument_model"], ont_model=meta[c["ont_run"]]["instrument_model"],
                               ont_bases=meta[c["ont_run"]]["base_count"]))
        bins = [[c for c in ok if tertile(c["ont_size"]) == t] for t in range(3)]
        picked, t = [], 0
        while len(picked) < want and any(bins):
            if bins[t % 3]:
                picked.append(bins[t % 3].pop(0))
            t += 1
        if len(picked) < want:
            print(f"WARNING: stratum {name!r} has only {len(picked)}/{want} eligible samples", file=sys.stderr)
        for c in ok:
            if c not in picked:
                skipped.append((c["Sample_ID"], c["ill_run"], c["ont_run"], name, "eligible, not needed (stratum full)"))
        selected += picked

    cols = ["Sample_ID", "run", "Scientific_Name", "stratum", "ill_run", "ill_size", "ill_etag", "ill_model",
            "ont_run", "ont_size", "ont_etag", "ont_model", "ont_bases", "Ref_FTP"]
    with open(out / "selection.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n"); w.writerow(cols)
        w.writerows([[c["ont_run"] if k == "run" else c[k] for k in cols] for c in selected])
    with open(out / "skipped.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["Sample_ID", "ill_run", "ont_runs", "stratum", "reason"]); w.writerows(skipped)
    total = sum(c["ill_size"] + c["ont_size"] for c in selected)
    json.dump({"seed": a.seed, "top_species": a.top_species, "per_species": a.per_species, "other": a.other,
               "pool_factor": a.pool_factor, "max_gb": a.max_gb, "min_ont_bases": a.min_ont_bases, "ont_size_tertile_cuts_bytes": cut,
               "selected": len(selected), "selected_bytes": total, "skipped": len(skipped),
               "paired_biosamples_available": sum(len(v) for v in by_species.values()),
               "master_list": str(Path(a.master_list).resolve()), "ont_runs_snapshot": "ena_ont_runs.snapshot.tsv",
               "purpose": "same-sample Illumina vs nanopore method pilot; not a prevalence sample"},
              open(out / "selection.json", "w"), indent=2)
    print(f"selected {len(selected)} samples ({total/1e9:.1f} GB, both platforms); skipped {len(skipped)}", file=sys.stderr)


if __name__ == "__main__":
    main()
