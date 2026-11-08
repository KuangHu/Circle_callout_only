#!/usr/bin/env python3
"""Seeded, stratified selection of bacterial ONT WGS BioSamples with their own assembly.

    python scripts/select_ont.py master_list.tsv ena_ont_runs.tsv --include PILOT_selection.tsv \
        --target 200 --out DIR

Cumulative de-duplicated target: every BioSample in --include (Sample_ID, ont_run; the
validated pilot) is kept with its ONT run; new BioSamples fill up to --target. New samples
do not need an Illumina run: the master_list supplies Sample_ID, Scientific_Name and the
BioSample's own assembly (Ref_FTP); the ENA snapshot supplies ONT runs (base_count,
read_count). Strata: top --top-species species (by number of eligible BioSamples) get
--per-species each, the rest form "other". Within a stratum: seeded draw, checks, then an even
pick across global tertiles of ONT bases (data volume).

Checks (failures -> skipped.tsv): ONT run in the SRA Open Data bucket (HEAD size + ETag),
object <= --max-gb; ENA metadata OXFORD_NANOPORE / WGS / GENOMIC with sample_accession ==
Sample_ID; base_count >= --min-ont-bases; Ref_FTP present.

Recorded per sample: instrument model, library name / selection / construction protocol
(kit, where submitters gave it), center, first_public, bases, reads, mean read length.
Coverage needs the assembly size: add_coverage() joins prepare_references.py output later.
Engineering/method cohort, not a prevalence sample.
"""
import argparse
import collections
import csv
import json
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sra_manifest import BUCKET, head  # noqa: E402
import select_runs  # noqa: E402

FIELDS = ("run_accession,sample_accession,instrument_platform,instrument_model,library_strategy,library_layout,"
          "library_source,library_selection,library_name,library_construction_protocol,base_count,read_count,"
          "center_name,first_public,tax_id,experiment_accession")
COLS = ["Sample_ID", "run", "Scientific_Name", "stratum", "origin", "ont_run", "ont_size", "ont_etag", "ont_model",
        "ont_bases", "ont_reads", "mean_read_len", "library_name", "library_selection", "library_construction_protocol",
        "center_name", "first_public", "experiment_accession", "Ref_FTP"]


def add_coverage(selection, references, out):
    """selection.tsv + prepare_references.py TSV (total_bases) -> coverage column."""
    size = {r["Sample_ID"]: int(r["total_bases"]) for r in csv.DictReader(open(references), delimiter="\t")
            if r["status"] == "OK" and r["total_bases"]}
    rows = list(csv.DictReader(open(selection), delimiter="\t"))
    with open(out, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]) + ["assembly_bp", "ont_coverage"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            bp = size.get(r["Sample_ID"])
            w.writerow(dict(r, assembly_bp=bp or "", ont_coverage=round(int(r["ont_bases"]) / bp, 1) if bp else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_list"); ap.add_argument("ont_runs")
    ap.add_argument("--include", help="pilot selection.tsv (Sample_ID, ont_run) kept as is")
    ap.add_argument("--out", required=True)
    ap.add_argument("--target", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261109)
    ap.add_argument("--top-species", type=int, default=15)
    ap.add_argument("--per-species", type=int, default=6)
    ap.add_argument("--pool-factor", type=int, default=3)
    ap.add_argument("--max-gb", type=float, default=3.0)
    ap.add_argument("--min-ont-bases", type=float, default=100e6)
    ap.add_argument("--add-coverage", nargs=2, metavar=("REFERENCES_TSV", "OUT_TSV"),
                    help="only join coverage onto DIR/selection.tsv and exit")
    a = ap.parse_args()
    out = Path(a.out)
    if a.add_coverage:
        add_coverage(out / "selection.tsv", *a.add_coverage); return
    rng = random.Random(a.seed)
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(a.ont_runs, out / "ena_ont_runs.snapshot.tsv")
    select_runs.FIELDS = FIELDS

    ont = collections.defaultdict(list)
    with open(a.ont_runs, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["base_count"] and float(row["base_count"]) >= a.min_ont_bases:
                ont[row["sample_accession"]].append(row)
    master = {}
    with open(a.master_list, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Sample_ID"] in ont and row["Ref_FTP"].strip() and row["Sample_ID"] not in master:
                master[row["Sample_ID"]] = row  # de-duplicate BioSamples (master_list has one row per run)
    pilot = []
    if a.include:
        for r in csv.DictReader(open(a.include), delimiter="\t"):
            pilot.append((r["Sample_ID"], r["ont_run"]))
    pilot_ids = {s for s, _ in pilot}
    by_species = collections.defaultdict(list)
    for sid, row in master.items():
        if sid not in pilot_ids:
            by_species[row["Scientific_Name"]].append(row)
    need = a.target - len(pilot)
    ranked = sorted(by_species, key=lambda s: (-len(by_species[s]), s))
    quotas = {s: a.per_species for s in ranked[:a.top_species]}
    quotas["other"] = max(0, need - sum(quotas.values()))
    strata = {s: by_species[s] for s in ranked[:a.top_species]}
    strata["other"] = [r for s in ranked[a.top_species:] for r in by_species[s]]

    skipped, candidates, limit = [], {}, a.max_gb * 1e9

    def check_run(row, run):
        try:
            size, etag = head(f"{BUCKET}/{run['run_accession']}/{run['run_accession']}")
        except LookupError as error:
            return None, f"{run['run_accession']} AWS: {error}"
        if size > limit:
            return None, f"{run['run_accession']} {size/1e9:.1f} GB > {a.max_gb}"
        return (run["run_accession"], size, etag), None

    for name, rows in strata.items():
        checked = []
        for row in rng.sample(rows, min(len(rows), quotas[name] * a.pool_factor)):
            runs = sorted(ont[row["Sample_ID"]], key=lambda r: r["run_accession"]); rng.shuffle(runs)
            pick, reasons = None, []
            for run in runs:
                pick, why = check_run(row, run)
                if pick:
                    break
                reasons.append(why)
            if not pick:
                skipped.append((row["Sample_ID"], "", name, "no usable ONT run: " + "; ".join(reasons)[:200])); continue
            checked.append(dict(row, stratum=name, origin="new", ont_run=pick[0], ont_size=pick[1], ont_etag=pick[2]))
        candidates[name] = checked
    pilot_rows = []
    for sid, run in pilot:
        row = master.get(sid)
        if row is None:
            sys.exit(f"pilot sample {sid} not eligible in master_list/ONT snapshot")
        size, etag = head(f"{BUCKET}/{run}/{run}")
        pilot_rows.append(dict(row, stratum="pilot", origin="pilot", ont_run=run, ont_size=size, ont_etag=etag))

    meta = select_runs.ena_meta([c["ont_run"] for rows in candidates.values() for c in rows] + [p["ont_run"] for p in pilot_rows])

    def problem(c):
        m = meta.get(c["ont_run"])
        if not m:
            return f"ENA: no metadata for {c['ont_run']}"
        if (m["instrument_platform"], m["library_strategy"], m["library_source"]) != ("OXFORD_NANOPORE", "WGS", "GENOMIC"):
            return f"ENA {c['ont_run']}: {m['instrument_platform']}/{m['library_strategy']}/{m['library_source']}"
        if m["sample_accession"] != c["Sample_ID"]:
            return f"ENA {c['ont_run']}: sample {m['sample_accession']}"
        return None

    def enrich(c):
        m = meta[c["ont_run"]]
        bases, reads = int(m["base_count"] or 0), int(m["read_count"] or 0)
        return dict(c, run=c["ont_run"], ont_model=m["instrument_model"], ont_bases=bases, ont_reads=reads,
                    mean_read_len=round(bases / reads) if reads else "", **{k: m.get(k, "") for k in (
                        "library_name", "library_selection", "library_construction_protocol", "center_name",
                        "first_public", "experiment_accession")})

    selected = []
    for p in pilot_rows:
        why = problem(p)
        if why:
            sys.exit(f"pilot sample {p['Sample_ID']} fails metadata check: {why}")
        selected.append(enrich(p))
    ok_all = {name: [enrich(c) for c in rows if not problem(c)] for name, rows in candidates.items()}
    for name, rows in candidates.items():
        for c in rows:
            why = problem(c)
            if why:
                skipped.append((c["Sample_ID"], c["ont_run"], name, why))
    vols = sorted(c["ont_bases"] for rows in ok_all.values() for c in rows)
    cut = [vols[len(vols) // 3], vols[2 * len(vols) // 3]] if vols else [0, 0]
    tertile = lambda b: 0 if b < cut[0] else (1 if b < cut[1] else 2)  # noqa: E731
    for name, ok in ok_all.items():
        bins = [[c for c in ok if tertile(c["ont_bases"]) == t] for t in range(3)]
        picked, t = [], 0
        while len(picked) < quotas[name] and any(bins):
            if bins[t % 3]:
                picked.append(bins[t % 3].pop(0))
            t += 1
        if len(picked) < quotas[name]:
            print(f"WARNING: stratum {name!r} has only {len(picked)}/{quotas[name]} eligible samples", file=sys.stderr)
        for c in ok:
            if c not in picked:
                skipped.append((c["Sample_ID"], c["ont_run"], name, "eligible, not needed (stratum full)"))
        selected += picked
    assert len({c["Sample_ID"] for c in selected}) == len(selected), "duplicate BioSample"
    with open(out / "selection.tsv", "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(selected)
    with open(out / "skipped.tsv", "w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["Sample_ID", "ont_run", "stratum", "reason"]); w.writerows(skipped)
    total = sum(c["ont_size"] for c in selected)
    json.dump({"seed": a.seed, "target": a.target, "pilot_included": len(pilot_rows), "new": len(selected) - len(pilot_rows),
               "top_species": a.top_species, "per_species": a.per_species, "quotas": quotas, "pool_factor": a.pool_factor,
               "max_gb": a.max_gb, "min_ont_bases": a.min_ont_bases, "ont_bases_tertile_cuts": cut,
               "eligible_biosamples": len(master), "selected": len(selected), "selected_bytes": total,
               "skipped": len(skipped), "master_list": str(Path(a.master_list).resolve()),
               "ont_runs_snapshot": "ena_ont_runs.snapshot.tsv",
               "purpose": "ONT-first circularization-junction cohort; not a prevalence sample"},
              open(out / "selection.json", "w"), indent=2)
    print(f"selected {len(selected)} BioSamples ({len(pilot_rows)} pilot + {len(selected) - len(pilot_rows)} new, "
          f"{total/1e9:.1f} GB); skipped {len(skipped)}", file=sys.stderr)


if __name__ == "__main__":
    main()
