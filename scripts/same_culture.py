#!/usr/bin/env python3
"""Data check that ONT and Illumina runs of one BioSample come from the same culture / DNA.

    same_culture.py --ont DIR [DIR ...] --illumina-depth DIR --out same_culture.tsv

Both platforms are aligned to the same per-BioSample assembly. Per contig depth is normalized
to the median depth of contigs >= 100 kb (the chromosome scale), so the normalized value is a
replicon copy number. A different culture or DNA prep typically shows as a replicon (>= 5 kb)
present in one platform (>= 0.25x) and absent from the other (< 0.05x), or as plasmid copy
numbers that disagree strongly. Metadata (same study / centre) cannot establish this.

Per pair: contigs compared, replicons (contigs >= 5 kb) present in one platform only, the
fraction of each platform's sequence that aligns, and a verdict:
  genome_mismatch          ONT aligned bases < 50% of the run (--selection gives run bases),
                           Illumina mapped reads < 90%, a contig >= 100 kb present in one
                           platform only, or one-platform-only contigs >= 5% of the assembly:
                           not the same genome/culture; excluded from cross-platform matching
  replicon_differs         a contig >= 20 kb present in one platform only (plasmid gain/loss)
  small_replicon_differs   only contigs < 20 kb differ (often ONT size selection losing small
                           plasmids rather than a different culture)
  consistent               none of the above
Spearman of log copy number and max copy-number fold are informational only: in fragmented
assemblies most contigs are ~1x chromosome pieces and their ranking is noise.
ONT depth comes from evidence/cresil/contig_depth.tsv inside <SID>_ont.tar.gz; Illumina from
<SID>.tsv (contig, length, aligned_bases, mean_depth; trailing #reads_seen / #reads_mapped).
A low Illumina mapped fraction to this assembly is reported too (contamination or a
different strain).
"""
import argparse
import csv
import io
import math
import statistics
import tarfile
from pathlib import Path


def norm(rows):
    big = [r[1] for r in rows.values() if r[0] >= 100000] or [r[1] for r in rows.values()]
    m = statistics.median(big) if big else 0
    return {c: (L, d / m if m else 0.0) for c, (L, d) in rows.items()}


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(x, y):
    if len(x) < 3:
        return None
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ont", nargs="+", required=True); ap.add_argument("--illumina-depth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--selection", help="select_ont.py table with Sample_ID and ont_bases")
    a = ap.parse_args()
    run_bases = {}
    if a.selection:
        run_bases = {r["Sample_ID"]: int(r["ont_bases"] or 0) for r in csv.DictReader(open(a.selection), delimiter="\t")}
    ont = {}
    for d in a.ont:
        for p in sorted(Path(d).glob("*_ont.tar.gz")):
            ont[p.name[:-len("_ont.tar.gz")]] = p
    out = []
    for ip in sorted(Path(a.illumina_depth).glob("*.tsv")):
        sid = ip.stem
        if not open(ip).readline().startswith("contig\t"):
            continue  # not a depth table (e.g. pairs.tsv)
        if sid not in ont:
            out.append({"sample": sid, "verdict": "no ONT archive"}); continue
        with tarfile.open(ont[sid]) as tar:
            m = next((m for m in tar if m.name.endswith("evidence/cresil/contig_depth.tsv")), None)
            text = tar.extractfile(m).read().decode() if m else ""
        o = {r["contig"]: (int(r["length"]), float(r["mean_depth"])) for r in csv.DictReader(io.StringIO(text), delimiter="\t")}
        lines = open(ip).read().splitlines()
        meta = dict(l[1:].split("\t", 1) for l in lines if l.startswith("#"))
        i = {r["contig"]: (int(r["length"]), float(r["mean_depth"]))
             for r in csv.DictReader([l for l in lines if not l.startswith("#")], delimiter="\t")}
        if not o or not i:
            out.append({"sample": sid, "verdict": "missing depth"}); continue
        on, il = norm(o), norm(i)
        common = [c for c in on if c in il]
        big = [c for c in common if on[c][0] >= 5000]
        only_ont = [c for c in big if on[c][1] >= 0.25 and il[c][1] < 0.05]
        only_ill = [c for c in big if il[c][1] >= 0.25 and on[c][1] < 0.05]
        both = [c for c in big if on[c][1] >= 0.25 and il[c][1] >= 0.25]
        rho = spearman([math.log(on[c][1] + 1e-3) for c in big], [math.log(il[c][1] + 1e-3) for c in big])
        fold = max((max(on[c][1], il[c][1]) / min(on[c][1], il[c][1]) for c in both), default=1.0)
        diff = only_ont + only_ill
        asm = sum(L for L, _ in on.values())
        ont_frac = sum(d * L for L, d in o.values()) / run_bases[sid] if run_bases.get(sid) else None
        ill_frac = int(meta["reads_mapped"]) / int(meta["reads_seen"]) if meta.get("reads_seen") else None
        if ((ont_frac is not None and ont_frac < 0.5) or (ill_frac is not None and ill_frac < 0.9)
                or any(on[c][0] >= 100000 for c in diff) or sum(on[c][0] for c in diff) >= 0.05 * asm):
            verdict = "genome_mismatch"
        elif any(on[c][0] >= 20000 for c in diff):
            verdict = "replicon_differs"
        elif diff:
            verdict = "small_replicon_differs"
        else:
            verdict = "consistent"
        out.append({"sample": sid, "contigs": len(common), "contigs_ge5kb": len(big),
                    "ont_median_depth_100kb": round(statistics.median([d for L, d in o.values() if L >= 100000] or [0]), 1),
                    "ill_median_depth_100kb": round(statistics.median([d for L, d in i.values() if L >= 100000] or [0]), 1),
                    "ont_aligned_base_fraction": "" if ont_frac is None else round(ont_frac, 3),
                    "ill_mapped_fraction": "" if ill_frac is None else round(ill_frac, 4),
                    "spearman_ge5kb": "" if rho is None else round(rho, 3), "max_copy_fold_ge5kb": round(fold, 2),
                    "ont_only_replicons": ";".join(f"{c}:{on[c][0]}:{on[c][1]:.2f}x" for c in only_ont),
                    "illumina_only_replicons": ";".join(f"{c}:{il[c][0]}:{il[c][1]:.2f}x" for c in only_ill),
                    "verdict": verdict})
    cols = list(dict.fromkeys(k for r in out for k in r))
    with open(a.out, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n", restval=""); w.writeheader(); w.writerows(out)
    print({v: sum(r["verdict"] == v for r in out) for v in sorted({r["verdict"] for r in out})})


if __name__ == "__main__":
    main()
