#!/usr/bin/env python3
"""Trace each junction-counter candidate through CReSIL identify_wgls -m linkage filters.

    cresil_trace.py TRIM.txt REF.fa.fai JUNCTIONS.tsv CRESIL_OUT_DIR OUT.tsv
                    [--mwindow 100000] [--minrsize 200] [--depth 5]

Re-implements the region steps of CReSIL 1.x linkage mode (identify_wgls.py) on CReSIL's
own trim.txt, to record WHERE a junction is lost, not to replace CReSIL:
  1 trim_junction_reads  reads in trim.txt (cresil trim keeps mapq >= 30 primary hits) with
                         segments within --end-tol of both junction ends
  2 region_length        breakpoints = start of + / end of - segments with order > 0 (every
                         segment after a read's first); regions = intervals between
                         breakpoints; CReSIL keeps only regions <= -mwindow (default 100 kb)
  3 read_covers_90pct    a region is kept only if one trimmed read segment covers >= 90% of it
  4 depth_minrsize       mean depth of all trimmed segments >= -depth and length >= -minrsize
  5 in_final_graph       an eccDNA_final.txt entry spans the junction's unit: the entry's
                         regions on this contig, coalesced across gaps <= --end-tol, give a
                         span whose ends agree with the junction (see span_matches). Entries
                         that only overlap the unit (e.g. a 1 kb cyclic sub-unit inside a 3 kb
                         plasmid) do not pass; they are listed in cresil_partial_ids
  6 cyclic               a spanning entry has ctc == True
excluded_at = first step that fails (steps 2-4 fail if any region of the unit fails, since
the cyclic graph then cannot close). Junctions: all whole_contig rows plus rows with
support_ge3 == yes.
"""
import argparse
import bisect
import collections
import csv
import re


def final_spans(path, gap=100):
    """eccDNA_final.txt -> [(contig, start0, end, id, cyclic)], one span per entry and contig
    after coalescing that entry's regions across gaps <= gap (CReSIL often splits one unit)."""
    spans = []
    try:
        rows = list(csv.DictReader(open(path), delimiter="\t"))
    except FileNotFoundError:
        return spans
    for r in rows:
        regs = collections.defaultdict(list)
        for reg in r["merge_region"].split(","):
            m = re.match(r"(.+):(\d+)-(\d+)_[+-]$", reg)
            regs[m.group(1)].append((int(m.group(2)) - 1, int(m.group(3))))
        for c, iv in regs.items():
            iv.sort()
            cur = list(iv[0])
            for x, y in iv[1:]:
                if x <= cur[1] + gap:
                    cur[1] = max(cur[1], y)
                else:
                    spans.append((c, cur[0], cur[1], r["id"], r["ctc"] == "True")); cur = [x, y]
            spans.append((c, cur[0], cur[1], r["id"], r["ctc"] == "True"))
    return spans


def span_matches(contig, start, end, mh, span, tol=50, frac=0.02):
    """Same rule as compare_platforms: ends within max(tol, frac x length), or shifted
    together along the junction's microhomology. ont_junctions.py stores the LEFT-most
    equivalent breakpoint, so equivalent positions lie 0..mh bp to the right only."""
    if span[0] != contig:
        return False
    ds, de = span[1] - start, span[2] - end
    lim = max(tol, frac * max(end - start, span[2] - span[1]))
    if abs(ds) <= lim and abs(de) <= lim:
        return True
    return abs(ds - de) <= tol and -tol <= ds <= mh + tol


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trim"); ap.add_argument("fai"); ap.add_argument("junctions"); ap.add_argument("cresil_out"); ap.add_argument("out")
    ap.add_argument("--mwindow", type=int, default=100000)
    ap.add_argument("--minrsize", type=int, default=200)
    ap.add_argument("--depth", type=float, default=5)
    ap.add_argument("--end-tol", type=int, default=100)
    a = ap.parse_args()
    lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) if l.strip()}
    targets = [r for r in csv.DictReader(open(a.junctions), delimiter="\t")
               if r["location"] == "whole_contig" or r["support_ge3"] == "yes"]
    contigs = {t["contig"] for t in targets}
    segs = collections.defaultdict(list)       # all trimmed segments (CReSIL aln_reads)
    allsegs = collections.defaultdict(list)    # with strand and read ID, for junction reads
    sup = collections.defaultdict(list)        # order > 0 (CReSIL trim_sup: breakpoints)
    with open(a.trim) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["ref"] not in contigs:
                continue
            s, e = int(r["r_start"]), int(r["r_end"])
            segs[r["ref"]].append((s, e))
            allsegs[r["ref"]].append((s, e, int(r["strand"]), r["readid"]))
            if int(r["order"]) > 0:
                sup[r["ref"]].append((s, e, int(r["strand"]), r["readid"]))
    final = final_spans(f"{a.cresil_out}/eccDNA_final.txt", a.end_tol)

    cov, bps, supsorted = {}, {}, {}
    for c in contigs:
        diff = [0] * (lengths[c] + 1)
        for s, e in segs[c]:
            diff[s] += 1; diff[min(e, lengths[c])] -= 1
        run, pref = 0, [0]
        for i in range(lengths[c]):
            run += diff[i]; pref.append(pref[-1] + run)
        cov[c] = pref
        # CReSIL: 5 bp breakpoint intervals ([s, s+5) for +, [e-5, e) for -), merged.
        iv = sorted({(s, s + 5) if strand == 1 else (e - 5, e) for s, e, strand, _ in sup[c]})
        merged = []
        for x, y in iv:
            if merged and x <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], y)
            else:
                merged.append([x, y])
        bps[c] = merged
        supsorted[c] = sorted((s, e) for s, e, _, _ in sup[c])

    rows = []
    for t in targets:
        c, start, end = t["contig"], int(t["start"]), int(t["end"])
        reads_s = {rid for s, e, _, rid in allsegs[c] if abs(s - start) <= a.end_tol}
        reads_e = {rid for s, e, _, rid in allsegs[c] if abs(e - end) <= a.end_tol}
        if t["location"] == "whole_contig":  # closure: segments near both contig ends
            reads_s = {rid for s, e, _, rid in allsegs[c] if s <= a.end_tol}
            reads_e = {rid for s, e, _, rid in allsegs[c] if e >= lengths[c] - a.end_tol}
        jreads = len(reads_s & reads_e)
        # Regions = gaps between merged breakpoint intervals inside the unit. Gaps shorter than
        # -minrsize are dropped by CReSIL without breaking a graph, so only longer ones decide.
        pieces, pos = [], start
        for x, y in bps[c]:
            if y <= start or x >= end:
                continue
            if x > pos:
                pieces.append((pos, x))
            pos = max(pos, y)
        if pos < end:
            pieces.append((pos, end))
        tiny = [p for p in pieces if p[1] - p[0] < a.minrsize]
        pieces = [p for p in pieces if p not in tiny]
        too_long = [p for p in pieces if p[1] - p[0] > a.mwindow]
        ss = supsorted[c]
        def covered90(p):
            lo = bisect.bisect_left(ss, (p[0] - a.mwindow * 2, 0))
            for s, e in ss[lo:]:
                if s >= p[1]:
                    break
                if min(e, p[1]) - max(s, p[0]) >= 0.9 * (p[1] - p[0]):
                    return True
            return False
        ok_len = [p for p in pieces if p not in too_long]
        not90 = [p for p in ok_len if not covered90(p)]
        mean_depth = lambda p: (cov[c][p[1]] - cov[c][p[0]]) / max(1, p[1] - p[0])  # noqa: E731
        lowdepth = [p for p in ok_len if p not in not90 and mean_depth(p) < a.depth]
        mh = int(t.get("microhomology_bp") or 0)
        hits = [f for f in final if span_matches(c, start, end, mh, f)]
        partial = [f for f in final if f not in hits and f[0] == c and min(f[2], end) - max(f[1], start) > 0]
        steps = [("trim_junction_reads", jreads > 0), ("region_length", not too_long), ("read_covers_90pct", not not90),
                 ("depth_minrsize", bool(pieces) and not lowdepth), ("in_final_graph", bool(hits)), ("cyclic", any(h[4] for h in hits))]
        excluded = next((name for name, ok in steps if not ok), "kept_cyclic")
        rows.append({"junction_id": t["junction_id"], "contig": c, "start": start, "end": end, "length": end - start,
                     "location": t["location"], "label": t["label"], "independent_molecules": t["independent_molecules"],
                     "trim_junction_reads": jreads, "regions": len(pieces), "gaps_below_minrsize": len(tiny),
                     "regions_over_mwindow": len(too_long), "max_region_bp": max((p[1] - p[0] for p in pieces), default=0),
                     "regions_not_90pct_covered": len(not90), "regions_low_depth_or_short": len(lowdepth),
                     "cresil_ids": ";".join(sorted({h[3] + ("(cyclic)" if h[4] else "") for h in hits})) or "none",
                     "cresil_partial_ids": ";".join(sorted({f"{h[3]}{'(cyclic)' if h[4] else ''}:{h[1]}-{h[2]}" for h in partial})) or "none",
                     "excluded_at": excluded})
    cols = ["junction_id", "contig", "start", "end", "length", "location", "label", "independent_molecules",
            "trim_junction_reads", "regions", "gaps_below_minrsize", "regions_over_mwindow", "max_region_bp", "regions_not_90pct_covered",
            "regions_low_depth_or_short", "cresil_ids", "cresil_partial_ids", "excluded_at"]
    with open(a.out, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
