#!/usr/bin/env python3
"""Count circle junctions directly from long-read alignments (PAF from minimap2 -c).

The long-read counterpart of Circle_finder's split reads, independent of CReSIL:

* A junction event is two consecutive alignments of one read (in read order) on the
  same contig and strand where the read jumps BACK along the reference: on the plus
  strand it reaches X (end of segment a) and continues at Y < X (start of segment b).
  The candidate circle is [Y, X). Minus-strand reads are mirrored.
* Events are clustered by (contig, start, end): an event joins a cluster when both
  endpoints are within --tolerance bp of the cluster's running median (ONT endpoint
  jitter otherwise splits one junction into several clusters).
* Per cluster:
    junction_reads     reads with at least one event
    multipass_reads    reads with >= 2 events
    beyond_unit_reads  reads that reach past one copy of [Y, X): flanking sequence
                       outside it, or segments covering the same bases twice. A
                       circle is linearized before sequencing, so a read from a single
                       circle never does this; a tandem duplication (Y..X Y..X) does
                       whenever a read is longer than the unit.
    expected_beyond_if_tandem
                       junction_reads x (fraction of mapped reads longer than the
                       unit + 2 x flank): how many beyond-unit reads a tandem copy
                       would have produced in this library.
* label (evidence categories, not proof of topology):
    replicon_junction  end->start of a whole contig (plasmid or closed chromosome);
                       --replicon-slack covers the terminal overlap that assemblers
                       leave on circular contigs (up to the largest k-mer, ~127 bp)
    tandem_or_mixed    >= --min-beyond beyond-unit reads AND >= --min-beyond-frac of
                       expected_beyond_if_tandem (a lone chimeric read does not count)
    circle_supported   >= --min-reads junction reads, beyond-unit reads below the
                       tandem threshold, although
                       a tandem copy would have given >= --min-expected-beyond
    junction_ambiguous enough junction reads, but reads too short to tell a circle
                       from a tandem copy (the same ambiguity as short reads)
    low_support        fewer than --min-reads junction reads

    ont_junctions.py ALN.paf REF.fa.fai OUT.tsv [--support-ids IDS.txt]
"""
import argparse
import bisect
import collections
import csv
import sys


def read_paf(path, min_mapq, min_aln):
    by_read = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ValueError(f"Malformed PAF line: {line[:80]!r}")
            tags = f[12:]
            if any(t.startswith("tp:A:") and t not in ("tp:A:P",) for t in tags):
                continue  # keep primary chains only (supplementary alignments are also tp:A:P)
            qs, qe, ts, te, mapq = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[11])
            if mapq < min_mapq or te - ts < min_aln:
                continue
            by_read[f[0]].append((qs, qe, f[5], ts, te, f[4], int(f[1])))
    return by_read


def events_for_read(segs, max_qgap, min_circle):
    """[(contig, start, end, beyond_bp)] for each back-jump in read order.

    beyond_bp: how far the read reaches past ONE copy of [start, end): flanking sequence
    outside it, or segments a and b covering the same reference bases. A read from a
    single (linearized) circle can never do either; a tandem duplication can.
    """
    segs = sorted(segs)
    out = []
    for a, b in zip(segs, segs[1:]):
        if a[2] != b[2] or a[5] != b[5] or abs(b[0] - a[1]) > max_qgap:
            continue
        if a[5] == "+":
            start, end = b[3], a[4]           # read: ... -> X=a.te, then Y=b.ts
            outside = max(start - a[3], b[4] - end, 0)
        else:
            start, end = a[3], b[4]           # mirrored: reference positions run downward
            outside = max(start - b[3], a[4] - end, 0)
        overlap = max(0, min(a[4], b[4]) - max(a[3], b[3]))
        if end - start >= min_circle:  # a genuine jump back (indel splits give end <= start)
            out.append((a[2], start, end, max(outside, overlap)))
    return out


def cluster(events, tolerance):
    """Cluster on (contig, start, end) against each cluster's running median endpoints.

    Sorted by (contig, start, end); an event joins the nearest recent cluster whose median
    start and end are both within tolerance. Returns clusters of event indices."""
    order = sorted(range(len(events)), key=lambda i: (events[i][1], events[i][2], events[i][3]))
    clusters = []
    for i in order:
        contig, s, e = events[i][1], events[i][2], events[i][3]
        best = None
        for c in reversed(clusters[-50:]):
            if c["contig"] == contig and abs(c["s"] - s) <= tolerance and abs(c["e"] - e) <= tolerance:
                d = abs(c["s"] - s) + abs(c["e"] - e)
                if best is None or d < best[0]:
                    best = (d, c)
        if best:
            c = best[1]
            c["members"].append(i); bisect.insort(c["ss"], s); bisect.insort(c["es"], e)
            c["s"], c["e"] = c["ss"][len(c["ss"]) // 2], c["es"][len(c["es"]) // 2]
        else:
            clusters.append({"contig": contig, "s": s, "e": e, "ss": [s], "es": [e], "members": [i]})
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf"); ap.add_argument("fai"); ap.add_argument("out")
    ap.add_argument("--support-ids")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-aln", type=int, default=200)
    ap.add_argument("--max-qgap", type=int, default=200)
    ap.add_argument("--min-circle", type=int, default=100)
    ap.add_argument("--tolerance", type=int, default=50)
    ap.add_argument("--flank", type=int, default=100)
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("--replicon-slack", type=int, default=150)
    ap.add_argument("--min-expected-beyond", type=float, default=3.0)
    ap.add_argument("--min-beyond", type=int, default=2, help="beyond-unit reads needed to call tandem_or_mixed")
    ap.add_argument("--min-beyond-frac", type=float, default=0.1,
                    help="...and at least this fraction of expected_beyond_if_tandem")
    a = ap.parse_args()
    lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) if l.strip()}
    events = []  # (read, contig, start, end, beyond)
    paf = read_paf(a.paf, a.min_mapq, a.min_aln)
    read_lengths = sorted(segs[0][6] for segs in paf.values())
    for read, segs in paf.items():
        for contig, s, e, beyond in events_for_read(segs, a.max_qgap, a.min_circle):
            L = lengths.get(contig)
            if L is not None and s <= a.replicon_slack and e >= L - a.replicon_slack:
                s, e = 0, L  # end->start of the whole contig: one cluster per replicon
            events.append((read, contig, s, e, beyond))
    rows = []
    support = set()
    for c in cluster(events, a.tolerance):
        members = [events[i] for i in c["members"]]
        per_read = collections.Counter(m[0] for m in members)
        beyond = {m[0] for m in members if m[4] > a.flank}
        starts = sorted(m[2] for m in members); ends = sorted(m[3] for m in members)
        start, end = starts[len(starts) // 2], ends[len(ends) // 2]
        L = lengths.get(c["contig"])
        if L is None:
            raise ValueError(f"Contig {c['contig']} not in {a.fai}")
        # Under a tandem duplication, junction reads longer than one unit (+ flank) would
        # reach beyond it. Expected number: junction reads x library fraction that long.
        need = end - start + 2 * a.flank
        frac_long = (len(read_lengths) - bisect.bisect_left(read_lengths, need)) / max(1, len(read_lengths))
        expected_beyond = len(per_read) * frac_long
        if start <= a.replicon_slack and end >= L - a.replicon_slack:
            label = "replicon_junction"
        elif len(beyond) >= a.min_beyond and len(beyond) >= a.min_beyond_frac * expected_beyond:
            label = "tandem_or_mixed"
        elif len(per_read) < a.min_reads:
            label = "low_support"
        elif expected_beyond >= a.min_expected_beyond:
            label = "circle_supported"
        else:
            label = "junction_ambiguous"  # reads too short to tell a circle from a tandem copy
        support.update(per_read)
        rows.append({"contig": c["contig"], "start": start, "end": end, "length": end - start,
                     "contig_length": L, "junction_reads": len(per_read), "junction_events": len(members),
                     "multipass_reads": sum(1 for n in per_read.values() if n >= 2),
                     "beyond_unit_reads": len(beyond), "expected_beyond_if_tandem": round(expected_beyond, 2),
                     "label": label})
    rows.sort(key=lambda r: (-r["junction_reads"], r["contig"], r["start"]))
    cols = ["contig", "start", "end", "length", "contig_length", "junction_reads", "junction_events",
            "multipass_reads", "beyond_unit_reads", "expected_beyond_if_tandem", "label"]
    with open(a.out, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    if a.support_ids:
        with open(a.support_ids, "w") as handle:
            handle.writelines(r + "\n" for r in sorted(support))
    counts = collections.Counter(r["label"] for r in rows)
    print(f"ont_junctions: {len(events)} events, {len(rows)} clusters, " +
          ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
