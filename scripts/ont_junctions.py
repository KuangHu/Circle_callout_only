#!/usr/bin/env python3
"""Genome-wide circularization-junction scan from long-read alignments (minimap2 -c PAF).

Independent of CReSIL and of any insertion or short-read call: every read is examined.

    ont_junctions.py ALN.paf REF.fa OUT.tsv [--events EVENTS.tsv] [--support-ids IDS.txt]

ALN.paf must come from `minimap2 -c -x map-ont --secondary=yes -N 5 -p 0.9`; secondary
alignments are used only to flag repeat / multi-mapping ambiguity. REF.fa: the BioSample's
own assembly (the same file the other callers use).

Junction event
    Two consecutive primary alignments of one read (in read order), same contig and strand,
    query gap <= --max-qgap, where the read jumps BACK on the reference: it reaches X (end of
    segment a) and continues at Y < X (start of segment b), X - Y >= --min-circle. Candidate
    unit [Y, X). Minus-strand reads are mirrored. Every event is written to --events with raw
    coordinates, anchors, mapping qualities, repeat flags and outside-flank lengths.
    No read has to cover the whole unit and no second pass is required.

Per event
    anchor_a / anchor_b   aligned reference bp of the segments before / after the junction
    unique_a / unique_b   mapq >= --min-mapq and no secondary alignment of the same read
                          bases scoring >= --repeat-score-frac of the primary
    insertion evidence    minimap2 often aligns a long read through a tandem duplication as ONE
                          alignment carrying the extra copy as a CIGAR insertion. A primary,
                          unique alignment with an insertion of length X-Y (+/- 5% or
                          --tolerance) placed within [Y, X] (+/- --tolerance), that continues
                          >= --flank bp outside [Y, X], is a linear tandem molecule: counted in
                          insertion_tandem_molecules and in tandem_flank_molecules. An
                          insertion fitting several nearby junctions is credited to one only
                          (most split-read events, then closest length): no double counting.
    outside_bp            bp of the segment reaching X (or starting at Y) aligned OUTSIDE
                          [Y, X): the same alignment continues from the unit into flanking
                          genome. unique_outside_bp counts it only on a unique segment; that
                          links the junction to unique outside sequence, i.e. a linear form
                          (tandem copy, or an integrated copy) exists.
    overlap_bp            bp covered by both segments (the read passes the same bases twice).
                          Recorded only: multimeric circles, rolling-circle products and
                          tandem arrays all produce it, so it proves no topology.

Clusters and equivalent breakpoints
    Events cluster when both endpoints are within --tolerance bp of a cluster's running
    median. Each cluster's junction is shifted to its left-most equivalent position using the
    reference: with microhomology (ref[Y-i] == ref[X-i], i = 1..h) junctions (Y-i, X-i) are
    indistinguishable. Clusters that coincide afterwards are merged; their raw medians and IDs
    stay in raw_clusters, and every event keeps its raw coordinates.

Independent support
    junction_reads          distinct read IDs (split-read events + matching insertion reads)
    independent_molecules   reads grouped when the reference positions where they start and
                            stop (the molecule's cut ends) agree within --dup-tolerance bp,
                            strand ignored: duplex template/complement pairs and duplicate
                            records count once. Thresholds use this number. All passes of one
                            read are one molecule, so it never exceeds junction_reads. On
                            replicons/circles shorter than the reads, unaligned partial end
                            passes put many chains' ends on the unit boundaries, so distinct
                            reads can merge: a conservative lower bound there.

Labels (evidence categories, never topology proof)
    tandem_supported    >= --min-tandem-reads independent molecules link the junction to
                        unique outside sequence (unique_outside_bp >= --flank). Does not
                        exclude circles of the same unit in the same culture.
    junction_supported  >= 1 independent molecule with both anchors unique and >= --min-anchor,
                        and not tandem_supported.
    ambiguous           junction reads only with non-unique (repeat / multi-mapping) or short
                        anchors.
    circle_compatible   column, not a label: "yes" when junction_supported, no molecule reaches
                        outside sequence, and the library's read lengths would have given
                        >= --min-expected-flank such molecules under a simple tandem; "no" when
                        any molecule reaches outside; "na" otherwise. Compatible != confirmed.
location: whole_contig (end -> start of a whole contig within --edge bp, which also absorbs the
terminal overlap assemblers leave on circular contigs), contig_edge (one end within --edge),
internal. support_ge3: independent_molecules >= 3 (primary-result threshold).
"""
import argparse
import bisect
import collections
import csv
import re
import sys

COLS = ["junction_id", "contig", "start", "end", "length", "contig_length", "location", "raw_clusters",
        "microhomology_bp", "junction_reads", "junction_events", "independent_molecules", "split_molecules",
        "anchored_molecules",
        "repeat_ambiguous_molecules", "tandem_flank_molecules", "insertion_tandem_molecules", "multipass_molecules",
        "expected_flank_if_tandem",
        "label", "circle_compatible", "support_ge3"]
EVENT_COLS = ["read", "evidence", "junction_id", "cluster_raw", "contig", "strand", "raw_start", "raw_end", "start", "end",
              "anchor_a", "anchor_b", "mapq_a", "mapq_b", "unique_a", "unique_b", "outside_bp",
              "unique_outside_bp", "overlap_bp", "molecule_start", "molecule_end"]


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def tag(fields, key, default=None):
    for t in fields:
        if t.startswith(key):
            return t[len(key):]
    return default


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_paf(path, min_aln, insertion_min=0):
    """{read: [primaries, secondaries, {query lengths}]}; primaries include supplementary (tp:A:P)."""
    reads = collections.defaultdict(lambda: ([], [], set()))
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ValueError(f"Malformed PAF line: {line[:80]!r}")
            tp = tag(f[12:], "tp:A:", "P")
            ts, te = int(f[7]), int(f[8])
            seg = {"qs": int(f[2]), "qe": int(f[3]), "contig": f[5], "ts": ts, "te": te, "strand": f[4],
                   "mapq": int(f[11]), "score": int(tag(f[12:], "AS:i:", tag(f[12:], "ms:i:", f[9]))),
                   "ins": []}
            if tp == "P" and insertion_min:
                pos = ts
                for n, op in CIGAR.findall(tag(f[12:], "cg:Z:", "")):
                    n = int(n)
                    if op in "MDN=X":
                        pos += n
                    elif op == "I" and n >= insertion_min:
                        seg["ins"].append((pos, n))
            entry = reads[f[0]]
            entry[2].add(int(f[1]))
            if tp == "P" and te - ts >= min_aln:
                entry[0].append(seg)
            elif tp == "S":
                entry[1].append(seg)
    return reads


def is_unique(seg, secondaries, min_mapq, frac):
    if seg["mapq"] < min_mapq:
        return False
    for s in secondaries:
        ov = min(seg["qe"], s["qe"]) - max(seg["qs"], s["qs"])
        if ov >= 0.5 * (seg["qe"] - seg["qs"]) and s["score"] >= frac * seg["score"]:
            return False
    return True


def events_for_read(prim, sec, o):
    prim = sorted(prim, key=lambda s: (s["qs"], s["qe"]))
    # Chains: consecutive segments on one contig and strand without a large read gap (the
    # condition for an event below). The molecule's cut ends are where the chain starts and stops
    # on the reference. Every pass of a multipass read shares them; using each pass's own segment
    # ends (or their min/max, always the whole unit) would make distinct multipass reads of a
    # small circle look like duplicates of one another.
    chain = [0]
    for a, b in zip(prim, prim[1:]):
        same = a["contig"] == b["contig"] and a["strand"] == b["strand"] and abs(b["qs"] - a["qe"]) <= o.max_qgap
        chain.append(chain[-1] + (not same))
    first, last = {}, {}
    for s, c in zip(prim, chain):
        first.setdefault(c, s); last[c] = s
    out = []
    for i, (a, b) in enumerate(zip(prim, prim[1:])):
        if chain[i] != chain[i + 1]:
            continue
        f, l = first[chain[i]], last[chain[i]]
        if a["strand"] == "+":
            start, end = b["ts"], a["te"]       # read: ... -> X = a.te, then continues at Y = b.ts
            at_x, at_y = a, b                   # segment ending at X / segment starting at Y
            mol = (f["ts"], l["te"])            # where the read starts and stops on the reference
        else:
            start, end = a["ts"], b["te"]       # mirrored: a starts (ref) at Y, b ends at X
            at_x, at_y = b, a
            mol = (f["te"], l["ts"])
        if end - start < o.min_circle:
            continue                            # indel split or forward jump, not a back-jump
        ua = is_unique(a, sec, o.min_mapq, o.repeat_score_frac)
        ub = is_unique(b, sec, o.min_mapq, o.repeat_score_frac)
        uniq = {id(a): ua, id(b): ub}
        out_x = max(0, start - at_x["ts"])      # segment reaching X began left of Y
        out_y = max(0, at_y["te"] - end)        # segment starting at Y runs past X
        out.append({"contig": a["contig"], "strand": a["strand"], "start": start, "end": end,
                    "anchor_a": a["te"] - a["ts"], "anchor_b": b["te"] - b["ts"],
                    "mapq_a": a["mapq"], "mapq_b": b["mapq"], "unique_a": ua, "unique_b": ub,
                    "outside_bp": max(out_x, out_y),
                    "unique_outside_bp": max(out_x if uniq[id(at_x)] else 0, out_y if uniq[id(at_y)] else 0),
                    "overlap_bp": max(0, min(a["te"], b["te"]) - max(a["ts"], b["ts"])),
                    "molecule_start": min(mol), "molecule_end": max(mol)})
    return out


def cluster(events, tolerance):
    order = sorted(range(len(events)), key=lambda i: (events[i]["contig"], events[i]["start"], events[i]["end"]))
    clusters = []
    for i in order:
        e = events[i]
        best = None
        for c in reversed(clusters[-50:]):
            if c["contig"] == e["contig"] and abs(c["s"] - e["start"]) <= tolerance and abs(c["e"] - e["end"]) <= tolerance:
                d = abs(c["s"] - e["start"]) + abs(c["e"] - e["end"])
                if best is None or d < best[0]:
                    best = (d, c)
        if best:
            c = best[1]
            c["members"].append(i); bisect.insort(c["ss"], e["start"]); bisect.insort(c["es"], e["end"])
            c["s"], c["e"] = c["ss"][len(c["ss"]) // 2], c["es"][len(c["es"]) // 2]
        else:
            clusters.append({"contig": e["contig"], "s": e["start"], "e": e["end"],
                             "ss": [e["start"]], "es": [e["end"]], "members": [i]})
    return clusters


def microhomology(seq, start, end, limit=200):
    """(left, right): junction (start, end) equals (start-i, end-i) for i <= left and
    (start+i, end+i) for i <= right, because the reference bases at both ends are identical."""
    left = 0
    while left < limit and start - left - 1 >= 0 and seq[start - left - 1] == seq[end - left - 1] != "N":
        left += 1
    right = 0
    while right < limit and end + right < len(seq) and seq[start + right] == seq[end + right] != "N":
        right += 1
    return left, right


def molecules(members, tol):
    """Group events into physical molecules. All events of one read are ONE molecule first (a
    read crossing the junction twice yields two events but is a single observation); reads are
    then merged when their reference start/stop (events_for_read: the chain's cut ends, shared
    by every pass) agree within tol (duplex strands, duplicate records)."""
    by_read = collections.defaultdict(list)
    for e in members:
        by_read[e["read"]].append(e)
    reads = sorted((((min(e["molecule_start"] for e in evs), max(e["molecule_end"] for e in evs)), name, evs)
                    for name, evs in by_read.items()), key=lambda t: (t[0], t[1]))
    reads = [(key, evs) for key, _, evs in reads]
    groups = []
    for key, evs in reads:
        for g in groups:
            if abs(g[0][0] - key[0]) <= tol and abs(g[0][1] - key[1]) <= tol:
                g[1].extend(evs); break
        else:
            groups.append((key, list(evs)))
    return [g[1] for g in groups]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf"); ap.add_argument("ref"); ap.add_argument("out")
    ap.add_argument("--events"); ap.add_argument("--support-ids")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-aln", type=int, default=200)
    ap.add_argument("--min-anchor", type=int, default=300)
    ap.add_argument("--max-qgap", type=int, default=200)
    ap.add_argument("--min-circle", type=int, default=100)
    ap.add_argument("--tolerance", type=int, default=50)
    ap.add_argument("--flank", type=int, default=100)
    ap.add_argument("--edge", type=int, default=150)
    ap.add_argument("--dup-tolerance", type=int, default=10)
    ap.add_argument("--repeat-score-frac", type=float, default=0.9)
    ap.add_argument("--min-tandem-reads", type=int, default=2)
    ap.add_argument("--min-expected-flank", type=float, default=3.0)
    o = ap.parse_args()

    seqs = read_fasta(o.ref)
    lengths = {k: len(v) for k, v in seqs.items()}
    reads = read_paf(o.paf, o.min_aln, insertion_min=o.min_circle)
    read_lengths = sorted(max(q) for _, _, q in reads.values() if q)
    dup_names = sum(1 for _, _, q in reads.values() if len(q) > 1)

    events = []
    for name, (prim, sec, _) in reads.items():
        for e in events_for_read(prim, sec, o):
            L = lengths[e["contig"]]
            e["raw_start"], e["raw_end"], e["read"] = e["start"], e["end"], name
            if e["start"] <= o.edge and e["end"] >= L - o.edge:
                e["start"], e["end"] = 0, L     # end -> start of a whole contig: one closure per replicon
            events.append(e)

    insertions = collections.defaultdict(list)  # contig -> sorted [(pos, len, read, ts, te, unique)]
    for name, (prim, sec, _) in reads.items():
        for seg in prim:
            for pos, n in seg["ins"]:
                insertions[seg["contig"]].append((pos, n, name, seg["ts"], seg["te"],
                                                  is_unique(seg, sec, o.min_mapq, o.repeat_score_frac)))
    for v in insertions.values():
        v.sort()
    loci = {}
    for cid, c in enumerate(cluster(events, o.tolerance), 1):
        seq, L = seqs[c["contig"]], lengths[c["contig"]]
        s, e = c["s"], c["e"]
        left, right = (0, 0) if (s == 0 and e == L) else microhomology(seq, s, e)
        key = (c["contig"], s - left, e - left)
        for m in c["members"]:
            events[m]["cluster_raw"] = cid
        l = loci.setdefault(key, {"raw": [], "members": [], "mh": 0})
        l["raw"].append((cid, s, e)); l["members"] += c["members"]; l["mh"] = max(l["mh"], left + right)

    # Each insertion is evidence for ONE junction: when several nearby junctions could explain
    # it (same tandem, slightly different split-read breakpoints), it goes to the one with the
    # most split-read events, then the closest unit length, so molecules are never double counted.
    owner = {}  # (contig, index in insertions[contig]) -> (rank key, jid)
    ordered = sorted(loci.items())
    for jid, (key, l) in enumerate(ordered, 1):
        contig, start, end = key
        L, unit = lengths[contig], end - start
        if start <= o.edge and end >= L - o.edge:
            continue
        col = insertions.get(contig, [])
        lo = bisect.bisect_left(col, (start - o.tolerance,))
        for i in range(lo, bisect.bisect_right(col, (end + o.tolerance + 1,))):
            pos, n, name, ts, te, uniq = col[i]
            if uniq and abs(n - unit) <= max(o.tolerance, 0.05 * unit) and (ts <= start - o.flank or te >= end + o.flank):
                rank = (-len(l["members"]), abs(n - unit), jid)
                if (contig, i) not in owner or rank < owner[(contig, i)][0]:
                    owner[(contig, i)] = (rank, jid)
    owned = collections.defaultdict(list)
    for (contig, i), (_, jid) in owner.items():
        owned[jid].append(insertions[contig][i])

    rows, support, extra_events = [], set(), []
    for jid, (key, l) in enumerate(ordered, 1):
        contig, start, end = key
        L = lengths[contig]
        members = [events[m] for m in l["members"]]
        for m in members:
            m["junction_id"] = jid
        ins_ev = [{"read": name, "contig": contig, "strand": ".", "evidence": "insertion",
                   "raw_start": pos, "raw_end": pos + n, "start": start, "end": end,
                   "anchor_a": pos - ts, "anchor_b": te - pos, "unique_a": True, "unique_b": True,
                   "outside_bp": max(start - ts, te - end, 0), "unique_outside_bp": max(start - ts, te - end, 0),
                   "overlap_bp": 0, "molecule_start": ts, "molecule_end": te, "junction_id": jid, "cluster_raw": ""}
                  for pos, n, name, ts, te, _ in sorted(owned.get(jid, []))]
        ins_mols = molecules(ins_ev, o.dup_tolerance)
        extra_events.extend(ins_ev)
        mols = molecules(members, o.dup_tolerance)
        all_mols = molecules(members + ins_ev, o.dup_tolerance)
        anchored = [g for g in mols if any(m["unique_a"] and m["unique_b"] and min(m["anchor_a"], m["anchor_b"]) >= o.min_anchor for m in g)]
        split_tandem = [g for g in mols if any(m["unique_outside_bp"] >= o.flank for m in g)]
        tandem = molecules([g[0] for g in split_tandem] + [g[0] for g in ins_mols], o.dup_tolerance)
        repeat = [g for g in mols if all(not (m["unique_a"] and m["unique_b"]) for m in g)]
        multipass = sum(1 for g in mols if len({m["read"] for m in g}) < len(g) or any(m["overlap_bp"] > o.flank for m in g))
        need = (end - start) + 2 * o.flank
        expected = len(all_mols) * (len(read_lengths) - bisect.bisect_left(read_lengths, need)) / max(1, len(read_lengths))
        if start <= o.edge and end >= L - o.edge:
            location = "whole_contig"
        elif start <= o.edge or end >= L - o.edge:
            location = "contig_edge"
        else:
            location = "internal"
        if len(tandem) >= o.min_tandem_reads:
            label = "tandem_supported"
        elif anchored:
            label = "junction_supported"
        else:
            label = "ambiguous"
        if tandem:
            compatible = "no"
        elif label == "junction_supported" and expected >= o.min_expected_flank:
            compatible = "yes"
        else:
            compatible = "na"
        names = sorted({m["read"] for m in members} | {m["read"] for m in ins_ev})
        support.update(names)
        rows.append({"junction_id": jid, "contig": contig, "start": start, "end": end, "length": end - start,
                     "contig_length": L, "location": location,
                     "raw_clusters": ";".join(f"{cid}:{s}-{e}" for cid, s, e in l["raw"]),
                     "microhomology_bp": l["mh"], "junction_reads": len(names), "junction_events": len(members),
                     "independent_molecules": len(all_mols), "split_molecules": len(mols),
                     "anchored_molecules": len(anchored),
                     "repeat_ambiguous_molecules": len(repeat), "tandem_flank_molecules": len(tandem),
                     "insertion_tandem_molecules": len(ins_mols),
                     "multipass_molecules": multipass, "expected_flank_if_tandem": round(expected, 2),
                     "label": label, "circle_compatible": compatible,
                     "support_ge3": "yes" if len(all_mols) >= 3 else "no"})
    rows.sort(key=lambda r: (-r["independent_molecules"], r["contig"], r["start"]))
    with open(o.out, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    if o.events:
        with open(o.events, "w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=EVENT_COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            for e in events:
                e.setdefault("evidence", "split")
            w.writeheader(); w.writerows(sorted(events + extra_events, key=lambda e: (e["junction_id"], e["evidence"], e["read"])))
    if o.support_ids:
        with open(o.support_ids, "w") as handle:
            handle.writelines(r + "\n" for r in sorted(support))
    counts = collections.Counter((r["location"], r["label"]) for r in rows)
    print(f"ont_junctions: {len(reads)} reads, {len(events)} events, {len(rows)} junctions, "
          f"read names with differing lengths {dup_names}; " +
          ", ".join(f"{k[0]}/{k[1]}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
