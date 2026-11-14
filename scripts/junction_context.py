#!/usr/bin/env python3
"""Genomic context of ONT junction calls: [upstream][circle][downstream] from the sample's own assembly.

    junction_context.py --archives DIR [DIR ...] --out PREFIX [--threads N] [--flank 1000]
    junction_context.py --junctions ont_junctions.tsv --reference REF.fa --sample SID --out PREFIX

Every ont_junctions.py call is an interval [start, end) (0-based, half-open) of the BioSample's
own assembly: reads were aligned to it, so the circle (or tandem unit) sequence IS ref[start:end].
A circle whose sequence is absent from the assembly is never called. This script only reads the
assembly; it changes no call, label or count, and needs no reads (archives keep the reference
path + SHA-256; the reference store is permanent).

Outputs
  PREFIX.tsv.gz   one row per junction (all candidates; filter on support_ge3 / label / location)
  PREFIX.fa.gz    one record per junction with a circle sequence: lowercase upstream flank,
                  UPPERCASE circle ref[start:end], lowercase downstream flank. Further records
                  (|copyN) give the same layout for every other near-identical copy of the circle
                  sequence in the assembly, oriented like the called locus.
  PREFIX.ge3.fa.gz  the same records for junctions with support_ge3 = yes only.
Columns (besides the ont_junctions identity/support columns)
  upstream / downstream  --flank bp left of start / right of end. Truncated at contig ends
                         (upstream_bp / downstream_bp give the actual length); contigs are never
                         wrapped, so contig_edge flanks may be short.
  boundary_repeat        ref[start:start+h] == ref[end:end+h]: bases present at BOTH boundaries
                         (direct repeat / microhomology). The breakpoint can be anywhere inside
                         it; e.g. the core/target bases of an IS110-family excision.
  circle_junction        last --junction-bp of the circle + first --junction-bp: the sequence
                         across the circle's joint, as junction reads show it.
  empty site             upstream + downstream (not stored: concatenate the two columns) is what
                         the genome would read after the unit is excised (one boundary_repeat left).
  genome_copies          near-identical copies of the circle sequence in the assembly (>= 90%
                         query coverage, >= 90% identity, minimap2), the called locus included.
                         > 1: a multi-copy element; the reads place the junction at the called
                         locus, but the source copy (and so which flanks) is not determined by
                         junction reads alone. Empty when not computed (whole_contig, circle
                         shorter than --min-copy-len or longer than --max-circle).
                         Overlapping hits are one locus (the hit at the called locus wins), so a
                         unit that also aligns shifted inside its own tandem array counts once.
  copy_loci              contig:start-end:strand of every copy (capped at --max-loci). Copy
                         boundaries are the alignment extended by the circle's unaligned ends
                         (strand-aware, clipped at contig ends), so every copy spans the full unit.
  called_locus_in_copies sanity check: "yes" when the search found the called locus itself.
whole_contig rows (replicon closures, a technical control) get no flanks or sequence.
"""
import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
import multiprocessing
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ID_COLS = ["junction_id", "contig", "start", "end", "length", "contig_length", "location", "label",
           "circle_compatible", "support_ge3", "independent_molecules", "junction_reads", "microhomology_bp"]
COLS = (["sample"] + ID_COLS + ["upstream_bp", "downstream_bp", "boundary_repeat", "upstream", "downstream",
        "circle_junction", "genome_copies", "copy_loci", "called_locus_in_copies", "in_fasta"])
COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")  # full IUPAC


def revcomp(s):
    return s.translate(COMP)[::-1]


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


def boundary_repeat(seq, start, end, limit=200):
    """Bases repeated at both boundaries: ref[start:start+h] == ref[end:end+h] (junctions are
    stored at their left-most equivalent position, so the repeat extends to the right)."""
    h = 0
    while h < limit and end + h < len(seq) and start + h < end and seq[start + h] == seq[end + h] != "N":
        h += 1
    return seq[start:start + h]


def layout(seq, start, end, flank):
    """(upstream, circle, downstream) with flanks truncated at the contig ends."""
    return seq[max(0, start - flank):start], seq[start:end], seq[end:end + flank]


def genome_copies(queries, reference, minimap2, min_cov=0.9, min_ident=0.9):
    """{query: [(contig, start, end, strand)]} near-identical copies of each query in reference."""
    if not queries:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "q.fa"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()))
        out = subprocess.run([minimap2, "-c", "-x", "asm20", "-s", "50", "--secondary=yes", "-N", "500", "-p", "0.1", "-t", "1",
                              str(reference), str(fa)], check=True, capture_output=True, text=True).stdout
    hits = {k: [] for k in queries}
    for line in out.splitlines():
        f = line.split("\t")
        qlen, qs, qe, matches, alen = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
        if (qe - qs) >= min_cov * qlen and alen and matches / alen >= min_ident:
            # The alignment may stop short of the circle's ends (>= min_cov is enough): project the
            # unaligned query ends onto the target (strand-aware) so the copy spans the whole unit.
            tlen, ts, te, strand = int(f[6]), int(f[7]), int(f[8]), f[4]
            left, right = (qs, qlen - qe) if strand == "+" else (qlen - qe, qs)
            hits[f[0]].append((f[5], max(0, ts - left), min(tlen, te + right), strand))
    return hits


def same_locus(c, s, e, contig, start, end):
    tol = max(50, int(0.05 * (end - start)))
    return c == contig and abs(s - start) <= tol and abs(e - end) <= tol


def collapse(hits, contig, start, end):
    """One locus per set of overlapping hits (a tandem array where the unit also aligns shifted
    by part of a repeat is one locus). The hit at the called locus wins its set, so a shifted
    placement inside the same array never replaces it."""
    own = lambda h: same_locus(h[0], h[1], h[2], contig, start, end)  # noqa: E731
    kept = []
    for h in sorted(hits, key=lambda h: (not own(h), abs(h[1] - start) + abs(h[2] - end) if own(h) else 0, h[0], h[1], -h[2])):
        if not any(k[0] == h[0] and h[1] < k[2] and k[1] < h[2] for k in kept):
            kept.append(h)
    return sorted(kept, key=lambda h: (h[0], h[1]))


def wrap(s, width=80):
    return "\n".join(s[i:i + width] for i in range(0, len(s), width)) + "\n"


def sample_context(sample, rows, reference, o):
    """(tsv_rows, [(support_ge3, fasta_record)]) for one sample's ont_junctions.tsv rows."""
    seqs = read_fasta(reference)
    queries, out, fasta = {}, [], []
    for r in rows:
        start, end = int(r["start"]), int(r["end"])
        if r["location"] != "whole_contig" and o.min_copy_len <= end - start <= o.max_circle:
            queries[r["junction_id"]] = seqs[r["contig"]][start:end]
    copies = genome_copies(queries, reference, o.minimap2) if o.minimap2 else {}
    for r in rows:
        row = {"sample": sample, **{k: r.get(k, "") for k in ID_COLS}}
        seq, start, end = seqs[r["contig"]], int(r["start"]), int(r["end"])
        if r["location"] == "whole_contig":
            out.append(dict(row, in_fasta="no (whole_contig)")); continue
        up, circle, down = layout(seq, start, end, o.flank)
        j = min(o.junction_bp, len(circle) // 2)
        row.update(upstream_bp=len(up), downstream_bp=len(down), boundary_repeat=boundary_repeat(seq, start, end),
                   upstream=up, downstream=down, circle_junction=circle[len(circle) - j:] + circle[:j])
        loci = collapse(copies[r["junction_id"]], r["contig"], start, end) if r["junction_id"] in copies else None
        if loci is not None:
            row.update(genome_copies=len(loci),
                       copy_loci=";".join(f"{c}:{s}-{e}:{st}" for c, s, e, st in loci[:o.max_loci]),
                       called_locus_in_copies="yes" if any(same_locus(c, s, e, r["contig"], start, end)
                                                           for c, s, e, _ in loci) else "no")
        if end - start <= o.max_circle:
            head = (f">{sample}|J{r['junction_id']} {r['contig']}:{start}-{end} length={end - start} "
                    f"location={r['location']} label={r['label']} independent_molecules={r['independent_molecules']} "
                    f"upstream={len(up)} downstream={len(down)} boundary_repeat={row['boundary_repeat'] or '-'} "
                    f"genome_copies={row.get('genome_copies', '') or 'na'}\n")
            ge3 = r.get("support_ge3") == "yes"
            fasta.append((ge3, head + wrap(up.lower() + circle + down.lower())))
            n = 1
            for c, s, e, st in (loci or []):
                if same_locus(c, s, e, r["contig"], start, end):
                    continue
                n += 1
                cu, cc, cd = layout(seqs[c], s, e, o.flank)
                text = cu.lower() + cc + cd.lower()
                fasta.append((ge3, f">{sample}|J{r['junction_id']}|copy{n} {c}:{s}-{e} strand={st} copy_of={r['contig']}:{start}-{end}\n"
                             + wrap(text if st == "+" else revcomp(text))))
            row["in_fasta"] = "yes"
        else:
            row["in_fasta"] = f"no (length > {o.max_circle})"
        out.append(row)
    return out, fasta


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def from_archive(args):
    path, o = args
    sample = Path(path).name[:-len(".tar.gz")]
    want = {"evidence/reference.json": None, "evidence/cresil/ont_junctions.tsv": None, "metadata/spec.json": None}
    with tarfile.open(path, "r:gz") as tar:
        for m in tar:
            if m.name in want:
                want[m.name] = tar.extractfile(m).read().decode()
                if all(v is not None for v in want.values()):
                    break
    platform = json.loads(want["metadata/spec.json"])["sample"]["platform"] if want["metadata/spec.json"] else None
    if platform == "illumina":
        return sample, [], [], None  # not an ONT archive: no junction scan to put in context
    if want["evidence/cresil/ont_junctions.tsv"] is None:
        return sample, [], [], "no ont_junctions.tsv"
    ref = json.loads(want["evidence/reference.json"])
    reference = ref.get("store_path") or ref["path"]
    if not os.path.isfile(reference) or sha256(reference) != ref["sha256"]:
        return sample, [], [], f"reference missing or SHA-256 mismatch: {reference}"
    rows = list(csv.DictReader(io.StringIO(want["evidence/cresil/ont_junctions.tsv"]), delimiter="\t"))
    out, fa = sample_context(sample, rows, reference, o)
    return sample, out, fa, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archives", nargs="+", help="directories of <SAMPLE>.tar.gz archives (Illumina ones are ignored)")
    ap.add_argument("--junctions"); ap.add_argument("--reference"); ap.add_argument("--sample")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--junction-bp", type=int, default=100)
    ap.add_argument("--max-circle", type=int, default=200000, help="longest circle written to FASTA / copy search")
    ap.add_argument("--min-copy-len", type=int, default=100)
    ap.add_argument("--max-loci", type=int, default=50)
    ap.add_argument("--minimap2", default="minimap2", help="minimap2 binary; '' skips the copy search")
    ap.add_argument("--threads", type=int, default=1)
    o = ap.parse_args()
    if bool(o.archives) == bool(o.junctions):
        ap.error("give --archives or --junctions/--reference/--sample")
    samples = problems = n = illumina = 0
    with gzip.open(o.out + ".tsv.gz", "wt", newline="", compresslevel=6) as t, \
            gzip.open(o.out + ".fa.gz", "wt", compresslevel=6) as f, gzip.open(o.out + ".ge3.fa.gz", "wt", compresslevel=6) as f3, \
            (multiprocessing.Pool(o.threads) if o.archives else contextlib.nullcontext()) as pool:
        w = csv.DictWriter(t, fieldnames=COLS, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        if o.junctions:
            rows = list(csv.DictReader(open(o.junctions), delimiter="\t"))
            results = iter([(o.sample, *sample_context(o.sample, rows, o.reference, o), "")])
        else:
            # every completed archive (any sample ID); the platform comes from the archive's spec
            paths = sorted(str(p) for d in o.archives for p in Path(d).glob("*.tar.gz") if not p.name.endswith(".partial.tar.gz"))
            results = pool.imap(from_archive, [(p, o) for p in paths], chunksize=2)
        for sample, rows, fa, problem in results:  # streamed: archives in sorted order
            if problem is None:
                illumina += 1; continue
            samples += 1; n += len(rows)
            if problem:
                problems += 1
                print(f"{sample}: {problem}", file=sys.stderr)
            w.writerows(rows)
            for ge3, rec in fa:
                f.write(rec)
                if ge3:
                    f3.write(rec)
    print(f"{samples} ONT samples ({problems} skipped; {illumina} Illumina archives ignored), {n} junctions -> "
          f"{o.out}.tsv.gz, {o.out}.fa.gz, {o.out}.ge3.fa.gz")


if __name__ == "__main__":
    main()
