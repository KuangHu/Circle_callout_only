#!/usr/bin/env python3
"""Deterministic simulated ONT controls with per-read truth for the junction counter and CReSIL.

    make_ont_controls.py --out DIR [--seed 20261108]

One reference (the "assembly") for all samples: a 300 kb chromosome carrying three copies of
a 1.3 kb IS-like element, plus three plasmid contigs - pS 5 kb, pM 60 kb (written with 100 bp
of terminal overlap, as assemblers leave on circular contigs) and pL 120 kb (above CReSIL's
100 kb linkage window). Samples (reads with ~5% random errors, lognormal lengths, both strands):
  linear   chromosome + plasmids only; internal junctions: none expected
  circle   + 2 kb circle (chr 50000-52000), 8 kb circle (chr 200000-208000) and a circle of
           IS copy 1 (junction inside a repeat: expected ambiguous)
  tandem   the sample genome carries a 3 kb tandem duplication (chr 150000-153000 twice) that
           the reference collapses to one copy; no circles
  mixed    tandem as above + circles of the same 3 kb unit
Circle molecules are linearized at a random point before sequencing, so their reads never
exceed one unit; plasmid reads likewise never exceed the plasmid length.

truth.json lists, per sample and junction (reference coordinates), the reads that cross it
with >= 300 bp on both sides ("evidence available"), and how many of those also reach
>= 100 bp of unique sequence outside the unit. Functional controls, not an accuracy benchmark.
"""
import argparse
import gzip
import hashlib
import json
import math
import random
from pathlib import Path

ANCHOR, FLANK = 300, 100


def manifest_entry(sid, fq, reference):
    """url: placeholder required by `init` (never fetched); path: the input for `stage-local`."""
    return {"sample_id": f"ontctl_{sid}", "accession": "synthetic-ont-control", "platform": "nanopore",
            "library": "wgs", "reference": str(reference), "peak_gb": 4,
            "files": [{"name": fq.name, "url": f"https://example.invalid/ont-controls/{sid}/{fq.name}",
                       "path": str(fq), "size_bytes": fq.stat().st_size, "md5": hashlib.md5(fq.read_bytes()).hexdigest()}]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True); ap.add_argument("--seed", type=int, default=20261108)
    a = ap.parse_args()
    out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(a.seed)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    comp = str.maketrans("ACGT", "TGCA")
    rc = lambda s: s.translate(comp)[::-1]  # noqa: E731

    IS = rnd(1300)
    chrom = list(rnd(300000))
    is_pos = [30000, 120000, 250000]
    for p in is_pos:
        chrom[p:p + len(IS)] = IS
    chrom = "".join(chrom)
    plasmids = {"pS": rnd(5000), "pM": rnd(60000), "pL": rnd(120000)}
    copies = {"pS": 4, "pM": 1, "pL": 1}
    ref = {"chr": chrom, "pS": plasmids["pS"], "pM": plasmids["pM"] + plasmids["pM"][:100], "pL": plasmids["pL"]}
    (out / "reference.fa").write_text("".join(f">{k}\n{v}\n" for k, v in ref.items()))

    T0, T1 = 150000, 153000
    units = {"c2k": ("chr", 50000, 52000), "c8k": ("chr", 200000, 208000),
             "cIS": ("chr", is_pos[0], is_pos[0] + len(IS)), "tan3k": ("chr", T0, T1)}

    def readlen():
        return max(500, min(60000, int(math.exp(rng.gauss(math.log(6000), 0.8)))))

    def mutate(s):
        out_, i = [], 0
        for ch in s:
            r = rng.random()
            if r < 0.03:
                out_.append(rng.choice("ACGT".replace(ch, "")))
            elif r < 0.04:
                out_.append(ch + rng.choice("ACGT"))
            elif r < 0.05:
                continue
            else:
                out_.append(ch)
        return "".join(out_)

    manifest, truth = [], {"seed": a.seed, "reference_contigs": {k: len(v) for k, v in ref.items()},
                           "units": {k: list(v) for k, v in units.items()},
                           "replicons": {"pS": [0, 5000], "pM": [0, 60100], "pL": [0, 120000]},
                           "notes": "pM reference has 100 bp terminal overlap; pL exceeds CReSIL -mwindow 100000; "
                                    "cIS junction lies inside a 3-copy repeat", "samples": {}}
    for sid in ("linear", "circle", "tandem", "mixed"):
        tandem = sid in ("tandem", "mixed")
        genome = chrom[:T1] + chrom[T0:T1] + chrom[T1:] if tandem else chrom
        records, spans = [], {}

        def emit(seq, name, crossings):
            seq = mutate(seq)
            if rng.random() < 0.5:
                seq = rc(seq)
            records.append(f"@{name}\n{seq}\n+\n{'5' * len(seq)}\n")
            for j, outside in crossings:
                spans.setdefault(j, []).append((name, outside))

        n = 0
        for _ in range(int(40 * len(genome) / 7000)):  # ~40x chromosome
            L = readlen(); p = rng.randrange(0, len(genome) - L)
            cross = []
            if tandem and p + ANCHOR <= T1 and p + L - ANCHOR >= T1:  # crosses copy1 -> copy2 junction
                cross.append(("tan3k", (T0 - p >= FLANK) or (p + L - (T1 + (T1 - T0)) >= FLANK)))
            emit(genome[p:p + L], f"g{n}", cross); n += 1
        for name, seq in plasmids.items():
            for _ in range(int(40 * copies[name] * len(seq) / 7000)):
                L = min(readlen(), len(seq)); p = rng.randrange(len(seq))
                frag = (seq + seq)[p:p + L]
                cross = [(name, False)] if p + ANCHOR <= len(seq) and p + L - ANCHOR >= len(seq) else []
                emit(frag, f"{name}_{n}", cross); n += 1
        circles = {"circle": {"c2k": 60, "c8k": 40, "cIS": 60}, "mixed": {"tan3k": 60}}.get(sid, {})
        for u, molecules in circles.items():
            contig, s, e = units[u]
            unit = chrom[s:e]
            for _ in range(molecules):
                L = min(readlen(), len(unit)); p = rng.randrange(len(unit))
                cross = [(u, False)] if p + ANCHOR <= len(unit) and p + L - ANCHOR >= len(unit) else []
                emit((unit + unit)[p:p + L], f"{u}_{n}", cross); n += 1
        rng.shuffle(records)
        folder = out / sid; folder.mkdir()
        fq = folder / "reads1.fastq.gz"
        with fq.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as h:
            h.write("".join(records).encode())
        manifest.append(manifest_entry(sid, fq, out / "reference.fa"))
        truth["samples"][sid] = {j: {"evidence_reads": len(v), "reads_reaching_outside": sum(o for _, o in v),
                                     "read_ids": sorted(r for r, _ in v)} for j, v in sorted(spans.items())}
    (out / "samples.json").write_text(json.dumps(manifest, indent=2))
    (out / "truth.json").write_text(json.dumps(truth, indent=2))
    print(json.dumps({s: {j: (v["evidence_reads"], v["reads_reaching_outside"]) for j, v in d.items()}
                      for s, d in truth["samples"].items()}))


if __name__ == "__main__":
    main()
