#!/usr/bin/env python3
"""ONT-first per-BioSample summary: junction counter vs CReSIL, with Illumina cross-support.

    compare_platforms.py --ont DIR [DIR ...] [--illumina DIR [DIR ...]] --out OUTDIR
                         [--tol 50] [--frac 0.02]

--ont dirs hold <SID>_ont.tar.gz (or <SID>.tar.gz) archives from the v3 pipeline; a later dir
wins for a repeated BioSample (reruns). --illumina dirs hold <SID>_ill.tar.gz or <SID>.tar.gz
Circle_finder archives; Illumina is optional cross-support, never a precondition.

Matching (same contig): endpoints within max(tol, frac x length) of each other, OR shifted
together along the junction's microhomology (|d_start - d_end| <= tol and |d_start| <=
microhomology_bp + tol), since such junctions are indistinguishable.

Outputs
  ont_junctions.tsv  every ONT junction, all labels, with CReSIL trace outcome, Illumina
                     split-read support (0 = Illumina present, no matching call; blank = no
                     Illumina run for this BioSample)
  per_sample.tsv     counts per BioSample
  summary.json       totals by location x label at several molecule thresholds, ONT-only vs
                     Illumina-supported, CReSIL agreement, CReSIL exclusion steps
Labels are evidence categories; nothing here is a confirmed circle.
"""
import argparse
import collections
import csv
import io
import json
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cresil_trace import final_spans, span_matches  # noqa: E402  same matching rule as the tracer

THRESHOLDS = (1, 2, 3, 5, 10)


def archive_map(dirs, suffix):
    out = {}
    for d in dirs:
        for p in sorted(Path(d).glob("*.tar.gz")):
            name = p.name[:-7]
            if name.endswith(".partial"):
                continue
            if name.endswith(suffix):
                out[name[: -len(suffix)]] = p
            elif not name.endswith(("_ill", "_ont")):
                out[name] = p  # archives named by bare BioSample
    return out


def read_members(path, names):
    got = {}
    with tarfile.open(path) as tar:
        for m in tar:
            for key, want in names.items():
                if m.name == want or m.name.endswith("/" + want):
                    got[key] = tar.extractfile(m).read().decode()
    return got


def tsv(text):
    return list(csv.DictReader(io.StringIO(text or ""), delimiter="\t"))


def matches(j, c, tol, frac):
    return span_matches(j["contig"], int(j["start"]), int(j["end"]), int(j.get("microhomology_bp") or 0), c, tol, frac)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ont", nargs="+", required=True); ap.add_argument("--illumina", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=int, default=50); ap.add_argument("--frac", type=float, default=0.02)
    a = ap.parse_args()
    ont = archive_map(a.ont, "_ont")
    ill = archive_map(a.illumina, "_ill")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows, per = [], []
    for sid, path in sorted(ont.items()):
        m = read_members(path, {"j": "evidence/cresil/ont_junctions.tsv", "t": "evidence/cresil/cresil_trace.tsv",
                                "e": "evidence/cresil/eccDNA_final.txt"})
        if "j" not in m:
            per.append({"sample": sid, "status": "not a v3 ONT archive"}); continue
        junc, trace = tsv(m["j"]), {r["junction_id"]: r for r in tsv(m.get("t"))}
        cres = []
        if "e" in m:
            ef = out / ".eccDNA_final.tmp"
            ef.write_text(m["e"]); cres = final_spans(ef, 100); ef.unlink()
        icalls = None
        if sid in ill:
            im = read_members(ill[sid], {"jt": "sample-ref.microDNA-JT.txt"})
            icalls = [(f[0], int(f[1]), int(f[2]), int(f[3])) for f in
                      (l.split("\t") for l in im.get("jt", "").splitlines()) if len(f) >= 4]
        for j in junc:
            hit_c = [c for c in cres if matches(j, c, a.tol, a.frac)]
            hit_i = [c for c in (icalls or []) if matches(j, c, a.tol, a.frac)]
            t = trace.get(j["junction_id"], {})
            excl = t.get("excluded_at", "")
            if excl in ("kept_cyclic", "cyclic") and not hit_c:
                excl = "in_final_graph"  # archives traced before the spanning rule: overlap only
            if excl == "kept_cyclic" and not any(c[4] for c in hit_c):
                excl = "cyclic"
            rows.append(dict(j, sample=sid,
                             cresil_match=";".join(sorted({c[3] + ("(cyclic)" if c[4] else "") for c in hit_c})) or "none",
                             cresil_excluded_at=excl,
                             illumina_split_reads="" if icalls is None else max([c[3] for c in hit_i], default=0)))
        sr = [r for r in rows if r["sample"] == sid]
        per.append({"sample": sid, "status": "ok", "illumina": "yes" if icalls is not None else "no",
                    "junctions": len(sr),
                    "internal_ge3": sum(r["location"] != "whole_contig" and int(r["independent_molecules"]) >= 3 for r in sr),
                    "whole_contig": sum(r["location"] == "whole_contig" for r in sr),
                    "cresil_cyclic": len({c[3] for c in cres if c[4]}),
                    "illumina_calls": "" if icalls is None else len(icalls),
                    "illumina_calls_ge3": "" if icalls is None else sum(c[3] >= 3 for c in icalls)})
    cols = ["sample"] + [k for k in (rows[0] if rows else {}) if k != "sample"]
    with open(out / "ont_junctions.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    pcols = list(dict.fromkeys(k for p in per for k in p))
    with open(out / "per_sample.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=pcols, delimiter="\t", lineterminator="\n", restval=""); w.writeheader(); w.writerows(per)

    def count(pred):
        return dict(sorted(collections.Counter(f"{r['location']}/{r['label']}" for r in rows if pred(r)).items()))
    ge = lambda r, t: int(r["independent_molecules"]) >= t  # noqa: E731
    with_ill = lambda r: r["illumina_split_reads"] != ""  # noqa: E731
    s = {"samples": sum(p["status"] == "ok" for p in per), "samples_with_illumina": sum(p.get("illumina") == "yes" for p in per),
         "by_min_independent_molecules": {f">={t}": count(lambda r, t=t: ge(r, t)) for t in THRESHOLDS},
         "primary_ge3": {"all": count(lambda r: ge(r, 3)),
                         "ont_only_where_illumina_exists": count(lambda r: ge(r, 3) and with_ill(r) and int(r["illumina_split_reads"]) == 0),
                         "illumina_supported_ge1": count(lambda r: ge(r, 3) and with_ill(r) and int(r["illumina_split_reads"]) >= 1),
                         "illumina_supported_ge3": count(lambda r: ge(r, 3) and with_ill(r) and int(r["illumina_split_reads"]) >= 3),
                         "cresil_cyclic_match": count(lambda r: ge(r, 3) and "(cyclic)" in r["cresil_match"]),
                         "cresil_excluded_at": dict(collections.Counter(f"{r['location']}/{r['cresil_excluded_at']}" for r in rows if ge(r, 3))),
                         "circle_compatible_yes": count(lambda r: ge(r, 3) and r["circle_compatible"] == "yes")},
         "match_rule": f"same contig; endpoints within max({a.tol} bp, {a.frac} x length) or shifted together within microhomology + {a.tol} bp",
         "caveat": "labels are evidence categories, not confirmed circles; whole_contig rows are replicon closures (technical control)"}
    json.dump(s, open(out / "summary.json", "w"), indent=2)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
