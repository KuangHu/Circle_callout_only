#!/usr/bin/env python3
"""Score make_ont_controls.py samples: truth junctions vs junction counter vs CReSIL.

    eval_ont_controls.py SIM_DIR RESULTS_DIR > report.json

RESULTS_DIR/<sample>/ holds ont_junctions.tsv, ont_junction_events.tsv, cresil_trace.tsv,
eccDNA_final.txt. Per truth junction: evidence_reads (simulated reads crossing it with
>= 300 bp each side), counter_reads (those reads among the counter's events at that
junction), label / location / circle_compatible, CReSIL outcome and excluded_at, and
outcome = detected | no_evidence (< 3 evidence reads) | missed (evidence, no counter call).
Also lists every counter junction with >= 3 molecules not explained by truth.
"""
import csv
import json
import sys
from pathlib import Path

EXPECT = {  # sample -> junction -> expected (label set, location)
    "linear": {"pS": ({"junction_supported"}, "whole_contig"), "pM": ({"junction_supported"}, "whole_contig"),
               "pL": ({"junction_supported"}, "whole_contig")},
    "circle": {"c2k": ({"junction_supported"}, "internal"), "c8k": ({"junction_supported"}, "internal"),
               "cIS": ({"ambiguous"}, "internal")},
    "tandem": {"tan3k": ({"tandem_supported"}, "internal")},
    "mixed": {"tan3k": ({"tandem_supported"}, "internal")},
}


def main():
    sim, res = Path(sys.argv[1]), Path(sys.argv[2])
    truth = json.load(open(sim / "truth.json"))
    units = dict(truth["units"], **{k: [k] + v for k, v in truth["replicons"].items()})
    report = {}
    for sample, juncs in truth["samples"].items():
        d = res / sample
        rows = list(csv.DictReader(open(d / "ont_junctions.tsv"), delimiter="\t"))
        events = list(csv.DictReader(open(d / "ont_junction_events.tsv"), delimiter="\t"))
        trace = {r["junction_id"]: r for r in csv.DictReader(open(d / "cresil_trace.tsv"), delimiter="\t")}
        out, explained = {}, set()
        for j in sorted(set(juncs) | set(EXPECT.get(sample, {}))):
            contig, s, e = units[j]
            t = juncs.get(j, {"evidence_reads": 0, "reads_reaching_outside": 0, "read_ids": []})
            hits = [r for r in rows if r["contig"] == contig and abs(int(r["start"]) - s) <= 150 and abs(int(r["end"]) - e) <= 150]
            hits.sort(key=lambda r: -int(r["independent_molecules"]))
            best = hits[0] if hits else None
            explained.update(r["junction_id"] for r in hits)
            ids = set(t["read_ids"])
            got = {ev["read"] for ev in events if best and ev["junction_id"] == best["junction_id"]}
            exp_labels, exp_loc = EXPECT.get(sample, {}).get(j, (set(), None))
            if best and int(best["independent_molecules"]) >= 3:
                outcome = "detected"
            elif t["evidence_reads"] < 3:
                outcome = "no_evidence"
            else:
                outcome = "missed"
            tr = trace.get(best["junction_id"]) if best else None
            out[j] = {"evidence_reads": t["evidence_reads"], "evidence_reaching_outside": t["reads_reaching_outside"],
                      "counter_reads_from_truth": len(ids & got), "counter_other_reads": len(got - ids),
                      "independent_molecules": int(best["independent_molecules"]) if best else 0,
                      "label": best["label"] if best else None, "location": best["location"] if best else None,
                      "circle_compatible": best["circle_compatible"] if best else None,
                      "tandem_flank_molecules": int(best["tandem_flank_molecules"]) if best else 0,
                      "duplicate_junction_rows": max(0, len(hits) - 1),
                      "outcome": outcome,
                      "label_as_expected": (best["label"] in exp_labels) if (best and exp_labels) else None,
                      "location_as_expected": (best["location"] == exp_loc) if (best and exp_loc) else None,
                      "cresil": tr["cresil_ids"] if tr else None, "cresil_excluded_at": tr["excluded_at"] if tr else None}
        extra = [{k: r[k] for k in ("contig", "start", "end", "location", "label", "independent_molecules", "circle_compatible")}
                 for r in rows if r["junction_id"] not in explained and int(r["independent_molecules"]) >= 3]
        report[sample] = {"truth_junctions": out, "unexplained_junctions_ge3": extra}
    json.dump(report, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
