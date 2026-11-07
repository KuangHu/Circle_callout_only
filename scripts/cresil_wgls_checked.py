#!/usr/bin/env python3
"""Pinned CReSIL WGLS: fail on subprocess errors; explicit no-candidate exits.

Compatible with the upstream Python 3.8 environment. Detection logic is unchanged.
"""
import json
from pathlib import Path
import subprocess
import sys
import traceback

NO_CANDIDATES = {
    "[ABORT] no potential merge region was detected",
    "[ABORT] no merge region was detected (lower -breakpoint, -depth; or region too short)",
    "[ABORT] no eccDNA was detected (lower -breakpoint, -depth; or region too short)"
}
# Pinned upstream filters columns of BedTool.to_dataframe() results without checking
# for emptiness. An empty frame has no columns, so zero breakpoints/regions raise
# KeyError at exactly these lines (depth and linkage modes) instead of exiting.
EMPTY_FRAME_LINES = {
    "df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name']>0]",
    "df_filtered_breaks_bdg = df_filtered_breaks_bdg[df_filtered_breaks_bdg['name'] >= breakpointmerge]",
    "df_count_merged_window = df_count_merged_window[df_count_merged_window['name'] > 1]",
    "df_rem_bg_aln_reads_merge = df_rem_bg_aln_reads_merge[['chrom','start','end']]",
}


def empty_region_error(error):
    if not isinstance(error, KeyError):
        return False
    frames = traceback.extract_tb(error.__traceback__)
    own = [i for i, f in enumerate(frames) if Path(f.filename).parts[-3:] == ("cresil", "cli", "identify_wgls.py")]
    if not own:
        return False
    last = frames[own[-1]]
    # The KeyError must be raised by pandas column lookup, not by upstream code itself.
    inner = frames[own[-1] + 1:]
    return (last.name == "main" and (last.line or "").strip() in EMPTY_FRAME_LINES
            and bool(inner) and all("pandas" in Path(f.filename).parts for f in inner))


def record_no_candidates(args, reason):
    out = Path(args.triminput).parent
    (out / "no_candidates.json").write_text(json.dumps({
        "status": "no_candidates_at_configured_thresholds", "upstream_exit": reason,
        "not_proof_of_absence": True
    }, indent=2))
    (out / "eccDNA_final.txt").write_text(
        "id\tmerge_region\tmerge_len\tnum_region\tctc\tnumreads\ttotalbase\tcoverage\n")
    print("circle-stream: normal no-candidate result; reason recorded", flush=True)


def strict_call(command, *args, **kwargs):
    # Upstream uses shell pipelines and ignores subprocess.call return codes.
    # Bash pipefail plus check=True prevents a failed aligner becoming a zero call.
    if kwargs.pop("shell", False):
        if not isinstance(command, str):
            raise TypeError("Unexpected upstream shell command")
        command = ["bash", "-euo", "pipefail", "-c", command]
    return subprocess.run(command, *args, check=True, **kwargs).returncode


def execute(module, args):
    original = subprocess.call
    subprocess.call = strict_call
    try:
        try:
            module.main(args)
        except SystemExit as error:
            if error.code in (None, 0):
                return
            reason = str(error.code).strip()
            if reason not in NO_CANDIDATES:
                raise
            record_no_candidates(args, reason)
        except KeyError as error:
            if not empty_region_error(error):
                raise
            record_no_candidates(args, "zero potential regions (upstream empty-frame KeyError)")
    finally:
        subprocess.call = original


if __name__ == "__main__":
    from cresil.cli import identify_wgls
    execute(identify_wgls, identify_wgls.argparser().parse_args(sys.argv[1:]))
