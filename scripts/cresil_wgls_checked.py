#!/usr/bin/env python3
"""Pinned CReSIL WGLS: fail on subprocess errors; explicit no-candidate exits.

Compatible with the upstream Python 3.8 environment. Detection logic is unchanged.
"""
import json
from pathlib import Path
import subprocess
import sys

NO_CANDIDATES = {
    "[ABORT] no potential merge region was detected",
    "[ABORT] no merge region was detected (lower -breakpoint, -depth; or region too short)",
    "[ABORT] no eccDNA was detected (lower -breakpoint, -depth; or region too short)"
}


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
            out = Path(args.triminput).parent
            (out / "no_candidates.json").write_text(json.dumps({
                "status": "no_candidates_at_configured_thresholds", "upstream_exit": reason,
                "not_proof_of_absence": True
            }, indent=2))
            (out / "eccDNA_final.txt").write_text(
                "id\tmerge_region\tmerge_len\tnum_region\tctc\tnumreads\ttotalbase\tcoverage\n")
            print("circle-stream: normal no-candidate result; reason recorded", flush=True)
    finally:
        subprocess.call = original


if __name__ == "__main__":
    from cresil.cli import identify_wgls
    execute(identify_wgls, identify_wgls.argparser().parse_args(sys.argv[1:]))
