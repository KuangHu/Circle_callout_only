#!/usr/bin/env bash
# Read-only capability checks; does not submit jobs or install packages.
set -euo pipefail
hostname -f
command -v python curl sbatch squeue sacct || true
python --version
python -m snakemake --version
conda run -n circle-short samtools --version
conda run -n cresil cresil identify_wgls --help
sacctmgr show association -p user="$USER"
printf '\nConfirm DTN long-running process policy, filesystem quotas and shared-directory visibility with site documentation/admins.\n'
