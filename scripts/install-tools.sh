#!/usr/bin/env bash
# Run on a site-approved installation/transfer host with internet access.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
conda env create -f "$repo/environment.yml"
conda env create -f "$repo/envs/short.yml"
conda env create -f "$repo/envs/cresil-upstream.yml"
mkdir -p "$repo/vendor"
git clone https://github.com/visanuwan/cresil.git "$repo/vendor/cresil"
git -C "$repo/vendor/cresil" checkout --detach 78356affa6adca6a41ea297c3d7d93f5ecb5ae6a
conda run -n cresil python -m pip install --no-deps "$repo/vendor/cresil"
conda list -n circle-control --explicit > "$repo/vendor/circle-control.explicit.txt"
conda list -n circle-short --explicit > "$repo/vendor/circle-short.explicit.txt"
conda list -n cresil --explicit > "$repo/vendor/cresil.explicit.txt"
printf 'Installed. Activate circle-control; copy and edit the two example JSON files.\n'
