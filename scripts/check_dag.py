#!/usr/bin/env python3
"""Parse the real Snakefile on Linux without running bioinformatics tools."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as temp:
    root=Path(temp)
    config=json.loads((repo/"config.example.json").read_text())
    spec={"repo":str(repo),"python":sys.executable,"config":config}
    path=root/"spec.json"
    path.write_text(json.dumps(spec))
    subprocess.run([sys.executable,"-m","snakemake","--snakefile",str(repo/"workflow"/"Snakefile"),
        "--directory",str(root),"--cores","8","--dry-run","--printshellcmds","--config",f"spec={path}","--","archive.ok"],check=True)
