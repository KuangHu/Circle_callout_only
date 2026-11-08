"""SRA object helpers: tool lookup, reference-sequence dependencies, offline conversion.

Runs submitted as aligned reads (e.g. BAM uploads) are stored by SRA as alignments
against external reference sequences ("refseqs"). fasterq-dump needs those sequences
to rebuild the reads and, if they are not next to the .sra file, fetches them over
the network. On LRC compute nodes that fetch can hang indefinitely. The 1-CPU fetch
job (which has working network access) therefore downloads every refseq a run needs
into the same directory as the .sra file, where sra-tools resolve them locally, and
the compute job converts with networking disabled and a time limit.
"""
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

REFSEQ_URL = "https://sra-download.ncbi.nlm.nih.gov/traces/refseq/{}"
SEQID = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.]*")


def sra_tool(name):
    """sra-tools live in the control environment next to the worker's Python."""
    path = Path(sys.executable).parent / name
    if not path.is_file():
        raise FileNotFoundError(f"{name} not found beside {sys.executable}; install sra-tools in the control env")
    return str(path)


def refseq_deps(sra_path):
    """[(seqid, location)] for every external reference sequence of an SRA object.

    Unaligned runs have none (align-info prints only a 'not a database' warning).
    location is 'remote::https://...' when sra-tools would fetch it over the network.
    """
    result = subprocess.run([sra_tool("align-info"), str(sra_path)], capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        raise RuntimeError(f"align-info failed ({result.returncode}): {result.stderr.strip()[:300]}")
    deps = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        seqid, location = line.split(",", 1)[0], line.rsplit(",", 1)[-1]
        if not SEQID.fullmatch(seqid):
            raise ValueError(f"Unexpected refseq ID from align-info: {seqid!r}")
        deps.append((seqid, location))
    return deps


def unresolved(deps):
    return [seqid for seqid, location in deps if location.startswith("remote::")]


def validate_vdb(path):
    result = subprocess.run([sra_tool("vdb-validate"), str(path)], capture_output=True, text=True, timeout=1800)
    text = result.stdout + result.stderr
    return result.returncode == 0 and "is consistent" in text, text


def _curl(url, dest, retries):
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        part.unlink(missing_ok=True)
        result = subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location", "--proto", "=https",
                                 "--proto-redir", "=https", "--connect-timeout", "30", "--speed-limit", "1024",
                                 "--speed-time", "120", "--max-time", "3600", "--output", str(part), url])
        if result.returncode == 0 and part.is_file() and part.stat().st_size > 0:
            ok, _ = validate_vdb(part)
            if ok:
                os.replace(part, dest)
                return
        if attempt + 1 < retries:
            time.sleep(min(30, 2 ** attempt))
    part.unlink(missing_ok=True)
    raise RuntimeError(f"Refseq download/vdb-validate failed after {retries} attempts: {url}")


def fetch_refseqs(sra_path, retries):
    """Download missing refseqs next to sra_path; return one record per dependency.

    Each refseq is itself a VDB table; with no published checksum it is accepted only
    if vdb-validate reports it consistent. Afterwards every dependency must resolve
    locally, or this raises (the sample is then held as DOWNLOAD_FAILED).
    """
    sra_path = Path(sra_path)
    records = []
    for seqid, location in refseq_deps(sra_path):
        dest = sra_path.parent / seqid
        url = REFSEQ_URL.format(seqid)
        if location.startswith("remote::"):
            if location != "remote::" + url:
                raise ValueError(f"Refseq {seqid} resolves to an unexpected location: {location}")
            _curl(url, dest, retries)
            fetched = True
        else:
            fetched = False
        records.append({"seqid": seqid, "url": url if fetched else None, "location_before": location,
                        "size_bytes": dest.stat().st_size if dest.is_file() else None})
    left = unresolved(refseq_deps(sra_path))
    if left:
        raise RuntimeError(f"Refseqs still unresolved after download: {left}")
    return records


def network_isolation():
    """argv prefix that runs a command without network access, if the node allows it."""
    try:
        ok = subprocess.run(["unshare", "-rn", "true"], capture_output=True, timeout=30).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        ok = False
    return ["unshare", "-rn"] if ok else []


def run_limited(argv, timeout_s):
    """Run argv in its own process group; kill the whole group on timeout."""
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        out, err = proc.communicate(timeout=timeout_s)
        return proc.returncode, out, err, False
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate()
        return proc.returncode, out, err, True
