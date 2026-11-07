import contextlib
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import time

REPO = Path(__file__).resolve().parents[1]
GB = 10**9


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


# S3 multipart ETag = md5(concat(md5(part_i))) + "-<parts>". NCBI's SRA ODP bucket
# uses 8 MiB parts; other common client defaults are tried when the count fits.
S3_PART_SIZES = tuple(m * 1024 * 1024 for m in (8, 16, 5, 15, 32, 64, 100))
ETAG = re.compile(r"[0-9a-f]{32}(?:-[1-9][0-9]{0,4})?")


def s3_etag_matches(path, etag):
    etag = etag.strip('"').lower()
    if not ETAG.fullmatch(etag):
        return False
    if "-" not in etag:
        return digest(path, "md5") == etag
    parts = int(etag.split("-")[1])
    size = Path(path).stat().st_size
    for part in S3_PART_SIZES:
        if -(-size // part) != parts:
            continue
        sums = []
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(part), b""):
                sums.append(hashlib.md5(block).digest())
        if f"{hashlib.md5(b''.join(sums)).hexdigest()}-{len(sums)}" == etag:
            return True
    return False


def is_sra(item):
    return item["name"].endswith(".sra")


def source_ok(path, item):
    """Size plus the source checksum: ENA MD5 for FASTQ, S3 ETag for SRA objects."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size != item["size_bytes"]:
        return False
    if item.get("md5"):
        return digest(path, "md5") == item["md5"].lower()
    return s3_etag_matches(path, item["s3_etag"])


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def workflow_fingerprint():
    files = list((REPO / "circlestream").glob("*.py")) + list((REPO / "scripts").glob("*.py")) + [REPO / "stream.py", REPO / "run_rule.py", REPO / "workflow" / "Snakefile", REPO / "third_party" / "Circle_finder" / "original.sh"]
    return fingerprint({str(p.relative_to(REPO)): digest(p) for p in sorted(files)})


@contextlib.contextmanager
def lock(path):
    """Crash-visible lock. Never steal automatically across hosts."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir()
    except FileExistsError:
        raise RuntimeError(f"Lock exists: {path}. Check its owner and running jobs before removing it.")
    try:
        write_json(path / "owner.json", dict(host=socket.gethostname(), pid=os.getpid(), time=time.time()))
        yield
    finally:
        (path / "owner.json").unlink(missing_ok=True)
        path.rmdir()


def safe_id(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", value):
        raise ValueError(f"Unsafe sample ID: {value!r}")
    return value


def validate_sample(sample):
    safe_id(sample["sample_id"])
    platform = sample["platform"]
    if platform not in ("nanopore", "illumina") or sample["library"] != "wgs":
        raise ValueError("v0.1 supports ordinary WGS only; enriched/RCA/ATAC need separate validated profiles")
    files = sample["files"]
    sra = any(is_sra(f) for f in files)
    if sra and len(files) != 1:
        raise ValueError("An SRA object is one file per run (paired reads are split during QC)")
    if not sra and len(files) != (2 if platform == "illumina" else 1):
        raise ValueError("Use one ONT FASTQ or one ordered R1/R2 pair per manifest entry")
    if len({f["name"] for f in files}) != len(files):
        raise ValueError("Duplicate file names")
    for item in files:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*(?:\.f(?:ast)?q(?:\.gz)?|\.sra)", item["name"]):
            raise ValueError("Input must be named .fastq[.gz], .fq[.gz] or <run>.sra, without directories")
        if not item["url"].startswith("https://") or any(x in item["url"] for x in ("REPLACE", "\n", "\r")):
            raise ValueError("Supply a public HTTPS URL (ENA FASTQ or NCBI SRA on AWS)")
        if not isinstance(item.get("size_bytes"), int) or item["size_bytes"] <= 0:
            raise ValueError("Source MD5 and positive compressed file size are required")
        if is_sra(item) and not item.get("md5"):
            if not ETAG.fullmatch(str(item.get("s3_etag", "")).strip('"').lower()):
                raise ValueError("SRA objects need the S3 ETag (or an MD5) as source checksum")
        elif not re.fullmatch(r"[0-9a-fA-F]{32}", str(item.get("md5", ""))):
            raise ValueError("Source MD5 and positive compressed file size are required")
    if sample["peak_gb"] * GB <= sum(f["size_bytes"] for f in files):
        raise ValueError("peak_gb must reserve raw + unpacked + mapping + caller + evidence workspace")
    ref = Path(sample["reference"])
    if not ref.is_absolute() or not ref.is_file() or ref.suffix == ".gz":
        raise ValueError(f"Reference must be an existing, uncompressed FASTA: {ref}")


def load_config(path):
    c = read_json(path)
    for key in ("state_dir", "scratch_dir", "archive_dir"):
        p = Path(c[key])
        if not p.is_absolute() or "CHANGE_ME" in str(p):
            raise ValueError(f"Configure an absolute {key}")
        c[key] = str(p.resolve())
    roots = [Path(c[k]) for k in ("state_dir", "scratch_dir", "archive_dir")]
    for i, a in enumerate(roots):
        for b in roots[i + 1:]:
            if a == b or a in b.parents or b in a.parents:
                raise ValueError("state, scratch, archive must be separate, non-nested directories")
    for key in ("scratch_budget_gb", "min_free_gb", "max_active_samples", "max_jobs", "poll_seconds", "download_retries"):
        if c[key] <= 0:
            raise ValueError(f"{key} must be positive")
    if c["slurm"]["cpus"] < 4:
        raise ValueError("CReSIL has internal minimum worker counts: allocate at least 4 CPUs")
    if c["caller"]["cresil_mode"] not in ("depth", "linkage"):
        raise ValueError("cresil_mode must be depth or linkage")
    if c.get("download_mode", "dtn") not in ("dtn", "slurm"):
        raise ValueError("download_mode must be dtn (DTN downloader) or slurm (1-CPU fetch jobs + node-local work)")
    if c.get("reference_store"):
        rs = Path(c["reference_store"])
        if not rs.is_absolute():
            raise ValueError("reference_store must be an absolute path")
        if any(rs == Path(c[k]) or Path(c[k]) in rs.parents for k in ("scratch_dir",)):
            raise ValueError("reference_store must not live inside scratch (it is permanent)")
    if local_mode(c):
        f = c.get("fetch")
        if not isinstance(f, dict) or not all(k in f for k in ("max_jobs", "batch_gb", "time", "memory_gb")):
            raise ValueError("slurm mode needs fetch: {max_jobs, batch_gb, time, memory_gb}")
        if f["max_jobs"] <= 0 or f["batch_gb"] <= 0 or f["memory_gb"] <= 0 or f.get("max_samples", 10) <= 0:
            raise ValueError("fetch limits must be positive")
    if c["caller"].get("short_read_policy", "reject") not in ("reject", "drop_pair"):
        raise ValueError("caller.short_read_policy must be reject or drop_pair")
    return c


def state_paths(c, sid):
    sid = safe_id(sid)
    return Path(c["state_dir"]) / sid, Path(c["scratch_dir"]) / sid


def local_mode(c):
    return c.get("download_mode", "dtn") == "slurm"


def work_root(c, sid):
    """Where a sample's working files live.

    dtn mode: the shared scratch directory (raw + work together).
    slurm mode: node-local disk of the running job; shared scratch holds only the
    downloaded source objects. Lustre metadata latency (seconds per operation on LRC)
    makes per-step small-file work on shared scratch the dominant cost.
    """
    if not local_mode(c):
        return state_paths(c, sid)[1]
    jid = os.environ.get("SLURM_JOB_ID", "")
    if not jid.isdigit():
        raise RuntimeError("Node-local work requires a Slurm allocation")
    return Path(os.environ.get("CS_LOCAL_ROOT", "/tmp")) / f"circle-stream-{jid}" / safe_id(sid)


def set_status(state, status, **extra):
    write_json(Path(state) / "status.json", dict(status=status, updated=time.time(), **extra))


def tree_size(path):
    total = 0
    for root, dirs, files in os.walk(path, followlinks=False):
        dirs[:] = [d for d in dirs if not (Path(root) / d).is_symlink()]
        for name in files:
            p = Path(root) / name
            try:
                if not p.is_symlink():
                    total += p.stat().st_size
            except FileNotFoundError:
                pass
    return total
