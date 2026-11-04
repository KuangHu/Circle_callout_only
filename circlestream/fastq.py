import gzip
import sqlite3
from pathlib import Path


def records(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as handle:
        while True:
            h = handle.readline()
            if not h:
                return
            s, plus, q = (handle.readline().rstrip("\r\n") for _ in range(3))
            h = h.rstrip("\r\n")
            if not h.startswith("@") or not plus.startswith("+") or not s or len(s) != len(q):
                raise ValueError(f"Malformed/truncated four-line FASTQ at {h[:80]!r} in {path}")
            if any(ord(x) < 33 or ord(x) > 126 for x in q):
                raise ValueError("Invalid Phred+33 quality character")
            if not set(s.upper()) <= set("ACGTNRYKMSWBDHV"):
                raise ValueError("Invalid DNA base")
            yield h[1:].split()[0], h, s, q


def pair_id(name):
    return name[:-2] if name.endswith(("/1", "/2")) else name


def materialize(inputs, output, platform):
    """Full validation, preserving original FASTQ; ONT also gets a FASTA for CReSIL."""
    from contextlib import ExitStack
    from itertools import zip_longest
    output.mkdir(exist_ok=True)
    stats = {"reads": 0, "bases": 0, "min_length": None, "max_length": 0}
    with ExitStack() as stack:
        id_path = output / "read_ids.sqlite"
        id_path.unlink(missing_ok=True)
        ids = sqlite3.connect(id_path)
        stack.callback(ids.close)
        ids.execute("CREATE TABLE ids (name TEXT PRIMARY KEY)")
        writers = [stack.enter_context((output / f"reads{i+1}.fastq").open("w", newline="\n")) for i in range(len(inputs))]
        fasta = stack.enter_context((output / "reads.fa").open("w", newline="\n")) if platform == "nanopore" else None
        streams = [records(p) for p in inputs]
        for group in zip_longest(*streams):
            if any(r is None for r in group):
                raise ValueError("R1/R2 record counts differ")
            if len(group) == 2 and pair_id(group[0][0]) != pair_id(group[1][0]):
                raise ValueError("R1/R2 IDs/order differ")
            name = pair_id(group[0][0]) if platform == "illumina" else group[0][0]
            try:
                ids.execute("INSERT INTO ids VALUES (?)", (name,))
            except sqlite3.IntegrityError:
                raise ValueError(f"Duplicate read/template ID: {name}; independent support cannot be counted reliably")
            for writer, (name, header, seq, qual) in zip(writers, group):
                if platform == "illumina" and len(seq) < 75:
                    raise ValueError("Circle_finder ordinary WGS profile requires reads >=75 bp; use a separate preprocessing profile")
                # Normalize pair suffixes only in the derived input, keeping original raw inputs until archive.
                outname = pair_id(name) if platform == "illumina" else name
                writer.write(f"@{outname}\n{seq}\n+\n{qual}\n")
                if fasta:
                    fasta.write(f">{name}\n{seq}\n")
                stats["reads"] += 1
                stats["bases"] += len(seq)
                stats["min_length"] = min(stats["min_length"] or len(seq), len(seq))
                stats["max_length"] = max(stats["max_length"], len(seq))
        ids.commit()
    if not stats["reads"]:
        raise ValueError("Empty FASTQ")
    return stats


def extract(inputs, wanted, output, platform):
    """Keep complete original records and all mates for a union of evidence read IDs."""
    found = set()
    counts = []
    for index, path in enumerate(inputs, 1):
        count = 0
        with gzip.open(output / f"support_reads{index}.fastq.gz", "wt", encoding="ascii") as handle:
            for name, h, s, q in records(path):
                key = pair_id(name) if platform == "illumina" else name
                if key in wanted:
                    handle.write(f"{h}\n{s}\n+\n{q}\n")
                    found.add(key)
                    count += 1
        counts.append(count)
    if found != wanted:
        raise ValueError(f"Could not recover {len(wanted-found)} supporting read IDs; cannot clean")
    if len(counts) == 2 and counts[0] != counts[1]:
        raise ValueError("Incomplete evidence mates")
    return counts
