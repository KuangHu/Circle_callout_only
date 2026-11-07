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


MIN_ILLUMINA_LENGTH = 75


def materialize(inputs, output, platform, short_read_policy="reject"):
    """Full validation, preserving original FASTQ; ONT also gets a FASTA for CReSIL.

    Illumina pairs with any mate shorter than MIN_ILLUMINA_LENGTH are rejected
    (default) or, with short_read_policy="drop_pair", dropped as whole pairs.
    Structural errors (truncation, quality length, pairing, duplicate IDs) always fail.
    """
    from contextlib import ExitStack
    from itertools import zip_longest
    if short_read_policy not in ("reject", "drop_pair"):
        raise ValueError(f"Unknown short_read_policy: {short_read_policy}")
    output.mkdir(exist_ok=True)
    stats = {"reads": 0, "bases": 0, "min_length": None, "max_length": 0}
    input_records = dropped = 0
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
            input_records += 1
            if platform == "illumina" and any(len(seq) < MIN_ILLUMINA_LENGTH for _, _, seq, _ in group):
                if short_read_policy == "reject":
                    raise ValueError("Circle_finder ordinary WGS profile requires reads >=75 bp; "
                                     "set caller.short_read_policy=drop_pair in a new run to drop short pairs")
                dropped += 1
                continue
            for writer, (name, header, seq, qual) in zip(writers, group):
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
        raise ValueError("Empty FASTQ" if not input_records else "No reads remain after short-read filtering")
    if platform == "illumina":
        stats.update(short_read_policy=short_read_policy, min_read_length=MIN_ILLUMINA_LENGTH,
                     input_pairs=input_records, dropped_pairs=dropped, retained_pairs=input_records - dropped,
                     retained_pair_fraction=(input_records - dropped) / input_records)
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
