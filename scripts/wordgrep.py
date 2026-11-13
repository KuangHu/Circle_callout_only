#!/usr/bin/env python3
"""Exact, fast replacement for `LC_ALL=C grep -w -F -f PATTERNS FILE`.

Circle_finder (upstream, unchanged in third_party/) filters its read tables with
`grep -w -Ff ids.txt table.txt`. GNU grep degrades badly when -w is combined with
hundreds of thousands of fixed patterns that are prefixes of each other (read IDs
such as ERR1.1, ERR1.12, ERR1.123): it re-tries shorter matches at every failed
word boundary, and single calls took up to 2 h on 1 GB runs. The adapter in
circlestream/rules.py runs this script in place of those lines.

Semantics reproduced exactly (C locale): a line is printed, in input order, if some
pattern occurs in it with a non-word byte (anything but [A-Za-z0-9_]) or the line
edge on both sides. Patterns cannot contain tabs here (they would make the per-field
search below incomplete), and an empty pattern matches differently in GNU grep; in
either case this script hands the work to grep itself so the result stays identical.

Set CS_WORDGREP_VERIFY=1 to also run GNU grep and fail on any byte difference
(validation only; it brings back grep's run time).
"""
import os
import re
import subprocess
import sys

WORD = bytearray(256)
for ch in b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_":
    WORD[ch] = 1
ALL_WORD = re.compile(rb"[A-Za-z0-9_]*").fullmatch
NON_WORD = re.compile(rb"[^A-Za-z0-9_]")


def grep_argv(patterns_path, table_path):
    return ["grep", "-w", "-F", "-f", patterns_path, table_path]


def load_patterns(path):
    data = open(path, "rb").read()
    if not data:
        return set()
    if data.endswith(b"\n"):
        data = data[:-1]
    return set(data.split(b"\n"))


def matches(line, pats, lens):
    for field in line.split(b"\t"):
        if field in pats:
            return True
        if ALL_WORD(field):
            continue  # only the whole field is word-bounded inside a pure-word field
        cuts = [m.start() for m in NON_WORD.finditer(field)]
        starts = [0] + [k + 1 for k in cuts]
        ends = cuts + [len(field)]
        for s in starts:
            for e in ends:
                if e > s and (e - s) in lens and field[s:e] in pats:
                    return True
    return False


def filter_lines(patterns_path, table_path, out):
    pats = load_patterns(patterns_path)
    if b"" in pats or any(b"\t" in p for p in pats):
        ref = subprocess.run(grep_argv(patterns_path, table_path), capture_output=True, env=dict(os.environ, LC_ALL="C"))
        out.write(ref.stdout)
        return ref.returncode
    lens = {len(p) for p in pats}
    found = 0
    if pats:
        with open(table_path, "rb") as handle:
            for raw in handle:
                line = raw[:-1] if raw.endswith(b"\n") else raw
                if matches(line, pats, lens):
                    out.write(line + b"\n")
                    found += 1
    return 0 if found else 1


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: wordgrep.py PATTERNS FILE > OUT")
    patterns_path, table_path = sys.argv[1:]
    if os.environ.get("CS_WORDGREP_VERIFY") == "1":
        import io
        buf = io.BytesIO()
        rc = filter_lines(patterns_path, table_path, buf)
        ref = subprocess.run(grep_argv(patterns_path, table_path), capture_output=True, env=dict(os.environ, LC_ALL="C"))
        same = ref.stdout == buf.getvalue() and ref.returncode == rc
        ours, theirs = buf.getvalue().count(b"\n"), ref.stdout.count(b"\n")
        verdict = "IDENTICAL" if same else "DIFFERENT"
        print(f"wordgrep-verify {verdict} {table_path} lines={ours} grep_lines={theirs} rc={rc} grep_rc={ref.returncode}",
              file=sys.stderr, flush=True)
        if not same:
            sys.exit(3)
        sys.stdout.buffer.write(buf.getvalue())
        sys.exit(0)
    rc = filter_lines(patterns_path, table_path, sys.stdout.buffer)
    sys.exit(0 if rc in (0, 1) else rc)  # no match is not an error (as in the adapter's grep wrapper)


if __name__ == "__main__":
    main()
