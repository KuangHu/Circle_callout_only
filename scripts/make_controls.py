#!/usr/bin/env python3
"""Small deterministic functional controls; not a realistic accuracy benchmark."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import random


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True)
    a=parser.parse_args()
    out=Path(a.out).resolve()
    out.mkdir(parents=True,exist_ok=False)
    rng=random.Random(1907)
    genome="".join(rng.choice("ACGT") for _ in range(30000))
    start,end=10000,12000
    circle=genome[start:end]
    ref=out/"reference.fa"
    ref.write_text(">chr1\n"+genome+"\n")
    rc=lambda s:s.translate(str.maketrans("ACGT","TGCA"))[::-1]
    manifest=[]
    for platform in ("nanopore","illumina"):
        for case in ("circle","linear","tandem"):
            sid=platform+"_"+case
            folder=out/sid;folder.mkdir()
            reads=[[],[]] if platform=="illumina" else [[]]
            def emit(fragment,i):
                if platform=="illumina":
                    for mate,seq in enumerate((fragment[:150],rc(fragment[-150:]))):
                        reads[mate].append(f"@r{i}/{mate+1}\n{seq}\n+\n{'I'*len(seq)}\n")
                else:
                    # Error-free ONT-like lengths exercise plumbing, not real ONT accuracy.
                    reads[0].append(f"@r{i}\n{fragment}\n+\n{'I'*len(fragment)}\n")
            length=350 if platform=="illumina" else 1000
            background=genome if case!="tandem" else genome[:end]+circle+genome[end:]
            nbackground=1800 if platform=="illumina" else 600
            for i in range(nbackground):
                pos=rng.randrange(len(background)-length+1)
                emit(background[pos:pos+length],i)
            if case=="circle":
                for i in range(300 if platform=="illumina" else 200):
                    pos=rng.randrange(len(circle))
                    emit((circle+circle)[pos:pos+length],nbackground+i)
            files=[]
            for i,records in enumerate(reads,1):
                p=folder/f"reads{i}.fastq.gz"
                # mtime=0 for byte-identical controls across runs.
                with p.open("wb") as raw:
                    with gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as handle:
                        handle.write("".join(records).encode())
                files.append({"name":p.name,"url":"https://example.invalid/controls/"+sid+"/"+p.name,
                    "size_bytes":p.stat().st_size,"md5":hashlib.md5(p.read_bytes()).hexdigest()})
            manifest.append({"sample_id":sid,"accession":"synthetic-functional-control","platform":platform,
                "library":"wgs","reference":str(ref),"peak_gb":4,"files":files})
    (out/"samples.json").write_text(json.dumps(manifest,indent=2))
    (out/"truth.json").write_text(json.dumps({"seed":1907,"reference_region_0based_halfopen":[start,end],
        "circle":"circular fragments added to linear background",
        "linear":"background only; a no-call is expected, not guaranteed by a caller",
        "tandem":"linear tandem insertion; junction can be indistinguishable from a circle",
        "limitations":"error-free reads, single reference, not a sensitivity/specificity or throughput benchmark"},indent=2))
    print(out/"samples.json")


if __name__=="__main__": main()
