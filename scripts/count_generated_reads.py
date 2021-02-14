#! /usr/bin/env python3
#
# Usage:
#  find *.fastq.gz | xargs python this.py
#
#$ -S $HOME/.pyenv/shims/python
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

import sys

from collections import Counter
from functools import partial

from Bio import SeqIO
import gzip

fastq_paths = sys.argv[1:]

for fastq_path in fastq_paths:
    print("Now processing: {}".format(fastq_path))
    if not fastq_path.endswith(".gz"):
        _open = open
    else:
        _open = partial(gzip.open, mode="rt")

    transcript_ids = []
    with _open(fastq_path) as f:
        for record in SeqIO.parse(f, "fastq"):
            derived_from = record.id.split(";")[0].split("/")[1].split("|")[0]
            transcript_ids.append(derived_from)

    transcript_id_counts = Counter(transcript_ids)

    with open("{}.counted.txt".format(fastq_path), "w") as f:
        for id, count in transcript_id_counts.most_common():
            f.write("{}\t{}\n".format(id, count))
