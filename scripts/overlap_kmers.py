#! /usr/bin/env python3
# $ -S $HOME/.pyenv/shims/python3
# $ -l epyc
# $ -l s_vmem=256G -l mem_req=256G
# $ -cwd
# $ -o ugelogs/
# $ -e ugelogs/

"""
Overlap kmers

Usage:
    overlap_kmers.py uniq <fasta>
    overlap_kmers.py agg <work_tsv>...

Options:
    <fastq>     : kmers in FASA format
    <work_tsv>  : Work TSV

"""

import os
import csv

import gzip
from collections import Counter
from itertools import combinations

import pandas as pd
from Bio import SeqIO
from docopt import docopt


def uniq(opt):
    path_input = opt['<fasta>']

    root, ext = os.path.splitext(path_input)
    path_output = "{}.uniqued.tsv".format(root)

    if ext == '.gz':
        open_ = gzip.open
        fmt = root.split('.')[-1]
    else:
        open_ = open
        fmt = ext[1:]

    if fmt == 'fa':
        fmt = 'fasta'

    f = open_(path_input, 'rt')

    kmers = {}

    for r in SeqIO.parse(f, fmt):
        seq_ = str(r.seq)

        try:
            kmers[seq_] += 1
        except KeyError:
            kmers[seq_] = 1

    f.close()

    with open(path_output, 'w') as f:
        for kv in kmers.items():
            print("\t".join(map(str, kv)), file=f)


def agg(opt):
    paths_input = opt['<work_tsv>']
    output_dir = os.path.dirname(paths_input[0])
    path_output = os.path.join(output_dir, "aggregated.tsv")
    n_sets = len(paths_input)

    def intersect(paths):
        kmers = {}

        for i, p in enumerate(paths):
            base_ = 2 ** i

            reader = pd.read_csv(
                p, delimiter='\t', chunksize=(10**6), header=None)

            for c in reader:
                for i in c[0].tolist():
                    try:
                        kmers[i] += base_
                    except KeyError:
                        kmers[i] = base_

        return kmers

    kmers_merged = intersect(paths_input)
    counts_kmer_merged = Counter(v for k, v in kmers_merged.items())
    bases = [2 ** i for i in range(n_sets)]

    combinations_ = sum([list(combinations(bases, i))
                         for i in range(1, n_sets + 1)], [])

    def subtotal(counts, base):
        try:
            print(base)
            return sum([v for k, v in counts.items() if (k & base) == base])
        except KeyError:
            return 0

    counts_hierarchical = {}
    for c in combinations_:
        counts_hierarchical[sum(c)] = subtotal(counts_kmer_merged, sum(c))

    def key_converted(key, n):
        converted = 'n'
        for i in range(n):
            base_ = 2 ** i

            if key & base_:
                converted += str(i + 1)

        return converted

    with open(path_output, 'w') as f:
        for k, v in counts_hierarchical.items():
            key_ = key_converted(k, n_sets)
            print('\t'.join([key_, str(v)]), file=f)


def main():
    opt = docopt(__doc__)

    if opt['uniq']:
        uniq(opt)

    if opt['agg']:
        agg(opt)


if __name__ == '__main__':
    main()
