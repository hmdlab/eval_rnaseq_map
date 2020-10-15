#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#
# Extract records from FASTQ/FASTA
#
#  Usage:
#    this.sh <list_target_ids.txt> <fastq/fasta>
#


import sys
import os

import gzip

from Bio import SeqIO


def main():
    path_list_target = sys.argv[1]
    path_fastq = sys.argv[2]

    f = open(path_list_target, 'r')
    target_ids = set([(l.strip()) for l in f])
    f.close()

    def match(id, db=target_ids):
        id_ = id.split('|')[0]
        if id_ in db:
            return True

        return False

    def puts(r, type='fastq'):
        if type in ['fa', 'fasta']:
            print(f">{r.description}\n{r.seq}")
            return

        q = ''.join(['I'] * len(r.seq))
        print(f"@{r.description}\n{r.seq}\n+\n{q}")
        return

    def extract(path):
        root, ext_ = os.path.splitext(path)
        if ext_ == '.gz':
            open_ = gzip.open
            ext = root.split('.')[-1]
        else:
            open_ = open
            ext = ext_[1:]
            if ext == 'fa':
                ext = 'fasta'

        handle = open_(path, 'rt')

        for r in SeqIO.parse(handle, ext):
            if match(r.id):
                puts(r, ext)

        handle.close()

    extract(path_fastq)


if __name__ == '__main__':
    main()
