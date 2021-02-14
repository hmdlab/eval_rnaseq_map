#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

"""
Create k-mer from transcriptome

Usage:
  create_kmers_from_transcriptome.py [options] <fasta>

Options:
  --k <INT>            : k-mer length [default: 100]
  --output-dir <PATH>  : Output directory [default: .]
  <fasta>              : Transcriptome FASTA file

"""

import sys
import os

from docopt import docopt

from Bio import SeqIO


def main():
    opt = docopt(__doc__)
    path_input = opt["<fasta>"]
    output_dir = opt["--output-dir"]
    length_kmer = int(opt["--k"])

    root, _ = os.path.splitext(os.path.basename(path_input))
    path_output = os.path.join(output_dir, f"{root}.fastq")

    sys.stdout = open(path_output, "w")

    with open(path_input) as f:
        for r in SeqIO.parse(f, "fasta"):
            id_seq = r.id
            seq = r.seq
            length = len(seq)

            if length > length_kmer:
                length_effective = length - length_kmer
                seqs_kmer = [
                    seq[i : i + length_kmer] for i in range(length_effective + 1)
                ]
                ids_kmer = [
                    f"{id_seq}|{i + 1}-{i + length_kmer}"
                    for i in range(length_effective + 1)
                ]
            else:
                seqs_kmer = [seq[0 : length + 1]]
                ids_kmer = [f"{id_seq}|{1}-{length}"]

            for i, s in zip(ids_kmer, seqs_kmer):
                q = "".join(["I"] * len(s))
                print(f"@{i}\n{s}\n+\n{q}")


if __name__ == "__main__":
    main()
