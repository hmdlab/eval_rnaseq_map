#! /usr/bin/env python3
#
# Convert FASTA to FASTQ
#
# Usage:
#   qsub -V this.sh <FASTA>
#
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=8G -l mem_req=8G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

import os
import sys
import re


def main():
    path_fasta = sys.argv[1]
    re_fa_header = re.compile("^>(.+)")
    path_fastq = f"{os.path.splitext(path_fasta)[0]}.fastq"

    fa_header = ""
    fa_seq = ""
    qa_header = "+\n"
    qa_seq = ""

    with open(path_fasta, "r") as fh:
        with open(path_fastq, "w") as fo:
            for buf in fh:
                buf.rstrip("rn")
                m = re_fa_header.match(buf)
                if m:
                    if fa_header != "":
                        print(fa_header, fa_seq, qa_header, qa_seq, sep="", file=fo)
                    fa_header = "@" + m.group(1) + "\n"
                    fa_seq = ""
                    qa_seq = ""
                else:
                    fa_seq = fa_seq + buf
                    qa_seq = qa_seq + "".join(["I"] * (len(buf) - 1))

        print(fa_header, fa_seq, qa_header, qa_seq, sep="", file=fo)

    sys.stderr.write(f"File wrote: {path_fastq}")


if __name__ == '__main__':
    main()
