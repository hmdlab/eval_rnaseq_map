#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=32G -l mem_req=32G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

"""
Calculate genomic coverage

Usage:
  calc_gcov <gtf>

Options:
  <gtf>  : GTF formatted gene annotation file

"""

import sys
import itertools

from docopt import docopt
import numpy as np
from gtfparse import read_gtf


def main():
    options = docopt(__doc__)

    gtf_path = options["<gtf>"]
    sys.stderr.write(gtf_path)

    gtf = read_gtf(gtf_path)
    chromosomes = list(gtf["seqname"].unique())

    for c in chromosomes:
        _gtf = gtf.query(f"seqname == '{c}' & feature == 'exon'")
        _gtf = _gtf[["start", "end"]].drop_duplicates()
        _ranges = []

        for s, e in zip(_gtf["start"], _gtf["end"]):
            _ranges = np.unique(np.append(_ranges, range(s, e + 1)))

        sys.stdout.write(f"{c}\t{len(_ranges)}\n")


if __name__ == "__main__":
    main()
