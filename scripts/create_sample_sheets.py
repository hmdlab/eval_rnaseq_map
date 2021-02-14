#! /usr/bin/env python
#
# Create sample sheets
#
# Usage:
#  this.py <input_dir>
#

import sys
import os
from random import sample, seed


def output(ctrl, case, input_dir, output_prefix):
    group = (["CTRL"] * len(ctrl)) + (["CASE"] * len(case))
    sample = [f"{g}_{str(i).zfill(2)}" for i, g in zip((ctrl + case), group)]
    fastq1 = [
        os.path.abspath(
            os.path.join(input_dir, "sample_{}_1.fastq.gz".format(str(i).zfill(2)))
        )
        for i in (ctrl + case)
    ]
    fastq2 = [
        os.path.abspath(
            os.path.join(input_dir, "sample_{}_2.fastq.gz".format(str(i).zfill(2)))
        )
        for i in (ctrl + case)
    ]

    header = ["sample", "fastq1", "fastq2", "group"]

    with open(f"sample_sheet_{output_prefix}.tsv", "w") as f:
        f.write("\t".join(header))
        f.write("\n")
        for r in zip(sample, fastq1, fastq2, group):
            f.write("\t".join(r))
            f.write("\n")


if __name__ == "__main__":
    input_dir = sys.argv[1]

    try:
        _seed = sys.argv[2]
        seed(_seed)
    except IndexError:
        pass

    ctrl = list(range(1, 25))
    case = list(range(25, 49))

    # #1
    # NOTE: Choose main samples
    nrep = 3
    ctrl_ = sorted(sample(ctrl, nrep))
    case_ = sorted(sample(case, nrep))

    output(ctrl_, case_, os.path.join(input_dir, f"conditions/D20M"), "MAIN")

    # 2 - 5
    lengths = [35, 50, 75]

    for l in lengths:
        output(ctrl_, case_, os.path.join(input_dir, f"conditions/L{l}"), f"L{l}")

    # 6 - 9
    nreads = ["5M", "10M", "40M"]

    for nr in nreads:
        output(ctrl_, case_, os.path.join(input_dir, f"conditions/D{nr}"), f"D{nr}")

    # 10 - 14
    nreps = [2, 4, 8, 12]

    for nrep in nreps:
        ctrl_ = sorted(sample(ctrl, nrep))
        case_ = sorted(sample(case, nrep))

        output(ctrl_, case_, os.path.join(input_dir, f"conditions/D20M"), f"R{nrep}")
