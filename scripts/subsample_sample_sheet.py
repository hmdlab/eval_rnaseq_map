#! /usr/bin/env python
#
# Subsample
#
# Usage:
#  this.py <sample_sheet>
#

import sys
import os
from random import sample, seed
from copy import deepcopy

import pandas as pd


def subsample(sample_sheet, population, n_rep, n_trial=1) -> list:
    list_ = []
    for t in range(n_trial):
        sampled = {k: sorted(sample(v, n_rep)) for k, v in population.items()}

        indexes = [v for value in sampled.values() for v in value]

        list_.append(sample_sheet.iloc[indexes])

    return list_


def subsample_mock(sample_sheet, population, n_rep, n_trial=1) -> list:
    list_ = []

    for t in range(n_trial):
        dict_ = {}
        _sample_sheet = deepcopy(sample_sheet)
        for k, v in population.items():
            sampled = sorted(sample(v, 2 * n_rep))
            ctrl = sorted(sample(sampled, n_rep))
            case = set(sampled) - set(ctrl)

            _sample_sheet.loc[ctrl, "group"] = "CTRL"
            _sample_sheet.loc[case, "group"] = "CASE"

            dict_[k] = _sample_sheet.iloc[sampled]

        list_.append(dict_)

    return list_


if __name__ == "__main__":
    sample_sheet_path = sys.argv[1]

    root, ext = os.path.splitext(os.path.basename(sample_sheet_path))

    # NOTE: Conditions
    seed(12345)
    n_rep = 12
    n_trial = 3

    sample_sheet = pd.read_table(sample_sheet_path, dtype=str)
    population = {
        g: sample_sheet.query(f"group == '{g}'").index.tolist()
        for g in set(sample_sheet.group)
    }

    # NOTE: Subsample and output
    for i, s in enumerate(subsample(sample_sheet, population, n_rep, n_trial)):
        s.sort_values(["sample"]).sort_values(["group"], ascending=False).to_csv(
            f"{root}_sub_{n_rep}_{i}", sep="\t", index=False
        )

    # NOTE: Generate mock comparisions
    for i, d in enumerate(subsample_mock(sample_sheet, population, n_rep, n_trial)):
        for k, v in d.items():
            v.sort_values(["sample"]).sort_values(["group"], ascending=False).to_csv(
                f"{root}_mocksub_{n_rep}_{i}_{k}{ext}", sep="\t", index=False
            )
