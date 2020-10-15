#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=96G -l mem_req=96G
#S -l d_rt=256:00:00 -l s_rt=256:00:00
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

"""
Evaluation alignment

Usage:
    eval_align_result.py calc [options] <aligned_tsv> 
    eval_align_result.py agg [options] <result_sqlite>... 

Options:
    --template-dir <PATAH>  : Template file (required run on qsub)
    <aligned_tsv>           : Aligned work TSV file
    <result_sqlite>...      : Caluculated result SQLite file(s)

"""

import os
import sys
from glob import glob
import subprocess

from jinja2 import Environment, FileSystemLoader
from docopt import docopt


def calc(opt):
    template_dir = opt['--template-dir'] if opt['--template-dir'] else os.path.dirname(__file__)

    path_tsv = opt['<aligned_tsv>']

    path_work_dir = os.path.dirname(os.path.abspath(path_tsv))

    paths_fastq = sorted([path for path in glob(f"{path_work_dir}/*") for p in ['.fastq', '.mate1', '.mate2'] if p in os.path.splitext(path)[1]])

    path_tsv_unaligned = f"{paths_fastq[0]}.unaligned.tsv"

    if os.path.exists(path_tsv_unaligned):
        os.remove(path_tsv_unaligned)

    for p in paths_fastq:
        cmd = f"cat {p} | grep '@' | cut -f1 | tr -d '@' >> {path_tsv_unaligned}"
        _ = subprocess.run(cmd, shell=True, capture_output=False)

    path_output = f"{path_tsv}.sqlite"
    path_sql = f"{path_tsv}.sqlite.sql"

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('eval_align_calc.sqlite.sql.j2')

    rendered = template.render(
        tsv_aligned=path_tsv,
        tsv_unaligned=path_tsv_unaligned
        )

    with open(path_sql, "w") as f:
        f.write(rendered)

    cmd = f"sqlite3 {path_output} < {path_sql}"
    sys.stderr.write(cmd)
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    sys.stdout.write(proc.stdout.decode())
    sys.stderr.write(proc.stderr.decode())

    for f in [path_tsv_unaligned, path_sql]:
        os.remove(f)


def agg(opt):
    template_dir = opt['--template-dir'] if opt['--template-dir'] else os.path.dirname(__file__)

    paths_sqlite = opt['<result_sqlite>']
    path_output = os.path.join(
        os.path.dirname(paths_sqlite[0]),
        "eval_align.merged.sqlite"
    )

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('eval_align_merge.sqlite.sql.j2')

    if os.path.exists(path_output):
        os.remove(path_output)

    for p in paths_sqlite:
        rendered = template.render(
            src_database=p
            )

        path_sql = f"{p}.sql"

        with open(path_sql, "w") as f:
            f.write(rendered)

        cmd = f"sqlite3 {path_output} < {path_sql}"
        sys.stderr.write(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)

        sys.stdout.write(proc.stdout.decode())
        sys.stderr.write(proc.stderr.decode())

        for f in [path_sql]:
            os.remove(f)

    path_sql = os.path.join(template_dir, 'eval_align_agg.sqlite.sql')
    cmd = f"sqlite3 {path_output} < {path_sql}"
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    sys.stdout.write(proc.stdout.decode())
    sys.stderr.write(proc.stderr.decode())


def main():
    opt = docopt(__doc__)

    if opt['calc']:
        calc(opt)

    if opt['agg']:
        agg(opt)


if __name__ == '__main__':
    main()
