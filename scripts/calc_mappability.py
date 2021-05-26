#! /usr/bin/env python3
#$ -S $HOME/.pyenv/shims/python3
#$ -l s_vmem=256G -l mem_req=256G
#$ -l d_rt=480:00:00 -l s_rt=480:00:00
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

"""
Calculate mappability

Usage:
    calc_mappability.py split [options] <bam>
    calc_mappability.py calc [options] <annotation_sqlite> <aligned_tsv>
    calc_mappability.py agg [options] <annotation_sqlite> <unaligned_fastq> <result_sqlite>...

Options:
    --template-dir <PATAH>  : Template file (required run on qsub)
    --output-dir <PATH>     : Output directory [default: .]
    --line-size <INT>       : Split line size [default: 100000000]
    <annotation_sqlite>     : Annotation SQLite file
    <aligned_tsv>           : Aligned work TSV file
    <unaligned_fastq>       : Unaligned FASTQ
    <result_sqlite>...      : Result of calculation SQLite file(s)

"""

import os
import sys
import subprocess

from jinja2 import Environment, FileSystemLoader
from docopt import docopt


def split(opt):
    path_bam = opt['<bam>']
    output_dir = opt['--output-dir']
    line_size = opt['--line-size']

    root, _ = os.path.splitext(os.path.basename(path_bam))
    output_subdir = os.path.join(output_dir, f"{root}")

    if os.path.exists(output_subdir):
        sys.stderr.write('File already exists.')
        return 1

    os.makedirs(output_subdir)
    prefix_output = os.path.join(output_dir, f"{root}", f"{root}.aligned")

    cmd = f"samtools view {path_bam} | cut -f1,3 | split -l {line_size} - {prefix_output}. --additional-suffix=.tsv"
    sys.stderr.write(cmd)
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    sys.stdout.write(proc.stdout.decode())
    sys.stderr.write(proc.stderr.decode())


def calc(opt):
    template_dir = opt['--template-dir'] if opt['--template-dir'] else os.path.dirname(__file__)
    output_dir = opt['--output-dir']

    path_annotation = opt['<annotation_sqlite>']
    path_tsv = opt['<aligned_tsv>']

    root, _ = os.path.splitext(os.path.basename(path_tsv))
    path_output = os.path.join(output_dir, f"{root}.sqlite")
    path_sql = os.path.join(output_dir, f".{root}.sql")

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('calc_mappability_calc.sqlite.sql.j2')

    rendered = template.render(
        tsv_aligned=path_tsv,
        annotation=path_annotation
        )

    with open(path_sql, "w") as f:
        f.write(rendered)

    cmd = f"sqlite3 {path_output} < {path_sql}"
    sys.stderr.write(cmd)
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    sys.stdout.write(proc.stdout.decode())
    sys.stderr.write(proc.stderr.decode())

    for f in [path_sql]:
        os.remove(f)


def agg(opt):
    template_dir = opt['--template-dir'] if opt['--template-dir'] else os.path.dirname(__file__)
    output_dir = opt['--output-dir']

    paths_sqlite = opt['<result_sqlite>']
    path_annotation = opt['<annotation_sqlite>']
    path_fastq = opt['<unaligned_fastq>']

    root = os.path.basename(os.path.dirname(paths_sqlite[0]))
    path_output = os.path.join(output_dir, f"{root}.aligned.merged.sqlite")

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('calc_mappability_merge.sqlite.sql.j2')

    if os.path.exists(path_output):
        os.remove(path_output)

    for p in paths_sqlite:
        rendered = template.render(
            src_database=p
            )

        root, _ = os.path.splitext(os.path.basename(p))
        path_sql = os.path.join(output_dir, f".{root}.sql")

        with open(path_sql, "w") as f:
            f.write(rendered)

        cmd = f"sqlite3 {path_output} < {path_sql}"
        sys.stderr.write(cmd)
        proc = subprocess.run(cmd, shell=True, capture_output=True)

        sys.stdout.write(proc.stdout.decode())
        sys.stderr.write(proc.stderr.decode())

        for f in [path_sql]:
            os.remove(f)

    root, _ = os.path.splitext(os.path.basename(path_output))
    path_tsv_unaligned = os.path.join(output_dir, f".{root}.unaligned.tsv")
    path_sql = os.path.join(output_dir, f".{root}.sql")

    cmd = f"cat {path_fastq} | grep '@' | cut -f1 | tr -d '@' > {path_tsv_unaligned}"
    _ = subprocess.run(cmd, shell=True, capture_output=False)

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('calc_mappability_agg.sqlite.sql.j2')

    rendered = template.render(
        annotation=path_annotation,
        tsv_unaligned=path_tsv_unaligned
        )

    with open(path_sql, "w") as f:
        f.write(rendered)

    cmd = f"sqlite3 {path_output} < {path_sql}"
    proc = subprocess.run(cmd, shell=True, capture_output=True)

    sys.stdout.write(proc.stdout.decode())
    sys.stderr.write(proc.stderr.decode())

    for f in [path_tsv_unaligned, path_sql]:
        os.remove(f)


def main():
    opt = docopt(__doc__)

    if opt['split']:
        split(opt)

    if opt['calc']:
        calc(opt)

    if opt['agg']:
        agg(opt)


if __name__ == '__main__':
    main()
