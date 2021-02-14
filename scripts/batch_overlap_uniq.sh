#! /bin/bash
#
# Align using bowtie2 (SR) merged kmers map to other transcriptome
#
# Usage:
#   qsub -V -t 1-5 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=128G -l mem_req=128G
#$ -pe def_slot 1
#$ -l d_rt=1488:00:00 -l s_rt=1488:00:00
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/


unique () {
    cmd="python scripts/overlap_kmers.py uniq \
    $1"
    echo $cmd
    eval $cmd
}

kmer_dir="results/kmers/100"
output_dir="results/mappabilities/100"

file_roots=(
    FANTOM_CAT.lv3_robust
    GCF_000001405.39_GRCh38.p13_transcripts.formatted
    gencode.v31.transcripts.formatted
    mitranscriptome.v2
    NONCODEv5_human
)

kmer="${kmer_dir}/${file_roots[$((SGE_TASK_ID-1))]}.fastq.gz"

unique $kmer
