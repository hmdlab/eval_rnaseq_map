#! /bin/bash
#
# Align using bowtie2 (SR) kmers in each annotation map to itself
#
# Usage:
#   qsub -V -t 1-5 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=24G -l mem_req=24G
#$ -pe def_slot 8
#$ -l d_rt=1488:00:00 -l s_rt=1488:00:00
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

align () {
  cmd="bowtie2 \
    -p 4 \
    -x ${index_dir}/$1 \
    -U ${kmer_dir}/$2 \
    --un ${output_dir}/$1.${2/.fastq.gz/.unaligned.fastq} \
    -S ${output_dir}/$1.${2/.fastq.gz/.sam} \
    --norc \
    --sensitive \
    --all \
    --seed 12345"
  echo $cmd
  eval $cmd

  cmd="samtools view -S -b ${output_dir}/$1.${2/.fastq.gz/.sam} > ${output_dir}/$1.${2/.fastq.gz/.bam}"
  echo $cmd
  eval $cmd && rm ${output_dir}/$1.${2/.fastq.gz/.sam}
}

index_dir="assets/references/indexes/bowtie2"
kmer_dir="results/kmers/075"
output_dir="results/mappabilities/075"

file_roots=(
GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted
GCF_000001405.39_GRCh38.p13_transcripts.formatted
gencode.v31.basic.transcripts.formatted
gencode.v31.transcripts.formatted
NONCODEv5_human
)

index=${file_roots[$((SGE_TASK_ID-1))]}
kmer=${file_roots[$((SGE_TASK_ID-1))]}.fastq.gz

align $index $kmer
