#! /bin/bash
#
# Merge kmers
#
# Usage:
#  qsub -t 1-5 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=8G -l mem_req=8G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

inputs=(
    gencode.v31.transcripts.formatted.fastq \
    GCF_000001405.39_GRCh38.p13_transcripts.formatted.fastq \
    FANTOM_CAT.lv3_robust.fastq \
    mitranscriptome.v2.fastq \
    NONCODEv5_human.fastq
)

work_dir="results/kmers/100"
input=${inputs[$((SGE_TASK_ID-1))]}
output=${input/.fastq/.mod.fastq}

cmd="cat ${work_dir}/${input} | sed -e 's/^@/@$((SGE_TASK_ID-1))/g' >> ${work_dir}/${output}"
echo $cmd
eval $cmd
