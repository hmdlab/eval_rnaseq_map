#! /bin/bash
#
# Calculate mappability
#
# Usage:
#   qsub -V -t 1-5 this.sh [subcommand]
#
#$ -S /bin/bash
#$ -l s_vmem=256G -l mem_req=256G
#$ -l d_rt=480:00:00 -l s_rt=480:00:00
#$ -l intel
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/


split () {
  cmd="python scripts/calc_mappability.py split \
    --output-dir ${src_dir} \
    --template-dir ${template_dir} \
    $1"
  echo $cmd
  eval $cmd
}


# NOTE: Run on login node -> n x qsub
calc () {
  cmd="qsub -V scripts/calc_mappability.py calc \
    --output-dir $(dirname $2) \
    --template-dir ${template_dir} \
    $1 \
    $2"
  echo $cmd
  eval $cmd
}


agg () {
  cmd="python scripts/calc_mappability.py agg \
    --output-dir $1 \
    --template-dir ${template_dir} \
    $2 \
    $3 \
    $4"
  echo $cmd
  eval $cmd
}

src_dir=$2

annotations=(
share/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted_curated.sqlite \
share/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted.sqlite \
share/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
share/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
share/assets/references/grch38/annotations/noncode/NONCODEv5_hg38.lncAndGene.sqlite
)

file_roots=(
GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted.GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted \
GCF_000001405.39_GRCh38.p13_transcripts.formatted.GCF_000001405.39_GRCh38.p13_transcripts.formatted \
gencode.v31.basic.transcripts.formatted.gencode.v31.basic.transcripts.formatted \
gencode.v31.transcripts.formatted.gencode.v31.transcripts.formatted \
NONCODEv5_human.NONCODEv5_human
)

annotation=${annotations[$((SGE_TASK_ID-1))]}
root=${src_dir}/${file_roots[$((SGE_TASK_ID-1))]}
alignment=${src_dir}/${file_roots[$((SGE_TASK_ID-1))]}.bam
unaligned_read=${src_dir}/${file_roots[$((SGE_TASK_ID-1))]}.unaligned.fastq
template_dir=scripts


if [ -z $1 ]; then
  echo No arguments specified
  exit 1
fi

if [ $1 = 'split' ]; then
  split $alignment
fi

if [ $1 = 'calc' ]; then
  len=$((${#file_roots[@]} ))
  for (( i=0; i<${#file_roots[@]}; i++ ));
  do
    inputs=$(find ${src_dir}/${file_roots[$i]} -path "*.aligned.*.tsv" | sort)
    for j in ${inputs[@]}
    do
      calc ${annotations[$i]} ${j}
    done
  done
fi

if [ $1 = 'agg' ]; then
  inputs=$(find ${root} -path "*.aligned.*.sqlite" | sort)
  agg $src_dir $annotation $unaligned_read "${inputs[@]}"
fi
