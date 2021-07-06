#! /bin/bash
#
# Resample reads using seqkit
#
# Usage:
#   find '*.fastq' | sort | xargs qsub -t 1-n this.sh
#
#$ -S /bin/bash
#$ -l month -l medium
#$ -l s_vmem=256G -l mem_req=256G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/

inputs=($@)

num_reads=(5000000 10000000 20000000 40000000)
seed=12345

input=${inputs[$((SGE_TASK_ID-1))]}
output_file=$(basename $input)

for n in ${num_reads[@]}
do
  if [ ! -e ${n} ]; then
    mkdir -p ${n}
  fi

  cmd="zcat ${input} | seqkit sample -n ${n} -s ${seed} -o ./${n}/${output_file}"
  echo $cmd
  eval $cmd
done
