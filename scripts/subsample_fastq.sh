#! /bin/bash
#
# Subsampling and subsequence by seqkit
#
# Usage:
#  find *.fastq | sort | xargs qsub -t 1-n*4 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=128G -l mem_req=128G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

nreads=(5000000 10000000 20000000 40000000)
n_iter=${#nreads[@]}
idx_nreads=$(( (SGE_TASK_ID - 1) % n_iter ))
nread=${nreads[$idx_nreads]}

inputs=($@)
idx_inputs=$(( (SGE_TASK_ID - 1) / n_iter ))
input=${inputs[$idx_inputs]}

output_dir="$(dirname $input)/D$(( nread / 1000000 ))M"
mkdir -p $output_dir

output_file=$(basename $input)
seed=12345
cmd="seqkit sample -n ${nread} -s ${seed} -o ${output_dir}/${output_file} ${input}"
echo $cmd
eval $cmd
