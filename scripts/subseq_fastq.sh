#! /bin/bash
#
# Subsequence by seqkit
#
# Usage:
#  find *.fastq | sort | xargs qsub -t 1-n*3 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=64G -l mem_req=64G
#$ -cwd
#$ -o ugelogs/
#$ -e ugelogs/

lengths=(35 50 75)
n_iter=${#lengths[@]}
idx_lengths=$(( (SGE_TASK_ID - 1) % n_iter ))
length=${lengths[$idx_lengths]}

inputs=($@)
idx_inputs=$(( (SGE_TASK_ID - 1) / n_iter ))
input=${inputs[$idx_inputs]}

output_dir="$(dirname $input)/L${length}"
mkdir -p $output_dir

output_file=$(basename $input)
cmd="zcat ${input} | seqkit subseq -r 1:${length} | gzip -c > ${output_dir}/${output_file}"
echo $cmd
eval $cmd
