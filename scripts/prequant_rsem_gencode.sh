#! /bin/bash
#
# Quantificate expression using RSEM
#
# Usage:
#   find '*.bam' | sort | xargs qsub -t 1-n this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=48G -l mem_req=48G
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/


inputs=($@)
input=${inputs[$((SGE_TASK_ID-1))]}
output_dir="$(basename $(dirname ${input}))_rsem"

if [ ! -e ${output_dir} ]; then
  mkdir -p ${output_dir}
fi

# IMPORTANT: FASTQ files must sorted by name
strands=("--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 0" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1" "--forward-prob 1")
idx=$(((SGE_TASK_ID-1) * 2))

strand=${strands[idx]}

# GENCODE v31 comprehensive
index_prefix="$HOME/share/assets/references/grch38/indexes/rsem/gencode.v31.transcripts.formatted.ercc"

# anno_bam=$2              # STAR alignment to annotation
# paired_end="true"        # "true" if alignment was on paired-end data.
# read_strand="unstranded" # strandedness of read (forward, reverse, unstranded)
rnd_seed=12345             # Random seed.  ENCODE has been using 12345
ncpus=8                    # Number of cpus available.

# paired end, strandness fr
extra_flags="--paired-end"
extra_flags="${extra_flags} ${strand}"

cmd="rsem-calculate-expression \
    --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p $ncpus \
    --no-bam-output --ci-memory 30000 ${extra_flags} \
    ${input} \
    ${index_prefix} \
    ${output_dir}/rsem"

echo ${cmd}
${cmd}
