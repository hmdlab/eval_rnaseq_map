#! /bin/bash
#
# SAM to BAM
#
# Usage:
#   find *.sam | xargs -L1 -I{} this.sh {}
#

cmd="qsub -b y -V -l s_vmem=48G -l mem_req=48G -cwd -o ${1/.sam/.bam} -e ${1/.sam/.log} 'samtools view -S -b $1'"
echo $cmd
eval $cmd
