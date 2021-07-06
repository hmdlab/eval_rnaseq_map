#! /bin/bash
#
# UGE wrapper for polyester.R
#
#$ -S /bin/bash
#$ -l s_vmem=128G -l mem_req=128G
#$ -l epyc
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/


script_dir="$HOME/projects/pj04f_eval_rnaseqde_map/scripts"

inputs=($@)

rscript="Rscript"
cmd="${rscript} $script_dir/sim_polyester_from_cntmat.R ${inputs[@]}"
echo $cmd
$cmd
