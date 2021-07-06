#! /bin/bash
#
# Build bowtie2 index
#
# Usage:
#   qsub -V -t 1-10 this.sh
#
#$ -S /bin/bash
#$ -l s_vmem=32G -l mem_req=32G
#$ -l d_rt=1488:00:00 -l s_rt=1488:00:00
#$ -cwd
#$ -o ./ugelogs/
#$ -e ./ugelogs/


build () {
  cmd="bowtie2-build -f ${SEQ_DIR}/$1 ${INDEX_DIR}/${1/.fa/}"
  echo $cmd
  $cmd
}

transcriptomes=(
gencode.v31.basic.transcripts.formatted.fa
gencode.v31.transcripts.formatted.fa
gencode.v31_refseq.v109.20190607.transcripts.formatted.fa
FANTOM_CAT.lv3_robust.fa
GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted.fa
GCF_000001405.39_GRCh38.p13_transcripts.formatted.fa
GCF_000001405.39_GRCh38.p13_pc_transcripts.formatted.fa
GCF_000001405.39_GRCh38.p13_lncRNA_transcripts.formatted.fa
mitranscriptome.v2.fa
NONCODEv5_human_transcripts_formatted.fa
gencode.v31.lncRNA_transcripts.formatted.fa
gencode.v31.pc_transcripts.formatted.fa
)

SEQ_DIR="assets/references/sequences"
INDEX_DIR="assets/references/indexes/bowtie2"

for t in "${transcriptomes[@]}"
do
  build $t
done
