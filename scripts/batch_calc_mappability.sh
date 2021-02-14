#! /bin/bash
#
# Calculate mappability
#
# Usage:
#   qsub -V -t 1-10 this.sh [subcommand]
#
#$ -S /bin/bash
#$ -l s_vmem=256G -l mem_req=256G
#$ -l d_rt=1488:00:00 -l s_rt=1488:00:00
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


src_dir="results/mappabilities/100"

annotations=(
    shared/assets/references/grch37/annotations/fantomcat/FANTOM_CAT.lv3_robust.sqlite \
    shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted_curated.sqlite \
    shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted.sqlite \
    shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted.sqlite \
    shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13_genomic.formatted.sqlite \
    shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
    shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
    shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
    shared/assets/references/grch38/annotations/gencode_refseq/gencode.v31_refseq.v109.20190607.sqlite \
    shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite \
    shared/assets/references/grch37/annotations/mitranscriptome/mitranscriptome.v2.sqlite \
    shared/assets/references/grch38/annotations/noncode/NONCODEv5_hg38.lncAndGene.sqlite
)

file_roots=(
    FANTOM_CAT.lv3_robust.FANTOM_CAT.lv3_robust \
    GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted.GCF_000001405.39_GRCh38.p13_transcripts.curated.formatted \
    GCF_000001405.39_GRCh38.p13_transcripts.formatted.GCF_000001405.39_GRCh38.p13_transcripts.formatted \
    GCF_000001405.39_GRCh38.p13_pc_transcripts.formatted.GCF_000001405.39_GRCh38.p13_pc_transcripts.formatted \
    GCF_000001405.39_GRCh38.p13_lncRNA_transcripts.formatted.GCF_000001405.39_GRCh38.p13_lncRNA_transcripts.formatted \
    gencode.v31.basic.transcripts.formatted.gencode.v31.basic.transcripts.formatted \
    gencode.v31.lncRNA_transcripts.formatted.gencode.v31.lncRNA_transcripts.formatted \
    gencode.v31.pc_transcripts.formatted.gencode.v31.pc_transcripts.formatted \
    gencode.v31_refseq.v109.20190607.transcripts.formatted.gencode.v31_refseq.v109.20190607.transcripts.formatted \
    gencode.v31.transcripts.formatted.gencode.v31.transcripts.formatted \
    mitranscriptome.v2.mitranscriptome.v2 \
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
