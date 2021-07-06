#! /bin/bash

input_dir=$1
strandedness="fr"
gtf="share/assets/references/grch38/annotations/gencode/gencode.v31.annotation.gtf"

find $input_dir -name "*.bam" | \
  sort | \
  grep -v unmapped.bam | \
  # grep -v Aligned.toTranscriptome.out.bam | \
  # grep -v to_transcriptome.out.bam | \
  xargs -L1 -I{} sh -c "python scripts/eval_align_conv.py --strandedness ${strandedness} ${gtf} {} > {}.eval.stdout"
