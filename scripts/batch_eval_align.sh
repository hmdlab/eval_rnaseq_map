#! /bin/bash

strandedness="fr"
gtf="shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.gtf"

find results/test01_main -name "*.bam" | \
  sort | \
  grep -v unmapped.bam | \
  grep -v Aligned.toTranscriptome.out.bam | \
  grep -v to_transcriptome.out.bam | \
  xargs -L1 -I{} sh -c "qsub -V -o {}.eval.stdout -e {}.eval.stderr scripts/eval_align_conv.py --strandedness ${strandedness} ${gtf} {}"
