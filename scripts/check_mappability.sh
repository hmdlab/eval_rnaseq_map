#! /bin/bash

input=$1

cat $input | q -d'|' "select (select sum(c2) from - where c1 = 'ext.transcript_alignment_truths') as c1, (select max(c2) from - where c1 = 'transcript_alignment_truths;after') as c2, (select sum(c2) from - where c1 = 'ext.gene_alignment_truths') as c3, (select max(c2) from - where c1 = 'gene_alignment_truths;after') as c4;"
