#!/bin/bash

SRA=$1

UNSORTED=`scripts/v002_hisat2-sra.sh $SRA`
scripts/sort-bam.sh $UNSORTED > "$SRA".sorted.bam
scripts/samtools-index.sh "$SRA".sorted.bam
scripts/samtools-call.py  "$SRA".sorted.bam
scripts/annotate-vep.sh "$SRA".vcf

