#!/bin/bash


hisat2 -x $GRCH38_IDX_BASE_PATH --sra $1 --no-spliced-alignment --threads 8 | samtools view -bS > "$1".unsorted.bam

echo "$1".unsorted.bam
