#!/bin/bash


hisat2 -f -x $GRCH38_IDX_BASE_PATH --sra $1 --rg-id unit1 --no-spliced-alignment \
    --threads 8 | samtools view -T $GRCH38_PATH -hb > "$1".unsorted.bam
