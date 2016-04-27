#!/bin/bash

samtools mpileup -ugf $GRCH38_PATH ~/test_data/SRR1927992.sorted.bam | bcftools call -vmO z -o test.vcf.gz
