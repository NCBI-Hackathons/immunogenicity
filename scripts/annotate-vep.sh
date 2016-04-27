#!/bin/bash


perl $VEP_PATH -i $1 --format vcf --vcf --species homo_sapiens_refseq
