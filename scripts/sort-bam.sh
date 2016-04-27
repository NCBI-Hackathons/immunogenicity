#!/bin/bash

TEMP_DIR="/tmp/${1}_sort/"
mkdir -p $TEMP_DIR
samtools sort -T $TEMP_DIR  --reference $GRCH38_PATH -@ 8 $1
rm -rf $TEMP_DIR
