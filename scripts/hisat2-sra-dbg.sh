#!/bin/bash

hisat2 -f -x $GRCH38_IDX_BASE_PATH --sra $1 --rg-id unit1 --no-spliced-alignment \
    --threads 8
