#!/bin/bash

java -jar $PICARD_PATH AddOrReplaceReadGroups RGLB=LaneX RGPU=NONE RGSM=AnySampleName I=$1 O=$2 RGPL=illumina
