#/usr/bin/env bash

java -jar $GATK_PATH -T HaplotypeCaller -I $1 -R $GRCH38_PATH -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o output.vcf
