# Configures some environment variables for finding tools

export GATK_PATH="/home/ubuntu/software/GATK/GenomeAnalysisTK.jar"
export HISAT2_PATH="/home/ubuntu/software/hisat2-2.0.3-beta/hisat2"
export PICARD_PATH="/home/ubuntu/software/picard-tools-2.2.2/picard.jar"

export PATH=$PATH:/home/ubuntu/software/hisat2-2.0.3-beta/
export SNPEFF_PATH="/home/ubuntu/software/snpEff/snpEff.jar"
export VEP_PATH="/home/ubuntu/software/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl"

module load samtools
module load python2

export GRCH38_PATH="/home/ubuntu/app_data/Ensembl/grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
export GRCH38_IDX_BASE_PATH="/home/ubuntu/app_data/Ensembl/grch38/genome"
export GRCH38_GFF_PATH="/home/ubuntu/app_data/Ensembl/grch38/Homo_sapiens.GRCh38.84.gff3"

