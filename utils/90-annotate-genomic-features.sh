#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e
#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced.
set -x

# This script annotates a bed file with genomic features
# Only gene are overlapped in this version to get an accurate gene-distal list
# Make three categories
# 1. genic
# 2. gene_proximal < 2kb (ths is captured by 1 or 2kb up/down categories)
# 3. gene distal > 2kb (anything called intergenic)

### Inputs
# The input bed file should be 6 columns
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: name (or other desired info) eg UMR or ACR
# Column5: optional eg ID eg UMR_21
# Column6: optional eg size

# The input reference file should be 6 columns
# The input file should be parsed eg from a GFF3 to only include the genes/genomic features that are of interest
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: name of genomic feature (assumed to be a typeof gene)
# Column5: score
# Column6: strand
# Column7: gene ID

# Output
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: classification genic/proximal/distal annotation
# Column5: score
# Column6: strand

###
#code to make script work on both osx and linux.

if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

########## Args #################
sample_name=$1
bed_file=$2
feature_type=$3
reference_annotation=$4
reference_annotation_mame=$5
chrom_sizes_file=$6
summary_output_folder=$7

# debug
# cd ~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation_beds
# sample_name=Sbicolor_SRR3286309
# bed_file=Sbicolor_SRR3286309_UMRs_6col.bed
# feature_type=gene
# reference_annotation=~/ws/refseqs/sorghum/annotation/Sbicolor_454_v3.1.1.genes_sorted.bed
# reference_annotation_mame=gene_v1
# chrom_sizes_file=~/ws/refseqs/sorghum/Sbicolor_454_v3.0.1_sorted.chrom.sizes
# summary_output_folder=~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation

#basename
sample_prefix="$(basename $bed_file)"
echo $sample_prefix

outputFile="${sample_prefix%%.*}_olap_${reference_annotation_mame}.bed"
echo $outputFile

mkdir -p logs

########## Modules #################

module load R/3.3.2
module load bedtools

########## overlaps #################

bedtools closest \
-a $bed_file \
-b $reference_annotation \
-mdb all \
-t all \
-D b \
-g $chrom_sizes_file \
> $outputFile


########## R module to parse overlaps #################
R -f ~/gitrepos/springerlab_methylation/utils/90-annotate-genomic-features.R \
--args $sample_name $outputFile $reference_annotation_mame $summary_output_folder

# to run
# bash ~/gitrepos/springerlab_methylation/utils/90-annotate-genomic-features.sh \
# <>.bed \
# >> \
# logs/$(date +%Y%m%d-%H%M%S)_make_sites.log 2>&1

# to run a batch
# find . -name "*.bed." | parallel -j 20 bash ~/gitrepos/springerlab_methylation/utils/context_sites.sh {} >> $(date +%Y%m%d-%H%M%S)_make_sites.log 2>&1
