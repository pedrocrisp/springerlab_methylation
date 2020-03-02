#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e
#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced.
set -x

# This script overlaps a umr file with an ACR file
# UMRs are annotated as inaccessible iUMR and accessible aUMRs and the total of each reported in the summary file
# ACRs are annotated as methylated (mACR) and unmethylated (uACRs) and the total of each reported in the summary file
# Overlap include min 1bp overlap

### Inputs
# The input bed files should be 6 columns
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: name (or other desired info) eg UMR or ACR
# Column5: optional eg ID eg UMR_21
# Column6: optional eg size

# Output
# Column1: chromosome
# Column2: start
# Column3: stop
# Column4: name (or other desired info) eg UMR or ACR
# Column5: optional eg ID eg UMR_21
# Column6: optional eg size
# Column7: annotation eg i/aUMR or m/uACR

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
# sample_name=$1
umr_bed_file=$1
acr_bed_file=$2
chrom_sizes_file=$3
summary_output_folder=$4

# debug
# cd ~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation_beds
# sample_name=Sbicolor_SRR3286309
# umr_bed_file=Sbicolor_SRR3286309_UMRs_6col.bed
# acr_bed_file=~/ws/analysis/ongoing/mC_sorghum/Sorghum_7days_leaf_ACRs.bed
# chrom_sizes_file=~/ws/refseqs/sorghum/Sbicolor_454_v3.0.1_sorted.chrom.sizes
# summary_output_folder=~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation

tmp_umr_sample_prefix="$(basename $umr_bed_file)"
umr_sample_prefix="${tmp_umr_sample_prefix%%.*}"

tmp_acr_sample_prefix="$(basename $acr_bed_file)"
acr_sample_prefix="${tmp_acr_sample_prefix%%.*}"

# optput file names
umr_outputFile=${umr_sample_prefix}_olap_${acr_sample_prefix}
echo $umr_outputFile
acr_outputFile=${acr_sample_prefix}_olap_${umr_sample_prefix}
echo $acr_outputFile

#basename
#umr_sample_prefix="$(basename $umr_bed_file)"
#echo $umr_sample_prefix

#umr_outputFile="${umr_sample_prefix%%.*}_access"
#echo $umr_outputFile

#acr_sample_prefix="$(basename $acr_bed_file)"
#echo $acr_sample_prefix

#acr_outputFile="${acr_sample_prefix%%.*}_meth"
#echo $acr_outputFile

mkdir -p logs

########## Modules #################

module load R/3.3.2
module load bedtools

########## overlaps ACR annotation #################

# ACRs regions that overlap UMRs
bedtools closest \
-a $acr_bed_file \
-b $umr_bed_file \
-mdb all \
-t all \
-D b \
-g $chrom_sizes_file \
> tmp_${acr_outputFile}.bed

########## R module to parse ACR annotation output #################
R -f ~/gitrepos/springerlab_methylation/utils/91-overlap-UMR-ACR-R1.R \
--args $acr_outputFile $summary_output_folder

########## overlaps UMR annotation #################

# UMRs that overlap ACRs
bedtools closest \
-a $umr_bed_file \
-b $acr_bed_file \
-mdb all \
-t all \
-D b \
-g $chrom_sizes_file \
> tmp_${umr_outputFile}.bed

########## R module to parse UMR annotation output #################
R -f ~/gitrepos/springerlab_methylation/utils/91-overlap-UMR-ACR-R2.R \
--args $umr_outputFile $summary_output_folder

# to run
# bash ~/gitrepos/springerlab_methylation/utils/91-overlap-UMR-ACR.sh \
# Sbicolor_SRR3286309_UMRs_6col.bed \
# ~/ws/analysis/ongoing/mC_sorghum/Sorghum_7days_leaf_ACRs.bed \
# ~/ws/refseqs/sorghum/Sbicolor_454_v3.0.1_sorted.chrom.sizes \
# ~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation
# >> \
# logs/$(date +%Y%m%d-%H%M%S)_UMR_ACR_overlap.log 2>&1
