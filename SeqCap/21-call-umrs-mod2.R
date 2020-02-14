#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-20-7
#R script to call UMRs using a 100bp tile

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
reference_100bp_tiles <- args[1]
sample_to_crunch <- args[2]
annotation_suffix <- args[3]
chrom_sizes_path <- args[4]
coverage_filter_min <- args[5]
site_filter_min <- args[6]
MR_percent <- args[7]
UMR_percent <- args[8]

######## de bug
# args
#reference_100bp_tiles = "~/ws/refseqs/sorghum/sites/Sbicolor_454_v3.0.1_100bp_tiles_zBased_sites_counts.txt"
#sample_to_crunch = "Sbicolor_SRR3286309"

# filter args
#coverage_filter_min = 3
#site_filter_min = 2
#MR_percent = 0.4
#UMR_percent = 0.1
#annotation_suffix = paste0("_mC_domains_II", 
#                           "_cov_",coverage_filter_min, 
#                           "_sites_",site_filter_min, 
#                           "_MR_",MR_percent,
#                           "_UMR_",UMR_percent)
########
###########################
library(tidyverse)
library(ggthemes)
# library("seqinr")
old.scipen <- getOption("scipen")
options(scipen=999)
# library(wesanderson)
library(RColorBrewer)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

###########################
# 
###########################

folder_prefix = paste0(sample_to_crunch, annotation_suffix)

dir.create(folder_prefix, showWarnings = F)

out_dir = paste0(folder_prefix, "/mC_UMT_annotation")
dir.create(out_dir, showWarnings = F)

out_dir_beds = paste0(folder_prefix, "/mC_UMT_annotation_beds")
dir.create(out_dir_beds, showWarnings = F)

######

#############
# folders
umr_out_dir = paste0(folder_prefix, "/mC_UMR_annotation")
dir.create(out_dir, showWarnings = F)

umr_out_dir_beds = paste0(folder_prefix, "/mC_UMR_annotation_beds")
dir.create(out_dir_beds, showWarnings = F)

#############
## Module #2
#############

# read in ref
mC_domains_merge <- read_tsv(paste0(out_dir_beds,"/", sample_to_crunch, "_UMTs_sorted_merge.bed"), col_names =  c("chr","start", "end"))

mC_domains_merge <- mC_domains_merge %>%
  mutate(type = "UMT",
         size = end-start,
         location = paste0(chr, ":", start, ":", end))

UMT_summary <- mC_domains_merge %>% summarise(total_merged_UMT = n(), total_MB = sum(size)/1000000)
UMT_summary
write.table(UMT_summary, paste0(out_dir, "/", sample_to_crunch, "_UMTs_merge_total_MB_summary.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

write.table(mC_domains_merge, paste0(out_dir_beds, "/", sample_to_crunch, "_UMTs_merge_size.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

