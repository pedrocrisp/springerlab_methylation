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
## Module #3
#############

# read in ref
# read b_size as ccharacter because there are some dots
ND_overlaps <- read_tsv(paste0(out_dir_beds, "/NDs_Olap_UMTs.bed"), col_names =  c("chr","start", "end","B_chr", "b_start", "b_end", "b_type","b_size", "b_location", "distance"), 
                        cols(b_size = col_character()))

ND_overlaps %>% mutate(name = "ND", size = end-start) %>% summarise(max(size))
# 4.8 MB!

ND_overlaps <- ND_overlaps %>%
  filter(distance == 1) %>%
  group_by(chr, start, end) %>%
  mutate(hits = n())

ND_overlaps %>% group_by(hits) %>% summarise(n_hits = n())

# 1 = 107,967 (maize, 97,422)
# 2 = 196,624 (maize 73,890)
# So 107,967 of the merged ND tiles are adjacent to a UMT; while ~100k (x/2) are inbetween adjacent UMTs
# double maize on inbetweeners

ND_between_UMTs <- ND_overlaps %>%
  filter(hits == 2) %>%
  distinct(chr, start, end) %>%
  ungroup() %>%
  mutate(name = "ND", size = end-start,
         location = paste0(chr, ":", start, ":", end))

max(ND_between_UMTs$size)
# 100?

g <- ND_between_UMTs %>% filter(size < 1000) %>%
  ggplot(., aes(size)) +
  geom_density() +
  theme_minimal() +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/ND_inbetweens_size_density.pdf"), h = 3, w = 3)

print(paste0("ND between UMRs to be tested for merging with UMRs: #"))
ND_between_UMTs
# 55477

print(paste0("ND between UMRs to be tested for merging with UMRs: MB"))
ND_between_UMTs %>% summarise(MB = sum(size)/1000000)

ND_between_UMTs_summary <- ND_between_UMTs %>% summarise(total_merged_NDs_for_testing = n(), MB = sum(size)/1000000)
ND_between_UMTs_summary
write.table(ND_between_UMTs_summary, paste0(out_dir, "/", sample_to_crunch, "_NDs_between_UMTs_summary.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

write.table(ND_between_UMTs, paste0(out_dir_beds, "/", sample_to_crunch, "_NDs_between_UMTs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

