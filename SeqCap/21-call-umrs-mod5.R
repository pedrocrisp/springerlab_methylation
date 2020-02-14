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
## Module #5
#############

# read in ref
ND_overlaps <- read_tsv(paste0(out_dir_beds, "/", sample_to_crunch, "_UMTs_merge_NDs_pct_filtered.bed"), col_names =  c("chr","start", "end", "features", "sizes", "locations"),
                        cols(sizes = col_character()))

ND_overlaps_sizes <- ND_overlaps %>% 
  mutate(features = strsplit(as.character(features), ","),
         sizes = strsplit(as.character(sizes), ","),
         locations = strsplit(as.character(locations), ",")) %>% 
  unnest(features, sizes, locations)

ND_overlaps_sizes

# check size distribution before merging

size_distro <- ND_overlaps_sizes %>%
  filter(!features == "ND") %>%
  mutate(size_cat = ifelse(as.double(sizes) < 300, "small", 
                           ifelse(as.double(sizes) >=300 & as.double(sizes) <900, "med", "large"))) 
#size summary
size_distro %>% group_by(size_cat) %>% summarise(n = n()) %>% mutate(pct = n/sum(n)*100)

# how often are the mergers longer than 3?

ND_overlaps_sizes %>% group_by(chr, start, end) %>% summarise(n = n()) %>% group_by(n) %>% summarise(freq = n())

ND_overlaps_sizes %>% group_by(chr, start, end) %>% summarise(n = n()) %>% group_by(n) %>% summarise(freq = n()) %>% filter(!n %in% c(1, 3)) %>% summarise(total = sum(freq))

# features_metadata <- ND_overlaps_sizes %>% 
#   mutate(cat = ifelse(features == "ND", "ND", "UMR")) %>%
#   filter(cat == "UMR") %>%
#   group_by(chr, start, end) %>%
#   summarise(features = paste(features, collapse = "|"))

# calculate percent ND
ND_overlaps_sizes_pct <- ND_overlaps_sizes %>% 
  mutate(cat = ifelse(features == "ND", "ND", "UMR")) %>%
  group_by(chr, start, end, cat) %>% 
  summarise(total_size = sum(as.double(sizes))) %>%
  ungroup() %>%
  spread(key = cat, value = total_size) %>%
  replace(., is.na(.), 0) %>%
  mutate(pct_ND = ND/(ND + UMR) *100)

ND_overlaps_sizes_pct

# summarise percent of UMR tiles
ND_overlaps_sizes_pct_summary <- ND_overlaps_sizes_pct %>%
  mutate(cat = ifelse(pct_ND == 0 , 0, ifelse(pct_ND > 34, ">34", "<=34"))) %>%
  group_by(cat) %>%
  summarise(n = n())

ND_overlaps_sizes_pct_summary
# none over 34 now


# add sizes
ND_overlaps_sizes_pct_size <- ND_overlaps_sizes_pct %>%
  mutate(size = end - start) %>% 
  mutate(size_cat = ifelse(size < 300, "small", 
                           ifelse(size >=300 & size <900, "med", "large"))) %>%
  mutate(UMT_cat = ifelse(size_cat %in% c("large", "med"), "UMR", "small_UMT")) %>%
  group_by(UMT_cat) %>%
  mutate(ID = ifelse(size_cat == "small", "UMT", paste("UMR", 1:n(), sep = "_")))

#size summary
ND_overlaps_sizes_pct_size %>% group_by(size_cat) %>% summarise(n = n()) %>% mutate(pct = n/sum(n)*100)
# 56% small

ND_overlaps_sizes_pct_size %>% filter(size_cat == "small") %>% summarise(n = n())

# write output
write.table(ND_overlaps_sizes_pct_size, paste0(out_dir_beds, "/", sample_to_crunch, "_UMTs_merge_NDs_filtered_size.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# total UMR MB calculation

summary_out <- ND_overlaps_sizes_pct_size %>% group_by(size_cat) %>% summarise(MB = sum(size)/1000000)
write.table(summary_out, paste0(out_dir, "/", sample_to_crunch, "_UMTs_merge_NDs_filtered_size_summary.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

summary_out
# large - 86.8 (maize 95.9 MB)
# med - 36.1 (maize 27.2 MB)
# small - 18.0

ND_overlaps_sizes_pct_size %>% filter(!size_cat == "small") %>% summarise(MB = sum(size)/1000000)
# 118

############################
# UMR only - size filtered
UMR_only_out <- ND_overlaps_sizes_pct_size %>%
  filter(UMT_cat == "UMR")
# 112,872

# write
write.table(UMR_only_out, paste0(out_dir_beds, "/", sample_to_crunch, "_UMRs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

# 6 column version - chr, start, end, data type (umr/acr), ID, other meta data (eg category)
UMR_only_out <- ND_overlaps_sizes_pct_size %>%
  filter(UMT_cat == "UMR") %>%
  select(chr, start, end, UMT_cat, ID, size_cat)

# write
write.table(UMR_only_out,  paste0(out_dir_beds, "/", sample_to_crunch, "_UMRs_6col.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

