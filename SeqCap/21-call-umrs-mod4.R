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
## Module #4
#############

# read in ref
ND_overlaps <- read_tsv(paste0(out_dir_beds, "/", sample_to_crunch, "_UMTs_merge_NDs.bed"), col_names =  c("chr","start", "end", "features", "sizes", "locations"),
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
# about 20K over 34%

################## ################## ##################
################## ################## ##################
# run this chunk the forst time through to get black list

# filter to tile >34% ND
over_34_pct <- ND_overlaps_sizes_pct  %>%
  mutate(cat = ifelse(pct_ND == 0 , 0, ifelse(pct_ND > 34, ">34", "<=34"))) %>% 
  filter(cat == ">34") %>% 
  select(chr, start, end, cat) %>%
  right_join(ND_overlaps_sizes, by = c("chr", "start", "end")) %>%
  filter(cat == ">34") 
over_34_pct

# number of concatonated tiles
over_34_pct %>% group_by(chr, start, end) %>% summarise(n = n()) %>% group_by(n) %>% summarise(freq = n())

# plot size distributin of tiles
g <- over_34_pct %>% ungroup() %>% filter(features == "ND") %>% filter(as.double(sizes) < 2000)  %>%
  ggplot(., aes(as.double(sizes))) +
  geom_density() +
  theme_minimal() +
  text_size_theme_8

ggsave(plot = g, filename = paste0(out_dir, "/ND_inbetweens_size_density_over_34.pdf"), h = 3, w = 3)

# extract list of ND tiles that contribute to merges that are >34% to get blacklist

black_list <- over_34_pct %>% 
  filter(features == "ND") %>% # 2,120
  select(locations, cat) %>%
  separate(locations, into = c("chr", "start", "end"), sep = ":")

write.table(black_list,  paste0(out_dir_beds,"/", sample_to_crunch, "_NDs_between_UMTs_black_list.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

################## ################## ##################
################## ################## ##################

# # add sizes
# # merge back in features metadata
# ND_overlaps_sizes_pct_meta <- ND_overlaps_sizes_pct %>%
#   left_join(features_metadata, by = c("chr", "start", "end"))
# 
# ND_overlaps_sizes_pct %>% distinct(chr, start, end)
# 
# # write output
# write.table(ND_overlaps_sizes_pct_meta, "B73L_mC_domains_v1_annotation_v0_4_UMTs_merge_NDs_filtered.bed", sep = "\t", quote = F, row.names = F, col.names = F)

# Remove black listed regions from the original ND merge
# read in ref
ND_overlaps <- read_tsv(paste0(out_dir_beds, "/NDs_Olap_UMTs.bed"), col_names =  c("chr","start", "end","B_chr", "b_start", "b_end", "b_type","b_size", "b_location", "distance"), 
                        cols(b_size = col_character()))

black_list <- read_tsv(paste0(out_dir_beds,"/", sample_to_crunch, "_NDs_between_UMTs_black_list.bed"), col_names =  c("chr","start", "end","cat"))

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
  left_join(black_list, by = c("chr","start", "end")) %>%
  ungroup() %>%
  mutate(cat = ifelse(is.na(cat), "ok", cat)) %>%
  filter(!cat == ">34") %>%
  mutate(name = "ND", size = end-start,
         location = paste0(chr, ":", start, ":", end)) %>%
  select(-cat)

max(ND_between_UMTs$size)
# 100?

print(paste0("ND between UMRs passed pct filtering to be merged with UMRs: #"))
ND_between_UMTs
# 55477

print(paste0("ND between UMRs passed pct filtering to be merged with UMRs: MB"))
ND_between_UMTs %>% summarise(MB = sum(size)/1000000)

ND_between_UMTs_summary <- ND_between_UMTs %>% summarise(ND_passed_filtering_for_merging = n(), MB = sum(size)/1000000)
ND_between_UMTs_summary
write.table(ND_between_UMTs_summary, paste0(out_dir, "/", sample_to_crunch, "_NDs_between_UMTs_pct_filtered_summary.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

write.table(ND_between_UMTs, paste0(out_dir_beds, "/", sample_to_crunch, "_NDs_between_UMTs_pct_filtered.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

