#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-20-7
#R script to call the proportion of regions of interest (eg ChIP peaks) that overlap i or aUMRs
# Where a region overlaps both it will be classified as an aUMR because it was already partially open

# Notes
# post filtering may be required to remove organelles or other undesired contigs

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
# sample_name <- args[1]
inputFile <- args[1]
summary_output_folder <- args[2]

######## de bug
# args
# inputFile <- "DAP_seq_27_TFs_Galli2018_Ricci2019_olap_generations_merged_UMRs_access_6col"
# summary_output_folder="~/ws/analysis/ongoing/UMRs_maize/analysis/generations_merged_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation"
# 

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
# Module #1
###########################

overlaps <- read_tsv(paste0("tmp_", inputFile, ".bed"),
                     col_names = c("chr", "start", "end", "summit", "location", "cat",
                                   "b_chr", "b_start", "b_end", "b_feature", "b_ID", "UMR_type", "distance"), 
                     cols(chr = col_character(),
                          # Olap_file = col_character(), 
                          b_chr = col_character()))

overlaps %>% print(n = 10)
# 34,514 (acrs 32,481) - only a few extra multi overlaps to resolve!

# resovle multi-overlaps
# I think most of the multi-overlaps will be at distance zero - this pipeline is mainly aimed at resolving that. Its possible some will be due to a peak being exactally the same length between an upstrea and downstream gene, but I might just pick one for these because I think they will be rare
# get distinct tile + feature rows
overlaps_distinct_collapsed <- overlaps %>%
  mutate(methylation = ifelse(distance == 0, UMR_type, "not_UMR")) %>%
  select(chr, start, end, summit, location, cat, methylation) %>%
  distinct() %>%
  group_by(chr, start, end) %>%
  summarise(methylation2 = paste(methylation, collapse = "_"))
  distinct()
# methylation2 
  
  overlaps %>% group_by(chr, start, end) %>% summarise(n = n()) %>% group_by(n) %>% summarise(n2 = n())
  # about 6,000 regions out of 2.15 M are in here twice... FYI. I assue this is just by chance
  
# if overlap includes both iUMR(s) and aUMR(s) then call it aUMR
overlaps_distinct_collapsed_classified <- overlaps_distinct_collapsed %>%
  mutate(methylation_type = ifelse(grepl("aUMR", methylation2), "aUMR",
                                   ifelse(grepl("iUMR", methylation2), "iUMR", 
                                          ifelse(grepl("not_UMR", methylation2), "not_UMR", methylation2))))

overlaps_distinct_collapsed_classified %>% group_by(chr, start) %>% summarise(n = n()) %>% group_by(n) %>% summarise(n2 = n())

overlaps_distinct_collapsed_classified
# 21,267

write.table(overlaps_distinct_collapsed_classified, paste0(inputFile, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_classified %>% 
  group_by(methylation_type) %>% 
  summarise(n = n()) %>%
  mutate(percentage = n/sum(n)*100)

overlaps_distinct_collapsed_filtered_summary
# # A tibble: 2 x 3
#   methylation     n percentage
#   <chr>       <int>      <dbl>
# 1 non_UMR      1804       5.62
# 2 UMR         30307      94.4 

write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(summary_output_folder, "/", inputFile, "_summary.tsv"))

