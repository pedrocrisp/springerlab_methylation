#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-20-7
#R script to annotate UMRs using an ACR file

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
# inputFile <- "generations_merged_UMRs_access_6col_olap_DAP_seq_27_TFs_Galli2018_Ricci2019"
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
                     col_names = c("chr", "start", "end", "feature", "ID", "UMR_type",
                                   "b_chr", "b_start", "b_end", "b_summit", "b_location", "b_cat", "distance"), 
                     cols(chr = col_character(),
                          b_summit = col_character(), 
                          b_chr = col_character()))

overlaps 
# 34,514 (acrs 32,481) - only a few extra multi overlaps to resolve!

# resovle multi-overlaps
# I think most of the multi-overlaps will be at distance zero - this pipeline is mainly aimed at resolving that. Its possible some will be due to a peak being exactally the same length between an upstrea and downstream gene, but I might just pick one for these because I think they will be rare
# get distinct tile + feature rows
overlaps_distinct_collapsed <- overlaps %>%
  mutate(overlap = ifelse(distance == 0, UMR_type, "no_overlap")) %>%
  select(chr, start, end, feature, ID, UMR_type, overlap) %>%
  distinct()
# 32,620 

overlaps_distinct_collapsed %>% group_by(chr, start) %>% summarise(n = n()) %>% group_by(n) %>% summarise(n2 = n())

overlaps_distinct_collapsed
# 21,267

write.table(overlaps_distinct_collapsed, paste0(inputFile, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed %>% 
  group_by(overlap, UMR_type) %>% 
  summarise(n = n()) %>%
  group_by(UMR_type) %>%
  mutate(percentage = n/sum(n)*100)

overlaps_distinct_collapsed_filtered_summary
# # A tibble: 2 x 3
#   overlap     n percentage
#   <chr>       <int>      <dbl>
# 1 aUMR      21076       19.5
# 2 iUMR         57178      52.9
# 3 no_overlap  29761     27.6

write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(summary_output_folder, "/", inputFile, "_summary.tsv"))


