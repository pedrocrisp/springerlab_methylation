#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-02-08
#R script to overlap UMRs files

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
# sample_prefix <- "Sbicolor_SRR3286309_UMRs_6col"
# inputFile <- "Sbicolor_SRR3286309_UMRs_6col_access.bed"
# summary_output_folder="~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation"
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
                     col_names = c("chr", "start", "end", "feature", "ID", "optional",
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
  mutate(accessibility = ifelse(distance == 0, "shared_UMR", "unique_UMR")) %>%
  select(chr, start, end, feature, ID, optional, accessibility) %>%
  distinct()
# 32,620

overlaps_distinct_collapsed %>% group_by(chr, start) %>% summarise(n = n()) %>% group_by(n) %>% summarise(n2 = n())

overlaps_distinct_collapsed
# 21,267

write.table(overlaps_distinct_collapsed, paste0(inputFile, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed %>%
  group_by(accessibility) %>%
  summarise(n = n()) %>%
  mutate(percentage = n/sum(n)*100)

overlaps_distinct_collapsed_filtered_summary
# # A tibble: 2 x 3
#   methylation     n percentage
#   <chr>       <int>      <dbl>
# 1 non_UMR      1804       5.62
# 2 UMR         30307      94.4

write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(summary_output_folder, "/", inputFile, "_summary.tsv"))
