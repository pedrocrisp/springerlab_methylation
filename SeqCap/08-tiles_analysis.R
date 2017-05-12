#!/usr/bin/Rscript
##########

# Peter Crisp
# 2017-04-08
# R script to amend 100bp tile data with incorrect chr end tiles
# because 100bp runs past the end of the chr in most cases

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
coverage_filter <- args[3]
coverage_filter

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

#########
# CHH analysis
#sample = "F1-16_Index15_S5"
#coverage_filter = 2

#load a file
mC_tile_table <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.txt"), col_names = TRUE,
                         cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  C = col_integer(),
  CT = col_integer(),
  ratio = col_double(),
  sites_with_data = col_integer(),
  chh_sites = col_integer()
))

# coverage filter
mC_tile_table_CHHcov <-
mC_tile_table %>%
  unite(tile, chr, start,end, sep=":") %>%
  mutate(CHH_cov= CT/chh_sites) %>%
  select(tile, CHH_cov) %>%
  filter(CHH_cov >= coverage_filter)

mC_tile_table_CHHcov

paste0("Total tiles in genome 21,350,963")
paste0("For sample ", sample, " there are ", length(mC_tile_table_CHHcov$tile), " tiles with >= 2x CHH coverage")

write.table(mC_tile_table_CHHcov, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH_cov.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
