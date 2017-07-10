#!/usr/bin/Rscript
##########

# Peter Crisp
# 2017-04-08
# R script to make CHH coverage files for summary stats analysis

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
tile_list <- args[3]
tile_list
output_folder <- args[4]
output_folder

#debugging
#sample = "F4-16_F5-16_Index2_S35"
#data_folder = "analysis/tiles"
#tile_list = "analysis_02_tiles_SeqCap_meta_140_samples/chh_2x_cov_80_sample_tile_list.tsv"
#output_folder = "analysis/tiles_filtered"

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

#########
#list of tiles to use for filtering
# read_tsv(tile_list)
tile_list <- read_tsv(tile_list)
tile_list

#########
# filter CHH tile file
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

# tile filter
mC_tile_table_filtered <-
mC_tile_table %>%
  unite(tile, chr, start,end, sep=":") %>%
  mutate(cov= CT/chh_sites) %>%
  inner_join(tile_list, "tile") %>%
  select(tile, ratio, chh_sites, cov)

mC_tile_table_filtered

write.table(mC_tile_table_filtered, paste0(output_folder, "/", sample, "_BSMAP_out.txt.100.CHH_filtered.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

#########
# filter CG tile file
mC_tile_table <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.txt"), col_names = TRUE,
                         cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  C = col_integer(),
  CT = col_integer(),
  ratio = col_double(),
  sites_with_data = col_integer(),
  cg_sites = col_integer()
))

# tile filter
# also here there is an additional CG filter requiring there to be at least 6 CG sites
# CHECK IF THIS SHOULD ALSO REQUIRE CG COVERAGE FILTER
mC_tile_table_filtered <-
mC_tile_table %>%
  unite(tile, chr, start,end, sep=":") %>%
  mutate(cov= CT/cg_sites) %>%
  inner_join(tile_list, "tile") %>%
  select(tile, ratio, cg_sites, cov) %>%
  filter(cg_sites >= 6)

mC_tile_table_filtered

write.table(mC_tile_table_filtered, paste0(output_folder, "/", sample, "_BSMAP_out.txt.100.CG_filtered.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

#########
# filter CHG tile file
mC_tile_table <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.txt"), col_names = TRUE,
                         cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  C = col_integer(),
  CT = col_integer(),
  ratio = col_double(),
  sites_with_data = col_integer(),
  chg_sites = col_integer()
))

# tile filter
# also here there is an additional CHG filter requiring there to be at least 6 CG sites
# CHECK IF THIS SHOULD ALSO REQUIRE CHG COVERAGE FILTER
mC_tile_table_filtered <-
mC_tile_table %>%
  unite(tile, chr, start,end, sep=":") %>%
  mutate(cov= CT/chg_sites) %>%
  inner_join(tile_list, "tile") %>%
  select(tile, ratio, chg_sites, cov) %>%
  filter(chg_sites >= 6)

mC_tile_table_filtered

write.table(mC_tile_table_filtered, paste0(output_folder, "/", sample, "_BSMAP_out.txt.100.CHG_filtered.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
