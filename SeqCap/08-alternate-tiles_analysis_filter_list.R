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
minCHHs <- args[5]
minCHHs
minCHH_cov <- args[6]
minCHH_cov
minCHGs <- args[7]
minCHGs
minCHG_cov <- args[8]
minCHG_cov
minCGs <- args[9]
minCGs
minCG_cov <- args[10]
minCG_cov

#debugging
#sample = "PAC004_root_HC_mother"
#data_folder = "analysis/tiles"
#tile_list = FALSE
#output_folder = "analysis/tiles_filtered"
#minCHHs <- 4
#minCHH_cov <- 2
#minCHGs <- 4
#minCHG_cov <- 2
#minCGs <- 4
#minCG_cov <- 2


###########################
#setup
library(tidyverse)
library(purrr)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

#########
#list of tiles to use for filtering
# read_tsv(tile_list)
if(!tile_list == FALSE){
tile_list <- read_tsv(tile_list)
return(tile_list)
}
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
  unite(tile, chr, start,end, sep=":", remove = F) %>%
  mutate(cov= CT/chh_sites) %>%
  when(tile_list == FALSE ~ ., ~ inner_join(tile_list, "tile")) %>%
  filter(chh_sites >= as.double(minCHHs) & cov >= as.double(minCHH_cov)) %>%
  select(chr, start, C, CT, chh_sites)

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
  unite(tile, chr, start,end, sep=":", remove = F) %>%
  mutate(cov= CT/cg_sites) %>%
  when(tile_list == FALSE ~ ., ~ inner_join(tile_list, "tile")) %>%
  filter(cg_sites >= as.double(minCGs) & cov >= as.double(minCG_cov)) %>%
  select(chr, start, C, CT, cg_sites)

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
  unite(tile, chr, start,end, sep=":", remove = F) %>%
  mutate(cov= CT/chg_sites) %>%
  when(tile_list == FALSE ~ ., ~ inner_join(tile_list, "tile")) %>%
  filter(chg_sites >= as.double(minCHGs) & cov >= as.double(minCHG_cov)) %>%
  select(chr, start, C, CT, chg_sites)

mC_tile_table_filtered

write.table(mC_tile_table_filtered, paste0(output_folder, "/", sample, "_BSMAP_out.txt.100.CHG_filtered.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
