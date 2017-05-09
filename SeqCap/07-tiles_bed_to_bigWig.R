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
reference_tile_file <- args[3]
reference_tile_file

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

#reference used to make amendments
reference_tiles <- read_tsv(reference_tile_file, col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              start_zBased = col_integer(),
                              cg_sites = col_integer(),
                              chg_sites = col_integer(),
                              chh_sites = col_integer()
                            ))

#########
# fix bed CG
# note: using write.table, write_delim converts to scientific notation
# and cannot seem to disable that
broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.bed"), col_names = F,
                            cols(
                              X1 = col_character(),
                              X2 = col_integer(),
                              X3 = col_integer(),
                              X4 = col_number()
                            ))

colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
#remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph2 <-
fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, cg_sites) %>% arrange(chr, start)
write.table(fixed_bedGraph3, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
rm(fixed_bedGraph3)

#########
# fix bed CHG
# note: using write.table, write_delim converts to scientific notation
# and cannot seem to disable that
broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.bed"), col_names = F,
                            cols(
                              X1 = col_character(),
                              X2 = col_integer(),
                              X3 = col_integer(),
                              X4 = col_number()
                            ))
colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
#remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph2 <-
fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, chg_sites) %>% arrange(chr, start)
write.table(fixed_bedGraph3, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
rm(fixed_bedGraph3)

#########
# fix bed CHH
# note: using write.table, write_delim converts to scientific notation
# and cannot seem to disable that
broken_bedGraph <- read_tsv(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.bed"), col_names = F,
                            cols(
                              X1 = col_character(),
                              X2 = col_integer(),
                              X3 = col_integer(),
                              X4 = col_number()
                            ))
colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "sites_with_data", "C", "CT", "ratio")
#remove NAs - this is necessary because the output from the perl script is any tile with data in any context per sample
broken_bedGraph <- na.omit(broken_bedGraph)

#fix
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph2 <-
fixed_bedGraph %>% select(chr, start_zBased, end, ratio) %>% arrange(chr, start_zBased)
write.table(fixed_bedGraph2, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
rm(fixed_bedGraph2)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph3 <-
fixed_bedGraph %>% select(chr, start, end, C, CT, ratio, sites_with_data, chh_sites) %>% arrange(chr, start)
write.table(fixed_bedGraph3, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
rm(fixed_bedGraph3)
