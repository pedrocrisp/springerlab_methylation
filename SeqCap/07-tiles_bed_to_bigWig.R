#!/usr/bin/Rscript
##########

#Peter Crisp
#2017-04-08
#R script to 100bp tile data with incorrect chr end tiles
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
reference_tiles <- read_delim(reference_tile_file, delim ="\t")

#########
# fix bed CG
# note: using write.table, write_delim converts to scientific notation
# and cannot seem to disable that

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.bed"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "C", "CT", "ratio")
broken_bedGraph <- na.omit(broken_bedGraph)

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start, end, ratio) %>% arrange(chr, start)
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start_zBased, end, C, CT, ratio) %>% arrange(chr, start_zBased)
colnames(fixed_bedGraph) <- c("chr", "start", "end", "C", "CT", "ratio")
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

#########
#fix bed CHG

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.bed"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "C", "CT", "ratio")
broken_bedGraph <- na.omit(broken_bedGraph)

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start, end, ratio) %>% arrange(chr, start)
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start_zBased, end, C, CT, ratio) %>% arrange(chr, start_zBased)
colnames(fixed_bedGraph) <- c("chr", "start", "end", "C", "CT", "ratio")
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

#########
#fix bed CHH

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.bed"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start_zBased", "broken_end", "C", "CT", "ratio")
broken_bedGraph <- na.omit(broken_bedGraph)

#keep as BED format (zero-based coordinates)
# only retain chr, start, end, ratio; then sort to make it a bedGraph file
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start, end, ratio) %>% arrange(chr, start)
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.sorted.bg"), sep = "\t", quote = F, row.names = F, col.names = F)

#convert to one-based coordinate .txt file and sort
fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start_zBased"))  %>% select(chr, start_zBased, end, C, CT, ratio) %>% arrange(chr, start_zBased)
colnames(fixed_bedGraph) <- c("chr", "start", "end", "C", "CT", "ratio")
write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
