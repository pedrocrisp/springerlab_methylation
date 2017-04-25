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
###########################

#reference used to make amendments
reference_tiles <- read_delim(reference_tile_file, delim ="\t")

#########
#fix bed CG
# note: using write.table, write_delim converts to scientific notation
# and cannot seem to disable that

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.bg"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start", "broken_end", "ratio")

fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start"))  %>% select(chr, start, end, ratio)

write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CG.fixed.bg"), sep = "\t", quote = F, row.names = F, col.names = F)

#########
#fix bed CHG

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.bg"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start", "broken_end", "ratio")

fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start"))  %>% select(chr, start, end, ratio)

write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHG.fixed.bg"), sep = "\t", quote = F, row.names = F, col.names = F)

#########
#fix bed CHH

broken_bedGraph <- read_delim(paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.bg"), delim ="\t", col_names = F)
colnames(broken_bedGraph) <- c("chr", "start", "broken_end", "ratio")

fixed_bedGraph <-
merge(broken_bedGraph, reference_tiles, by = c("chr", "start"))  %>% select(chr, start, end, ratio)

write.table(fixed_bedGraph, paste0(data_folder, "/", sample, "_BSMAP_out.txt.100.CHH.fixed.bg"), sep = "\t", quote = F, row.names = F, col.names = F)
