#!/usr/bin/Rscript
##########

#Peter Crisp
#2019-12-12
#R script to create 100bp tile reference

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
genome <- args[1]
genome

###########################
#setup
library(data.table)
library(tidyverse)
old.scipen <- getOption("scipen")
options(scipen=999)
###########################

# use data table to make 100 bp tile reference
# genome = "Sbicolor_454_v3.0.1"
chrom.sizes_file <- read_tsv(paste0("../", genome, ".chrom.sizes"), col_names = c("chr", "size"))

# # This method creates a whole tile at the end of the chromosomes, eventhough it will often be a partial tile. 
# # This tile get omitted in the next step by bedtools getfasta because it is beyond the end of the chromosome
#  reference_tiles_table <- data.table(chrom.sizes_file)[, seq(1, size, by=100), by=chr]

# This method omits the partial tile at the end of the chromosome
reference_tiles_table <- data.table(chrom.sizes_file)[, seq(1, size - (size %% 100), by=100), by=chr]

### make into proper bed (zero based)
reference_tiles <- as.tibble(reference_tiles_table)
reference_tiles_2 <- reference_tiles %>% rename(start = V1) %>% 
  mutate(end = start+99, start_zBased = start-1)

reference_tiles_2

# # A tibble: 7,097,206 x 4
#    chr   start   end start_zBased
#    <chr> <dbl> <dbl>        <dbl>
#  1 Chr11     1   100            0
#  2 Chr11   101   200          100
#  3 Chr11   201   300          200
#  4 Chr11   301   400          300
#  5 Chr11   401   500          400
#  6 Chr11   501   600          500
#  7 Chr11   601   700          600
#  8 Chr11   701   800          700
#  9 Chr11   801   900          800
# 10 Chr11   901  1000          900

# write output using write.table because write_tsv coerces to scientific notation
write.table(reference_tiles_2, paste0(genome, "_100bp_tiles_zBased.txt"),  sep = "\t", quote = F, row.names = F)

# write out as a bed file too
reference_tiles_bed <- reference_tiles_2 %>%
  mutate(start = start_zBased) %>%
  select(chr, start, end) %>%
  mutate(name = paste0(chr, ":", start, "-", end))

#awk '{ if (NR==1622117) print $0 }' GDDH13_1-1_formatted_Pt_100pb_tiles_for_sites_calc.bed

# again use write.table
write.table(reference_tiles_bed, paste0(genome, "_100bp_tiles.bed"),  sep = "\t", quote = F, row.names = F, col.names = F)


############ 100 bp tile reference for calling sites (+2nt) #############

# tiles <- read_tsv("../ZmaysPH207_443_v1_0_100bp_tiles_zBased.txt", col_names = F, 
#                   cols_only(
#                     X1 = col_character(),
#                     X2 = col_integer(),
#                     X3 = col_integer()
#                   ))

reference_tiles_2 %>% distinct(chr)

tiles_for_sites <- reference_tiles_2 %>% 
  # slice(1:10) %>%
  group_by(chr) %>%
  mutate(new_start = ifelse(start_zBased >= 2, start_zBased - 2, 0),
         new_end = ifelse(end < max(end), end + 2, max(end)))

tiles_for_sites %>% filter(start == min(start) | end == max(end))

tiles_for_sites_slim <- tiles_for_sites %>% select(chr, new_start, new_end)

tiles_for_sites_slim %>% distinct(chr)

write.table(tiles_for_sites_slim, paste0(genome, "_100pb_tiles_for_sites_calc.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
