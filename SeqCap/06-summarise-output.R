#!/usr/bin/Rscript
##########

#Peter Crisp
#2017-04-08
#R script to summarise methylation output

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample <- args[1]
sample
data_folder <- args[2]
data_folder
outDir <- args[3]
outDir


###########################
#setup
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(limma)
###########################
# Process mC data}

#read in mC data
mC_data <- read.table(paste0(data_folder, "/", sample, "_BSMAP_out_ontarget_mC.txt"), sep="\t")

#give the columns names
colnames(mC_data) <- c("chr", "start", "stop", "context", "sites", "mC", "CT")

#dplyr oh yeah - melt, cast, sort, add columns
mC_data_processed <-
  mC_data %>%
    mutate(specific_target = as.character(paste0(chr, ":", start, ":",stop)), ratio = format(round(mC/CT,2), nsmall = 2)) %>%
      arrange(specific_target) %>%
        select(context:ratio) %>%
          melt(variable.name = "key", value.names = "value", id.vars = c("specific_target", "context")) %>%
            mutate(context_value=paste0(context, "_", key)) %>%
              dcast(specific_target~context_value, value.var = "value") %>%
                separate(specific_target, into = c("v4_Chr", "v4_start", "v4_end"), sep= '[:-]') %>%
                  mutate(v4spec = paste0(v4_Chr, ":", v4_start, "-", v4_end))

###########################
# process read count data}

#read in data: read counts
read_counts_data <- read.table(paste0(data_folder, "/", sample, "_specific_region_count_pairs_clipOverlap.txt"), sep="\t", stringsAsFactors=FALSE)
colnames(read_counts_data) <- c("v4spec", "reads")

print("how may regions have read data")
length(read_counts_data$v4spec)

print("how many regions have mC data")
length(mC_data_processed$v4spec)

print("check read counts match mC data")
read_counts_data[which(!read_counts_data$v4spec %in% mC_data_processed$v4spec),]
length(read_counts_data[which(!read_counts_data$v4spec %in% mC_data_processed$v4spec),1])


print("check the reverse matches?")
read_counts_data[which(!mC_data_processed$v4spec %in% read_counts_data$v4spec),]
length(read_counts_data[which(!mC_data_processed$v4spec %in% read_counts_data$v4spec),1])

#merge

final_data <- merge(mC_data_processed, read_counts_data, by="v4spec", all = TRUE)

###########################

# annotation
SeqCap_v2_annotation_metadata_v4_final <- read.csv("~/umn/refseqs/maize/SeqCap/Seqcap_ultimate_annotation_files/SeqCapEpi2_regions_annotation_v2_v4.csv")

final_data_annotated <- merge(SeqCap_v2_annotation_metadata_v4_final, final_data, by.x = "v4specfic",by.y ="v4spec", all = TRUE, suffixes = c("_filtered_data", "_ref_annotation"))

###########################
#write output

write.csv(final_data_annotated, paste0(outDir, "/", sample, "_ratio_annotated.csv"), row.names = FALSE)
