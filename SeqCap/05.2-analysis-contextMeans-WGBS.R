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
data_path <- args[2]
data_path
coverageFilter <- args[3]
coverageFilter
out_folder <- "analysis/contextMeans"
out_folder
loci_of_interst_file <- args[4]
loci_of_interst_file
loci_of_interest_name <- args[5]
loci_of_interest_name

#sample <- "BR3_3"
#data_path <- "analysis/BSMAPratio"
#coverageFilter <- 1
#out_folder <- "analysis/contextMeans"
#loci_of_interst_file <- "none"
#loci_of_interest_name <- "none"

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
library(ggplot2)
###########################

# mode funciton
getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

if (!loci_of_interst_file == "none") {
# loci_of_interst_file
loi_file <- read_tsv(loci_of_interst_file)
}


############## Context
############## CG, CHG. CHH
# a sample
# sample = "BR3_3"

context_means <- function(context){
# context = "CHH"

# read file, but only required columns
# specify that chromosome is a character column
context_file <- read_tsv(paste0(data_path, "/", sample, "_BSMAP_out_", context, ".bedGraph"),
                      col_names = c("chr", "start", "end", "ratio", "CT"),
                      cols_only(ratio = col_double(), CT = col_double())
                      )

# summarise
context_file_summary <- context_file %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_context_summary_", context, ".tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT >= coverageFilter) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_context_summary_", context, "_CT20.tsv"), sep = "\t", quote = F, row.names = F)

# Do a third iteration filtering based on a list of loci of interest
# Filters based on chr and start position

if (!loci_of_interst_file == "none") {

context_file_summary_CT20_loci_of_interest <- context_file %>%
  filter(CT_count >= coverageFilter) %>%
  inner_join(loi_file, by = c("chr", "start")) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20_loci_of_interest

write.table(context_file_summary_CT20_loci_of_interest, paste0(out_folder, "/", sample, "_context_summary_", context, "_CT20_", loci_of_interest_name,".tsv"), sep = "\t", quote = F, row.names = F)

}

rm(context_file)

}

### run the function for all contexts
walk(c("CG", "CHG", "CHH"), context_means)

#####################
########## subcontext

sub_context_means <- function(context){
# context = "CHH"

# read file, but only required columns
# specify that chromosome is a character column
context_file <- read_tsv(paste0(data_path, "/", sample, "_BSMAP_out_subcontext_", context, ".bedGraph"),
                      col_names = c("chr", "start", "end", "ratio", "CT"),
                      cols_only(ratio = col_double(), CT = col_double())
                      )

# summarise
context_file_summary <- context_file %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_subcontext_summary_", context, ".tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT >= coverageFilter) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_subcontext_summary_", context, "_CT20.tsv"), sep = "\t", quote = F, row.names = F)

# Do a third iteration filtering based on a list of loci of interest
# Filters based on chr and start position

if (!loci_of_interst_file == "none") {

context_file_summary_CT20_loci_of_interest <- context_file %>%
  filter(CT_count >= coverageFilter) %>%
  inner_join(loi_file, by = c("chr", "start")) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20_loci_of_interest

write.table(context_file_summary_CT20_loci_of_interest, paste0(out_folder, "/", sample, "_subcontext_summary_", context, "_CT20_", loci_of_interest_name,".tsv"), sep = "\t", quote = F, row.names = F)

}

rm(context_file)

}

### run the function for all contexts
walk(c("CG",
      "CAG",
      "CCG",
      "CTG",
      "CAA",
      "CAC",
      "CAT",
      "CCA",
      "CCC",
      "CCT",
      "CTA",
      "CTC",
      "CTT"), sub_context_means)
