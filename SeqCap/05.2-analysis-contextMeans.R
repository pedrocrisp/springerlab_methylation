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

# loci_of_interst_file
loi_file <- read_tsv(loci_of_interst_file)


# a sample
# sample = "US_1_Index1_S11"

# read file, but only required columns
# specify that chromosome is a character column
context_file <- read_tsv(paste0(data_path, "/", sample, "_BSMAP_out.txt"), col_names = c("chr", "start", "end", "strand", "context", "ratio", "eff_CT_count", "C_count", "CT_count", "rev_G_count", "rev_GA_count", "CI_lower", "CI_upper"),
                            cols_only(
                              chr = col_double(),
                              start = col_double(),
                              context = col_character(),
                              C_count = col_double(),
                              CT_count = col_double()
                            ))
# # coverageFilter
# coverage_summary <- context_file %>%
#   summarise(mean = mean(CT_count),
#             median = median(CT_count),
#             mode = getMode(CT_count),
#             sd=sd(CT_count),
#             n=n(),
#             q5= quantile(CT_count, .05),
#             q95= quantile(CT_count, .95)
#             )
#
# write.table(coverage_summary, paste0(out_folder, "/", sample, "_coverage_summary.tsv"), sep = "\t", quote = F, row.names = F)

# summarise
context_file_summary <- context_file %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_context_summary.tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT_count >= coverageFilter) %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_context_summary_CT20.tsv"), sep = "\t", quote = F, row.names = F)

# Do a third iteration filtering based on a list of loci of interest
# Filters based on chr and start position

context_file_summary_CT20_loci_of_interest <- context_file %>%
  filter(CT_count >= coverageFilter) %>%
  inner_join(loi_file, by = c("chr", "start")) %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20_loci_of_interest

write.table(context_file_summary_CT20_loci_of_interest, paste0(out_folder, "/", sample, "_context_summary_CT20_", loci_of_interest_name,".tsv"), sep = "\t", quote = F, row.names = F)


########## subcontext

context_file <- read_tsv(paste0(data_path, "/", sample, "_BSMAP_out_subcontext.txt"), col_names = c("chr", "start", "end", "strand", "context", "ratio", "eff_CT_count", "C_count", "CT_count", "rev_G_count", "rev_GA_count", "CI_lower", "CI_upper"),
                            cols_only(
                              context = col_character(),
                              C_count = col_double(),
                              CT_count = col_double()
                            ))

context_file_summary <- context_file %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_subcontext_summary.tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT_count >= coverageFilter) %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            mode = getMode(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_subcontext_summary_CT20.tsv"), sep = "\t", quote = F, row.names = F)
