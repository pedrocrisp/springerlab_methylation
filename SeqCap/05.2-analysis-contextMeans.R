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
data_path <- "analysis/BSMAPratio/"
data_path
out_folder <- "analysis_context"
out_folder

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
library(ggplot2)
###########################

# a sample
# sample = "US_1_Index1_S11"

# read file, but only required columns
# specify that chromosome is a character column
context_file <- read_tsv(paste0(data_path, "/", sample, "_BSMAP_out.txt"), col_names = c("chr", "start", "end", "strand", "context", "ratio", "eff_CT_count", "C_count", "CT_count", "rev_G_count", "rev_GA_count", "CI_lower", "CI_upper"),
                            cols_only(
                              context = col_character(),
                              C_count = col_double(),
                              CT_count = col_double()
                            ))
# summarise
context_file_summary <- context_file %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_context_summary.tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT_count >= 20) %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_context_summary_CT20.tsv"), sep = "\t", quote = F, row.names = F)

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
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary

write.table(context_file_summary, paste0(out_folder, "/", sample, "_subcontext_summary.tsv"), sep = "\t", quote = F, row.names = F)

# Do a second iteration with a coverage filter to compare
context_file_summary_CT20 <- context_file %>%
  filter(CT_count >= 20) %>%
  mutate(ratio = C_count/CT_count*100) %>%
  group_by(context) %>%
  summarise(mean = mean(ratio),
            median = median(ratio),
            sd=sd(ratio),
            n=n(),
            q5= quantile(ratio, .05),
            q95= quantile(ratio, .95)
            )

context_file_summary_CT20

write.table(context_file_summary_CT20, paste0(out_folder, "/", sample, "_subcontext_summary_CT20.tsv"), sep = "\t", quote = F, row.names = F)
