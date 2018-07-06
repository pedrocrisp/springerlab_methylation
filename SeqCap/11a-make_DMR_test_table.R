#!/usr/bin/Rscript
##########

# Peter Crisp
# 2018-07-06
# R script to make a table of all contrasts for a list of samples
# Make redundant comparison eg A vs B and B vs A (twice as big as necessary)

# To call
# R -f ~/gitrepos/springerlab_methylation/SeqCap/11a-make_DMR_test_table.R --args <sample_list.txt>

library(tidyverse)
library(tools)

############ read in file stack
  # ARGS
#sample_list_to_analyse <-"samples_HC_and_GD_minus3.txt"

args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample_list_to_analyse <- args[1]
sample_list_to_analyse

sample_list_subset <- read_tsv(sample_list_to_analyse, col_names = c("sample"))

# sample combinations
DE_test_combos <- as.tibble(t(combn(pull(sample_list_subset, sample), 2)))
colnames(DE_test_combos) <- c("sample1", "sample2")

DE_test_combos_inverse <- DE_test_combos
colnames(DE_test_combos_inverse) <- c("sample2", "sample1")

DE_test_combos_all <- bind_rows(DE_test_combos, DE_test_combos_inverse)

# add filenames and test name
DE_test_combos_all <- DE_test_combos_all %>%
  mutate(test_name = paste0(sample1, ".vs.", sample2)) %>%
  select(test_name, sample1, sample2)

DE_test_combos_all

write_tsv(DE_test_combos_all, "DMR_tests_combos_all_table.tsv")

DE_test_combos_all_names <- select(DE_test_combos_all, test_name)

write_tsv(DE_test_combos_all_names, "DMR_tests_combos_all_list.txt", col_names = F)
