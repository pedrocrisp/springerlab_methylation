#!/usr/bin/Rscript
##########

#Peter Crisp
#2017-04-11
#R script to make aggregate table from the individual summarise methylation output

####################################
#libraries
library(tidyverse)

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned

# list the folder with the individual sample output files output from the make ratios script, should end "_ratio_annotated.csv"
data_path <- args[1]
data_path
outPrefix <- args[2]
outPrefix

####################################
# read in dataset

# merge individual output files into a single final output table
# inspired by:
# http://stackoverflow.com/questions/39795136/how-can-i-manipulate-dataframe-columns-with-different-values-from-an-external-ve

# list the files
files <- dir(data_path, pattern = "*_ratio_annotated.csv") # get file names

# inspired by:
# #http://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R
# Read in as a list of dfs
data <- data_frame(filename = files) %>% # create a data frame
                                         # holding the file names
  mutate(file_contents = map(filename,          # read files into
           ~ read_csv(file.path(data_path, .))) # a new data column
        )

####################################
# tidy

# unnest the dfs
# works like rbind
#long-ish table, gather properly, then fully cast
big_data <- unnest(data)

# gather the data so that it can be spread with "sample"_"mC-type" as columns
big_data_spread <-
  big_data %>%
    separate(filename, c("sample", "suffix"), sep="_ratio_annotated") %>%
      select(-suffix)  %>%
        gather(variable, value, CG_CT:reads) %>%
          unite(temp, sample, variable) %>%
            spread(temp, value)

####################################
#write the final output table
write_csv(big_data_spread, paste0(outPrefix, "_mC_all_samples.csv"))

####################################
