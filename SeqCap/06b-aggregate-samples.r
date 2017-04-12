#!/usr/bin/Rscript
##########

#Peter Crisp
#2017-04-11
#R script to make aggregate table from the individual summarise methylation output

####################################
#Usage
#Rscript \
#~/gitrepos/springerlab_methylation/SeqCap/06b-aggregate-samples.R \
#OnTargetCoverage_annotated \
#Mei_final \
#~/umn/refseqs/maize/SeqCap/Seqcap_ultimate_annotation_files/SeqCapEpi2_regions_annotation_v2_v4.csv

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
annotationFile <- args[3]
annotationFile
keyfile <- args[4]

####################################
# read in dataset

# merge individual output files into a single final output table
# inspired by:
# http://stackoverflow.com/questions/39795136/how-can-i-manipulate-dataframe-columns-with-different-values-from-an-external-ve

# list the files
files <- dir(data_path, pattern = "*_ratio_annotated.csv") # get file names
files

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

# rename samples
key <- read_csv(keyfile)
big_data <- merge(key, big_data, by ="filename", all = TRUE)

# gather the data so that it can be spread with "sample"_"mC-type" as columns
#split off all metadata except SampleID and unique region variable - doesnt work otherwise, plus this is quicker
big_data_spread <-
  big_data %>%
    #separate(filename, c("sample", "suffix"), sep="_ratio_annotated") %>%
      select(sample_description, Region, CG_CT:reads)  %>%
        gather(variable, value, -(sample_description:Region)) %>%
        # gather(variable, value, CG_CT:reads) %>%
          unite(temp, sample_description, variable) %>%
            spread(temp, value)

####################################
## re annotate with metadata
#annotationFile = "~/umn/refseqs/maize/SeqCap/Seqcap_ultimate_annotation_files/SeqCapEpi2_regions_annotation_v2_v4.csv"
#get annotation
SeqCap_v2_annotation_metadata_v4_final <- read_csv(annotationFile)

#merge and sort by coordinate
final_data_annotated <- merge(SeqCap_v2_annotation_metadata_v4_final, big_data_spread, by ="Region", all = TRUE) %>% arrange(v4_Chr,	v4_start,	v4_end)

####################################
#write the final output table
write_csv(final_data_annotated, paste0(outPrefix, "_mC_all_samples.csv"))

####################################
