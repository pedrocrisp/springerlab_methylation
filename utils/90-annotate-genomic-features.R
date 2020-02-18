#!/usr/bin/Rscript
##########

#Peter Crisp
#2020-20-7
#R script to call UMRs using a 100bp tile

# Notes
# Currently this script removes the organelles by using filter(!chr %in% c("Mt", "Pt")) - this may not catch all organelles 
# post filtering may be required to remove organelles or other undesired contigs

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
sample_prefix <- args[1]
inputFile <- args[2]
reference_annotation_mame <- args[3]
summary_output_folder <- args[4]
# coverage_filter_min <- as.double(args[5])
# site_filter_min <- as.double(args[6])
# MR_percent <- as.double(args[7])
# UMR_percent <- as.double(args[8])

######## de bug
# args
# sample_prefix <- "Sbicolor_SRR3286309_UMRs_6col"
# inputFile <- "Sbicolor_SRR3286309_UMRs_6col_olap_gene_v1.bed"
# reference_annotation_mame= "gene_v1"
# summary_output_folder="~/ws/analysis/ongoing/mC_sorghum/analysis/Sbicolor_SRR3286309_mC_domains_II_cov_3_sites_2_MR_0.4_UMR_0.1/mC_UMT_annotation"
# 

########
###########################
library(tidyverse)
library(ggthemes)
# library("seqinr")
old.scipen <- getOption("scipen")
options(scipen=999)
# library(wesanderson)
library(RColorBrewer)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

###########################
# Module #1
###########################

overlaps <- read_tsv(inputFile,
                     col_names = c("chr", "start", "end", "feature", "ID", "optional",
                                   "b_chr", "b_start", "b_end", "b_feature_type", "b_score", "b_strand", "b_ID","distance"), 
                     cols(chr = col_character(),
                          # Olap_file = col_character(), 
                          b_chr = col_character()))

overlaps 
# 34,514 (acrs 32,481) - only a few extra multi overlaps to resolve!

# resovle multi-overlaps
# I think most of the multi-overlaps will be at distance zero - this pipeline is mainly aimed at resolving that. Its possible some will be due to a peak being exactally the same length between an upstrea and downstream gene, but I might just pick one for these because I think they will be rare
# get distinct tile + feature rows
overlaps_distinct <- overlaps %>%
  select(chr, start, end, b_feature_type, distance) %>%
  distinct()
# 32,620 

overlaps %>%
  select(chr, start, end) %>%
  distinct()
# 32,481
# no multi overlaps because not overlapping features in this set and nother fell eqi-distant from 2 features I guess

# colapse, I take mean distance for tiles with multioverlaps. I think this is ok because gene-gene and TE-TE will called ambiguous because they are half way and for gene-TE if its distance is 0 then I will call it geneic TE and mean(0,0 = 0) and if distance is > 0 I'll just call it ambiguous. This scheme ensures no tile is double counted for distribution purposes.
overlaps_distinct_collapsed <- overlaps_distinct %>%
  group_by(chr, start, end) %>%
  arrange(chr, start, b_feature_type) %>%
  summarise(feature = paste(b_feature_type, collapse = "-"), distance2 = mean(distance))

overlaps_distinct_collapsed
# 32,481

overlaps_distinct_collapsed %>% group_by(feature) %>% summarise(n = n()) %>% print(n = 35)

# # A tibble: 14 x 2
#    feature                                      n
#    <chr>                                    <int>
#  1 .                                          135
#  2 lincRNA_gene                              1573
#  3 lincRNA_gene-nonSyntenic_gene                2
#  4 lincRNA_gene-syntenic_gene                  16
#  5 lincRNA_gene-tRNA_gene                       3
#  6 miRNA_gene                                 160
#  7 miRNA_gene-syntenic_gene                     2
#  8 nonSyntenic_gene                          5062
#  9 nonSyntenic_gene-syntenic_gene              39
# 10 nonSyntenic_gene-syntenic_gene-tRNA_gene     1
# 11 nonSyntenic_gene-tRNA_gene                  12
# 12 syntenic_gene                            24348
# 13 syntenic_gene-tRNA_gene                     63
# 14 tRNA_gene                                 1065

############ Gene overlaps 


# sRNA data annotation rules
# miRNA > syntenic_gene > non_syntenic > TE > ncRNA
# annotation does not need to have all these feature types, if it has none of these feature type then it will be called other_feature

# This distance does now take into account gene strand...

overlaps_distinct_collapsed_filtered <- overlaps_distinct_collapsed %>% ungroup() %>% 
  # slice(1:10000) %>%
  mutate(classification = ifelse(grepl("miRNA_gene", feature) & distance2 == 0, "miRNA_gene",
                                 ifelse(grepl("syntenic_gene", feature) & distance2 == 0, "syntenic_gene",
                                        ifelse(grepl("nonSyntenic_gene", feature) & distance2 == 0, "nonSyntenic_gene",
                                               ifelse(grepl("tRNA_gene", feature) & distance2 == 0, "tRNA_gene",
                                                      ifelse(grepl("lincRNA_gene", feature) & distance2 == 0, "lincRNA_gene",
                                                             ifelse(grepl("gene", feature) & distance2 == 0, "gene",
                                                                    ifelse(distance2 == 0, "other_feature",
                                                             ifelse(distance2 < 0 & distance2 >= -1000, "1kb_upstream",
                                                                    ifelse(distance2 < -1000 & distance2 >= -2000, "2kb_upstream",
                                                                           ifelse(distance2 > 0 & distance2 <= 1000, "1kb_downstream",
                                                                                  ifelse(distance2 > 1000 & distance2 <= 2000, "2kb_downstream", 
                                                                                          "intergenic"))))))))))))

overlaps_distinct_collapsed_filtered
# 32,481

#######
# write 100 bp annotation?
overlaps_distinct_collapsed_filtered_out <- overlaps_distinct_collapsed_filtered %>%
  mutate(score = ".",
         strand = ".") %>%
  select(chr, start, end, classification, score, strand)

write.table(overlaps_distinct_collapsed_filtered_out, paste0(sample_prefix, "_annotated_", reference_annotation_mame, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

########
#summarise
overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered %>%
  group_by(classification) %>%
  summarise(n = n())

overlaps_distinct_collapsed_filtered_summary <- overlaps_distinct_collapsed_filtered_summary %>% 
  mutate(percentage = n/sum(n)*100,
         Mb = n * 100 / 1000000)

overlaps_distinct_collapsed_filtered_summary
# # A tibble: 10 x 4
#    classification          n percentage     Mb
#    <chr>               <int>      <dbl>  <dbl>
#  1 1kb_downstream_gene  3008      9.26  0.301 
#  2 1kb_upstream_gene    4168     12.8   0.417 
#  3 2kb_downstream_gene  1120      3.45  0.112 
#  4 2kb_upstream_gene    1383      4.26  0.138 
#  5 intergenic           9139     28.1   0.914 
#  6 lincRNA_gene          520      1.60  0.052 
#  7 miRNA_gene             94      0.289 0.0094
#  8 nonSyntenic_gene     1207      3.72  0.121 
#  9 syntenic_gene       11369     35.0   1.14  
# 10 tRNA_gene             473      1.46  0.0473


write_tsv(overlaps_distinct_collapsed_filtered_summary, paste0(summary_output_folder, "/", sample_prefix, "_annotated_", reference_annotation_mame, "_summary.tsv"))

####### pull distal

distal_only_bed <- overlaps_distinct_collapsed_filtered_out %>% filter(classification == "intergenic") %>% 
  mutate(name = sample_prefix) %>% select(chr, start, end, name, score, strand)
# 8,972
# write

distal_only_bed %>% distinct(chr)

write.table(distal_only_bed, paste0(sample_prefix, "_annotated_", reference_annotation_mame, "_distal.bed"), sep = "\t", quote = F, row.names = F, col.names = F)



