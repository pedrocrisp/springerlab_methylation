#!/usr/bin/Rscript
##########

# Peter Crisp
# 2018-07-09
# R script to call DMRs using DSS, supplying 100 bp tiles

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
contrast <- args[1]
contrast
DMR_contrasts_table_file <- args[2]
DMR_contrasts_table_file
path_to_data_files <- args[3]
path_to_data_files

####### libs
library(tidyverse)
library(DSS)
library(bsseq)
library(limma)
library(tidygenomics)
library(ggthemes)
old.scipen <- getOption("scipen")
options(scipen=999)

#######  args

# contrast = "PAC004_root_HC_mother.vs.PAC012_gen_ST_0"
# DMR_contrasts_table_file = "DMR_tests_combos_all_table.tsv"
# path_to_data_files = "analysis/tiles_filtered_4C_2x"

####### set up

outFolder <- paste0(path_to_data_files, "_DSS_DMRs")
dir.create(outFolder)

DMR_contrasts_table <- read_tsv(DMR_contrasts_table_file)

sample1 = pull(filter(DMR_contrasts_table, test_name == contrast), sample1)
sample1

sample2 = pull(filter(DMR_contrasts_table, test_name == contrast), sample2)
sample2

print(paste0("Calling DMR tiles for ", sample1, " vs ", sample2))

#######  #######
# CG
#######  #######

context = "CG"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

####### 2 Perform statistical test for DML by calling DMLtest function
# t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T)
# proc.time() - t1

# user  system elapsed
# 508.827   0.139 509.181

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# merge adjacent DMRs
# not the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)
dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)

# dmls_bed %>% select(chr, start, end, diff) %>% arrange(chr, start)
# dmls_bed %>% select(chr, start, end, diff, cluster_id) %>% arrange(cluster_id)

dmls_bed_merged <- dmls_bed %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

# no directionality to DMR 59,731; if hypo/hyper is considered 60,134

dmls_bed_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_summarised

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest))
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut of ther columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain contect sprcific tile lists for comparison

CG_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CG_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CG_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  #######
# CHG
#######  #######

context = "CHG"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

####### 2 Perform statistical test for DML by calling DMLtest function
t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T)
proc.time() - t1

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# merge adjacent DMRs
# not the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)
dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)

dmls_bed_merged <- dmls_bed %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

dmls_bed_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_summarised

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest))
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut of ther columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain contect sprcific tile lists for comparison

CHG_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CHG_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CHG_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  #######
# CHH
#######  #######

context = "CHH"

#######  read in the files
sample1_data <- read_tsv(file.path(path_to_data_files, paste0(sample1, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample1_data <- sample1_data %>% select("chr", "pos", "N", "X", "sites")
# sample1_data <- sample1_data %>% filter(chr == "Chr01")
sample1_data

sample2_data <- read_tsv(file.path(path_to_data_files, paste0(sample2, "_BSMAP_out.txt.100.", context, "_filtered.txt")), col_names = c("chr", "pos", "X", "N", "sites"), skip = 1)
sample2_data <- sample2_data %>% select("chr", "pos", "N", "X", "sites")
# sample2_data <- sample2_data %>% filter(chr == "Chr01")
sample2_data

#######  1. create an object of BSseq class (requires bsseq Bioconductor package)

BSobj <- makeBSseqData(list(sample1_data, sample2_data), c(sample1, sample2) )
BSobj

####### 2 Perform statistical test for DML by calling DMLtest function
t1 <- proc.time()
dmlTest <- DMLtest(BSobj, group1=sample1, group2=sample2, smoothing = F, equal.disp = T)
proc.time() - t1

####### 3. Call DMRs

dmls <- callDML(dmlTest, p.threshold=0.01, delta = 0.1)
dmls_tbl <- as.tibble(dmls)
dmls_tbl

# merge adjacent DMRs
# not the bed must be ordered - tidy genomics does not check!
dmls_bed <- dmls_tbl %>% mutate(chr = as.character(chr), start = as.double(pos), end = pos+100, type = ifelse(diff > 0, "hyper", "hypo")) %>%
  select(chr, start, end, diff, type)  %>%
  arrange(chr, start)
dmls_bed_hyper <- filter(dmls_bed, type == "hyper")
dmls_bed_hypo <- filter(dmls_bed, type == "hypo")

dmls_bed_hyper <- genome_cluster(dmls_bed_hyper, by=c("chr", "start", "end"), max_distance = 1)
dmls_bed_hypo <- genome_cluster(dmls_bed_hypo, by=c("chr", "start", "end"), max_distance = 1)

dmls_bed <- bind_rows(dmls_bed_hyper, dmls_bed_hypo)

dmls_bed_merged <- dmls_bed %>% group_by(cluster_id, type) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), mean_diff = mean(diff), max_diff = max(diff), min_diff = min(diff)) %>%
  mutate(length = end - start) %>%
  ungroup()
dmls_bed_merged

dmls_bed_summarised <- dmls_bed_merged %>%
  group_by(type) %>%
  summarise(n = n(), mean_length = mean(length), max_length = max(length), mean_diff = mean(mean_diff))
dmls_bed_summarised

####### 4. Make results files and write out tables

dmlTest_calls <- as.tibble(dmlTest) %>% left_join(dmls_tbl, by = colnames(dmlTest)) 
dmlTest_calls

### write DMR calls (large) file - only retain mu1 and mu2 cut of ther columns to make file smaller
dmlTest_calls_out <- dmlTest_calls %>% select(chr, pos, mu1, mu2)
write_tsv(dmlTest_calls_out, paste0(outFolder, "/", contrast, "_", context, "_dmlTest_calls.tsv"))

### write only DMR tiles
write_csv(dmls_tbl, paste0(outFolder, "/", contrast, "_", context, "_DMRs.csv"))

### write summary of merged DMR tiles
write_csv(dmls_bed_summarised, paste0(outFolder, "/", contrast, "_", context, "_DMR_merged_summary.csv"))

####### 5. Retain contect sprcific tile lists for comparison

CHH_tiles_tested <- dmlTest_calls %>% mutate(tile = paste0(chr, "_", pos)) %>% select(tile)

CHH_DMRs <- dmls_tbl %>% mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, 1, -1)) %>% select(tile, DMR)

CHH_mu <- dmls_tbl %>%
  mutate(tile = paste0(chr, "_", pos), DMR = ifelse(diff > 0, "hyper", "hypo")) %>%
  mutate(DMR_type = paste0(context, ".", DMR)) %>%
  select(tile, DMR_type, mu1, mu2)

#######  ####### #######  ####### #######  ####### #######  ####### #######  ####### #######  #######
# Some analysis
#######  ####### #######  ####### #######  ####### #######  ####### #######  ####### #######  #######

outFolder_analysis <- paste0(outFolder, "/Analysis_", contrast)
dir.create(outFolder_analysis)

# Overlap between tiles that passed filtering

CG_tiles_tested <- CG_tiles_tested %>% mutate(context = "CG", pass = 1)
CHG_tiles_tested <- CHG_tiles_tested %>% mutate(context = "CHG", pass = 1)
CHH_tiles_tested <- CHH_tiles_tested %>% mutate(context = "CHH", pass = 1)

tile_lists <- bind_rows(CG_tiles_tested, CHG_tiles_tested, CHH_tiles_tested)

table_venn <- tile_lists %>%
  mutate(pass = as.integer(pass)) %>%
  spread(key = context, value = pass, fill = 0)

table_venn_summary <- table_venn %>%
  group_by(CG, CHG, CHH) %>%
  summarise(Counts = n()) %>%
  ungroup() %>%
  mutate(tiles_pct = round(Counts/sum(Counts)*100, digits = 1))

table_venn_summary

ref_venn <- as.tibble(expand.grid(0:1, 0:1, 0:1)) %>% mutate(Counts = 0, tiles_pct = 0)
colnames(ref_venn) <- colnames(table_venn_summary)

ref_venn

table_venn_summary_all <- ref_venn %>% anti_join(table_venn_summary, by = c(colnames(table_venn_summary)[1:3])) %>% bind_rows(table_venn_summary)

table_venn_summary_counts <- table_venn_summary_all %>%
  select(-tiles_pct)

venn_data <- as.matrix(table_venn_summary_counts)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_tiles_assessed_context_overlap.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "Tiles assessed per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

table_venn_summary_2 <- table_venn_summary_all %>%
  mutate(Counts = tiles_pct) %>%
  select(-tiles_pct)

venn_data <- as.matrix(table_venn_summary_2)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_tiles_assessed_context_overlap_pct.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "Tiles assessed per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

# Overlap between DMRs

CG_DMRs <- CG_DMRs %>% mutate(context = "CG")
CHG_DMRs <- CHG_DMRs %>% mutate(context = "CHG")
CHH_DMRs <- CHH_DMRs %>% mutate(context = "CHH")

tile_lists <- bind_rows(CG_DMRs, CHG_DMRs, CHH_DMRs)
tile_lists

table_venn <- tile_lists %>%
  mutate(DMR = as.integer(DMR)) %>%
  spread(key = context, value = DMR, fill = 0)

venn_data <- vennCounts(table_venn[,-1], include=c("both"))
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_both.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "all DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hyper DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()


venn_data <- vennCounts(table_venn[,-1], include=c("up"))
venn_data
pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_up.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "hyper DMRs per context",
            include=c("up"),
            counts.col=c('red'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hyper DMRs per context",
            include=c("up"),
            counts.col=c('red'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

venn_data <- vennCounts(table_venn[,-1], include=c("down"))
venn_data
pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_down.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "hypo DMRs per context",
            include=c("down"),
            counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)

vennDiagram(venn_data, names = "",
            # main= "hypo DMRs per context",
            include=c("down"),
            counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

venn_data <- vennCounts(table_venn[,-1], include=c("both"))
venn_data

venn_data_2 <- venn_data
class(venn_data_2) <- "matrix"

venn_data_pct <- as.tibble(venn_data_2) %>%
  mutate(Counts = round(Counts/sum(Counts)*100, digits = 1))

venn_data_pct

venn_data <- as.matrix(venn_data_pct)
class(venn_data) <- "VennCounts"
venn_data

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_context_overlap_both_pct.pdf"), width = 4, height = 4)
vennDiagram(venn_data,
            main= "DMRs per context",
            # include=c("down"),
            # counts.col=c('blue'),
            show.include=T,
            cex=c(0.75)
            # names=c("Excess-light 1 hour","Drought 7 days")
)
dev.off()

###### Density plot of mC levels for DMRs

DMR_lists <- bind_rows(CG_mu, CHG_mu, CHH_mu) %>% gather(key = sample, value = mC, -DMR_type, -tile)
DMR_lists

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_mC_values_distribution.pdf"), width = 6, height = 4)
ggplot(DMR_lists, aes(mC, colour = DMR_type)) +
  geom_density() +
  theme_minimal() +
  labs(title = "DMRs mC values distribution") +
  scale_color_ptol()
dev.off()

DMR_lists_2 <- DMR_lists %>% separate(DMR_type, into = c("context", "direction"))

pdf(paste0(outFolder_analysis, "/", contrast,"_DMRs_mC_values_distribution_facet.pdf"), width = 10, height = 7)
ggplot(DMR_lists_2, aes(mC, colour = sample)) +
  geom_density() +
  theme_minimal() +
  labs(title = "DMRs mC values distribution") +
  scale_color_ptol() +
  facet_grid(direction ~ context)
dev.off()
