#!/usr/bin/Rscript
##########

#Peter Crisp
#2019-12-12
#R script to count number of sites per 100bp tile

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
genome <- args[1]
genome

###########################
#setup
library(tidyverse)
# disable scientific notation
old.scipen <- getOption("scipen")
options(scipen=999)
library(ggplot2)
library(ggthemes)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

# text sizes
text_size_theme_8_labels <- theme(axis.text=element_text(size=8),
                                  axis.title=element_text(size=8),
                                  axis.text.x=element_text(angle = 45, hjust = 1),
                                  legend.title=element_text(size=8),
                                  legend.text=element_text(size=8),
                                  panel.background = element_rect(fill = "transparent") # bg of the panel
                                  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
                                  , panel.grid.major = element_blank() # get rid of major grid
                                  , panel.grid.minor = element_blank() # get rid of minor grid
                                  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
                                  , legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg
# magic geom_text conversion ratio
# https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
label_size = 25.4/72 * 8

###########################

#########  3. Merge in C sites data ##########

# genome = "Sbicolor_454_v3.0.1"
reference_tiles <- read_tsv(paste0(genome, "_100bp_tiles_zBased.txt"), col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              start_zBased = col_integer()
                            ))

unique(reference_tiles$chr)

reference_tiles_2 <- reference_tiles
# %>% filter(chr %in% c(1:10, "Pt", "Mt"))

reference_tiles_2
# A tibble: 7,097,206 x 4
#    chr   start   end start_zBased
#    <chr> <int> <int>        <int>
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

# unique(reference_tiles_2$chr)

##### read in site files and parse

# read in CG/CHG/CHH sites reference file created above
reference_tiles_C_sites <- read_tsv(paste0(genome, "_100bp_tiles_sites.tsv"), col_names = TRUE,
                                    cols(
                                      seqID = col_character(),
                                      patternName = col_character(),
                                      frequency = col_double()
                                    ))

# parse to make a bed file
reference_tiles_C_sites_bed <- reference_tiles_C_sites %>%
  #head(n=1000) %>%
  separate(seqID, into = c("chr", "start", "end"), sep = c(":|-")) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  mutate(start = ifelse(start >= 1, start+2, 0), end = end-2) %>%
  spread(key = patternName, value = frequency) %>%
  arrange(chr, start) %>%
  rename(cg_sites = CG,
         chg_sites = CHG,
         chh_sites = CHH
  )

reference_tiles_C_sites_bed

# # A tibble: 6,253,575 x 6
#    chr   start   end cg_sites chg_sites chh_sites
#    <chr> <dbl> <dbl>    <dbl>     <dbl>     <dbl>
#  1 Chr00     0   100       12         8        34
#  2 Chr00   100   200       NA        NA        18
#  3 Chr00   200   300        4         3         8
#  4 Chr00   300   400       10         2        26
#  5 Chr00   400   500       10         4        26
#  6 Chr00   500   600       12         2        27
#  7 Chr00   600   700        8         2        25
#  8 Chr00   700   800        4         4        23
#  9 Chr00   800   900       NA         4        26
# 10 Chr00   900  1000       NA         4        31

# unique(reference_tiles_C_sites_bed$chr)

#merge
reference_tiles_sites <- left_join(reference_tiles_2, reference_tiles_C_sites_bed, by = c("chr", "start_zBased" ="start", "end"))
reference_tiles_sites

# # A tibble: 7,097,206 x 7
#    chr   start   end start_zBased cg_sites chg_sites chh_sites
#    <chr> <int> <dbl>        <dbl>    <dbl>     <dbl>     <dbl>
#  1 Chr11     1   100            0        2         6        30
#  2 Chr11   101   200          100        6         3        29
#  3 Chr11   201   300          200        2         3        32
#  4 Chr11   301   400          300        4         3        33
#  5 Chr11   401   500          400        2         6        27
#  6 Chr11   501   600          500       NA         8        31
#  7 Chr11   601   700          600        2         5        21
#  8 Chr11   701   800          700        2         8        26
#  9 Chr11   801   900          800        2         6        25
# 10 Chr11   901  1000          900        2        11        37

# replace NA with 0 and re order
# which(is.na(reference_tiles_sites))

reference_tiles_sites <- reference_tiles_sites %>%
  replace(., is.na(.), 0) %>%
  mutate(chr = factor(chr, levels = unique(chr))) %>%
  arrange(chr, start)

unique(reference_tiles_sites$chr)

#write the one-based txt file
write.table(reference_tiles_sites, paste0(genome, "_100bp_tiles_zBased_sites_counts.txt"), sep = "\t", quote = F, row.names = F)

##############

#write a bed file (zero based and no head)
reference_tiles_2_bed <- select(reference_tiles_sites, chr, start_zBased, end, cg_sites, chg_sites, chh_sites)

write.table(reference_tiles_2_bed, paste0(genome, "_100bp_tiles_zBased_sites_counts.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

##############
##############
##############
## 4. Ref inc Ns - Reformat and merge in with whole tile file

#########  4. Merge in C sites data ##########
# genome = "Bdistachyon_314_v3.0"

reference_tiles <- read_tsv(paste0(genome, "_100bp_tiles_zBased.txt"), col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              start_zBased = col_integer()
                            ))

# unique(reference_tiles$chr)

reference_tiles_2 <- reference_tiles
# %>% filter(chr %in% c(1:10, "Pt", "Mt"))

reference_tiles_2
# A tibble: 7,097,206 x 4
#    chr   start   end start_zBased
#    <chr> <int> <int>        <int>
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

# unique(reference_tiles_2$chr)

##### read in site files and parse

# read in CG/CHG/CHH sites reference file created above
reference_tiles_C_sites <- read_tsv(paste0(genome, "_100bp_tiles_sites_Ns.tsv"), col_names = TRUE,
                                    cols(
                                      seqID = col_character(),
                                      patternName = col_character(),
                                      frequency = col_double()
                                    ))

# parse to make a bed file
reference_tiles_C_sites_bed <- reference_tiles_C_sites %>%
  #head(n=1000) %>%
  separate(seqID, into = c("chr", "start", "end"), sep = c(":|-")) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  mutate(start = ifelse(start >= 1, start+2, 0), end = end-2) %>%
  spread(key = patternName, value = frequency) %>%
  arrange(chr, start) %>%
  rename(cg_sites = CG,
         chg_sites = CHG,
         chh_sites = CHH
  )

reference_tiles_C_sites_bed

# # A tibble: 7,097,154 x 7
#    chr   start   end cg_sites chg_sites chh_sites     N
#    <chr> <dbl> <dbl>    <dbl>     <dbl>     <dbl> <dbl>
#  1 Chr00     0   100       12         8        34    NA
#  2 Chr00   100   200       NA        NA        18    NA
#  3 Chr00   200   300        4         3         8    NA
#  4 Chr00   300   400       10         2        26    NA
#  5 Chr00   400   500       10         4        26    NA
#  6 Chr00   500   600       12         2        27    NA
#  7 Chr00   600   700        8         2        25    NA
#  8 Chr00   700   800        4         4        23    NA
#  9 Chr00   800   900       NA         4        26    NA
# 10 Chr00   900  1000       NA         4        31    NA

# unique(reference_tiles_C_sites_bed$chr)

# replace NAs with 0

reference_tiles_C_sites_bed <- reference_tiles_C_sites_bed %>% replace(., is.na(.), 0)

reference_tiles_C_sites_bed %>% group_by(N) %>% summarise(tile_Ns = n())

# # A tibble: 101 x 2
#        N tile_Ns
#    <dbl>   <int>
#  1     0 6247638
#  2     2     225
#  3     4      38
#  4     6      40
#  5     8      39
#  6    10      43
#  7    12      38
#  8    14      49
#  9    16      37
# 10    18      44

reference_tiles_C_sites_bed %>% filter(N >100)

# # A tibble: 847,188 x 7
#    chr    start    end cg_sites chg_sites chh_sites     N
#    <chr>  <dbl>  <dbl>    <dbl>     <dbl>     <dbl> <dbl>
#  1 Chr00 133900 134000        0         0         0   196
#  2 Chr00 134000 134100        0         0         0   200
#  3 Chr00 134100 134200        0         0         0   200
#  4 Chr00 134200 134300        0         0         0   200
#  5 Chr00 134300 134400        0         0         0   200
#  6 Chr00 134400 134500        0         0         0   200
#  7 Chr00 134500 134600        0         0         0   200
#  8 Chr00 134600 134700        0         0         0   200
#  9 Chr00 134700 134800        0         0         0   200
# 10 Chr00 134800 134900        0         0         0   200

# There are 0.85 million 100 bp tiles in Apple that have > 50% N
# 6.2 million tiles have 0 Ns out of 7.1 million

#merge
reference_tiles_sites <- left_join(reference_tiles_2, reference_tiles_C_sites_bed, by = c("chr", "start_zBased" ="start", "end"))
reference_tiles_sites

# # A tibble: 7,097,206 x 8
#    chr   start   end start_zBased cg_sites chg_sites chh_sites     N
#    <chr> <int> <dbl>        <dbl>    <dbl>     <dbl>     <dbl> <dbl>
#  1 Chr11     1   100            0        2         6        30     0
#  2 Chr11   101   200          100        6         3        29     0
#  3 Chr11   201   300          200        2         3        32     0
#  4 Chr11   301   400          300        4         3        33     0
#  5 Chr11   401   500          400        2         6        27     0
#  6 Chr11   501   600          500        0         8        31     0
#  7 Chr11   601   700          600        2         5        21     0
#  8 Chr11   701   800          700        2         8        26     0
#  9 Chr11   801   900          800        2         6        25     0
# 10 Chr11   901  1000          900        2        11        37     0

# where are the tiles with no data?
reference_tiles_sites %>% filter(is.na(N))

# # A tibble: 71 x 8
#    chr      start      end start_zBased cg_sites chg_sites chh_sites     N
#    <fct>    <int>    <dbl>        <dbl>    <dbl>     <dbl>     <dbl> <dbl>
#  1 Chr11 28924701 28924800     28924700       NA        NA        NA    NA
#  2 Chr11 36902801 36902900     36902800       NA        NA        NA    NA
#  3 Chr11 38966101 38966200     38966100       NA        NA        NA    NA
#  4 Chr11 43059701 43059800     43059700       NA        NA        NA    NA
#  5 Chr10 22963101 22963200     22963100       NA        NA        NA    NA
#  6 Chr10 33836601 33836700     33836600       NA        NA        NA    NA
#  7 Chr10 33836701 33836800     33836700       NA        NA        NA    NA
#  8 Chr10 33836801 33836900     33836800       NA        NA        NA    NA
#  9 Chr10 33836901 33837000     33836900       NA        NA        NA    NA
# 10 Chr10 33837001 33837100     33837000       NA        NA        NA    NA

# grep -B 2 -A 3 'Chr11:28924698-28924802' GDDH13_1-1_formatted_Pt_100bp_tiles.fa
# >Chr11:28924698-28924802
# TAAAAAATAAATAAATAAAAATAAATAAAATTAAATTAAATTAAATTAAATTAAAAAATAATAAATAATAAATAAATATAAATAAATAATAAATAATAAATAAA

# Ok going to make the assumption that if a base is not an C or an N, then it must be an A,T, or G...
# 71 tiles with no definitely no Cs

# re arrange and replace NAs with 0
reference_tiles_sites <- reference_tiles_sites %>%
  replace(., is.na(.), 0) %>%
  mutate(chr = factor(chr, levels = unique(chr))) %>%
  arrange(chr, start)

# unique(reference_tiles_sites$chr)

#write the one-based txt file
write.table(reference_tiles_sites, paste0(genome, "_100bp_tiles_zBased_sites_counts_Ns.txt"), sep = "\t", quote = F, row.names = F)

##############

#write a bed file (zero based and no head)
reference_tiles_2_bed <- select(reference_tiles_sites, chr, start_zBased, end, cg_sites, chg_sites, chh_sites, N)

write.table(reference_tiles_2_bed, paste0(genome, "_100bp_tiles_zBased_sites_counts_Ns.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

##############
##############
##############
# summary plots

outFolder = "sites_analysis/"

dir.create(outFolder)

############## summary stats

summary_distro <- reference_tiles_2_bed %>%
  # slice(1:10000) %>%
  mutate(tile = 1:n()) %>%
  select(tile, cg_sites:N) %>%
  gather(key = motif, value = sites, -tile) %>%
  group_by(motif, sites) %>%
  summarise(freq = n())

summary_distro

write.table(summary_distro, paste0(genome, "_sites_counts_Ns_summary_distro.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


summary_distro %>% filter(!motif == "N")

g <- ggplot(summary_distro, aes(x = sites, y = freq, colour = motif)) +
  geom_line() +
  theme_minimal()

# g

pdf(paste0(outFolder, "/summary_distro_line.pdf"), h = 3, w = 8)
print(g)
dev.off()

plot_data <- summary_distro %>%
  # filter(!motif == "N") %>%
  mutate(bin = ceiling(sites / 10)) %>%
  mutate(sites_range = ifelse(bin == 0, paste0("0:0"), paste0(bin*10-9, ":", bin*10))) %>%
  group_by(bin) %>%
  group_by(motif, sites_range) %>%
  summarise(freq = sum(freq)) %>%
  ungroup()

ftr_lvl <- c(paste0(c(0, seq(1, 200, by = 10)), ":", c(seq(0, 200, by = 10))))

# ftr_lvl <- plot_data %>% select(bin, sites_range) %>% distinct(bin, sites_range) %>% pull(sites_range)

plot_data <- plot_data %>%
  mutate(sites_range = factor(sites_range, levels = ftr_lvl))

g <- ggplot(plot_data, aes(x = motif, y = freq, fill = sites_range)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "B", begin = 0.2) +
  theme_minimal()

# g

pdf(paste0(outFolder, "/summary_distro_stack_bar.pdf"), h = 4, w = 5)
print(g)
dev.off()


reference_tiles <- read_tsv(paste0(genome, "_100bp_tiles_zBased_sites_counts.txt"), col_names = TRUE,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer(),
                              cg_sites = col_integer(),
                              chg_sites = col_integer(),
                              chh_sites = col_integer()
                            ))

reference_tiles[is.na(reference_tiles)] <- 0

# reference_tiles %>% group_by(cg_sites) %>% summarise(cg_Cs = n())  %>%
# ungroup() %>% mutate(cum_percent = cumsum(n_sites)/sum(n_sites)*100) %>%
# mutate(remaining_percent = 100 -lag(cum_percent)) %>% print(n = nrow(.))

gather_Cs <- reference_tiles %>%
  mutate(tile_ID = paste0(chr, "_", start)) %>%
  select(tile_ID, cg_sites, chg_sites, chh_sites) %>%
  gather(key = context, value = Cs, -tile_ID)

gather_Cs

gather_Cs_summary <- gather_Cs %>% group_by(context, Cs) %>% summarise(n_Cs = n()) %>% ungroup()

gather_Cs_summary

g <- ggplot(gather_Cs_summary, aes(x = factor(Cs), y = n_Cs, fill = context)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_minimal() +
  scale_fill_fivethirtyeight() +
  # xlim(c(0, 50)) +
  labs(x = "number of Cs", y = "Number of tiles") +
  NULL

# print(g)

pdf(paste0(outFolder, "Cs_per_context_per_tile.pdf"), h = 4, w = 5)
print(g)
dev.off()

######################
# check distribution of zero site among tiles

reference_tiles_zeros <- reference_tiles %>%
  mutate(CG = ifelse(cg_sites > 0, "CG", ""),
         CHG = ifelse(chg_sites > 0, "CHG", ""),
         CHH = ifelse(chh_sites > 0, "CHH", ""))

reference_tiles_zeros_summary <- reference_tiles_zeros %>%
  mutate(tile_ID = paste0(chr, "_", start),
         tile_sites = paste0(CG, ":", CHG, ":", CHH)) %>%
  group_by(tile_sites) %>%
  summarise(n = n())

reference_tiles_zeros_summary <- reference_tiles_zeros_summary %>%
  mutate(percentage = round(n / sum(n)*100, digits = 2))

reference_tiles_zeros_summary_plot <- reference_tiles_zeros_summary %>% arrange(desc(percentage)) %>% mutate(tile_sites = factor(tile_sites, levels = tile_sites))

g <- ggplot(reference_tiles_zeros_summary_plot, aes(x = tile_sites, y = percentage)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(y = "Percent of 100bp tiles in genome",
       x = NULL,
       fill="Methylation Domain") +
  geom_text(aes(label=percentage, vjust=-0.25), size = label_size) +
  text_size_theme_8_labels

g

ggsave(paste0(outFolder, "/zero_site_distribution.pdf"), h=4, w = 3)


######################
# check distribution of zero site among tiles


reference_tiles_zeros <- reference_tiles %>%
  mutate(CG = ifelse(cg_sites >= 2, "CG", ""),
         CHG = ifelse(chg_sites >= 2, "CHG", ""),
         CHH = ifelse(chh_sites >= 2, "CHH", ""))

reference_tiles_zeros_summary <- reference_tiles_zeros %>%
  mutate(tile_ID = paste0(chr, "_", start),
         tile_sites = paste0(CG, ":", CHG, ":", CHH)) %>%
  group_by(tile_sites) %>%
  summarise(n = n())

reference_tiles_zeros_summary <- reference_tiles_zeros_summary %>%
  mutate(percentage = round(n / sum(n)*100, digits = 2))

reference_tiles_zeros_summary_plot <- reference_tiles_zeros_summary %>% arrange(desc(percentage)) %>% mutate(tile_sites = factor(tile_sites, levels = tile_sites))

g <- ggplot(reference_tiles_zeros_summary_plot, aes(x = tile_sites, y = percentage)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(y = "Percent of 100bp tiles in genome",
       x = NULL,
       fill="Methylation Domain") +
  geom_text(aes(label=percentage, vjust=-0.25), size = label_size) +
  text_size_theme_8_labels

# g

ggsave(paste0(outFolder, "/two_or_more_site_distribution.pdf"), h=4, w = 3)


######################
# check distribution of zero site among tiles


reference_tiles_zeros <- reference_tiles %>%
  mutate(CG = ifelse(cg_sites >= 4, "CG", ""),
         CHG = ifelse(chg_sites >= 4, "CHG", ""),
         CHH = ifelse(chh_sites >= 4, "CHH", ""))

reference_tiles_zeros_summary <- reference_tiles_zeros %>%
  mutate(tile_ID = paste0(chr, "_", start),
         tile_sites = paste0(CG, ":", CHG, ":", CHH)) %>%
  group_by(tile_sites) %>%
  summarise(n = n())

reference_tiles_zeros_summary <- reference_tiles_zeros_summary %>%
  mutate(percentage = round(n / sum(n)*100, digits = 2))

reference_tiles_zeros_summary_plot <- reference_tiles_zeros_summary %>% arrange(desc(percentage)) %>% mutate(tile_sites = factor(tile_sites, levels = tile_sites))

g <- ggplot(reference_tiles_zeros_summary_plot, aes(x = tile_sites, y = percentage)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(y = "Percent of 100bp tiles in genome",
       x = NULL,
       fill="Methylation Domain") +
  geom_text(aes(label=percentage, vjust=-0.25), size = label_size) +
  text_size_theme_8_labels

# g

ggsave(paste0(outFolder, "/four_or_more_site_distribution.pdf"), h=4, w = 3)
