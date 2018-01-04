# SeqCap pipeline

Pipeline for processing SeqCap data.

--------

## Installation Notes
This pipeline requires a number of pieces of software, some are loaded as modules, some are expected to be in your path, depending on the step.

It is designed to run on UMN MSI computer clusters.

Bsmap is a special case. It looks for bsmap in ~/software/bsmap-2.74 and puts this dir in your path. I also edited the script sam2bam.sh to use the system version of samtools; however, ultimately I dont use this file because I resorted to output in sam format because I fix the sam/bam due to a discordant pairs issues, then manually run samtools to sort and index bams.

To do things:
 - soft code more steps
 - described dependancies and required files structures

## Example pipeline execution - Server Steps

### Step 1 Trim reads

Use trimgalore to trim the reads.

usage="USAGE:
bash 01-trim_galore_qsub.sh <sample_list.txt>"

```
#01-trim_galore
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/01-trim_galore_qsub.sh \
samples.txt
```
**Methods summary**
Reads were trimmed and QC'ed with trim_galore version 0.4.3, powered by cutadapt v1.8.1 and fastqc v0.11.5.

### Step 2 Mapping using bsmap

Map using bsmap. This is preferred over bismark for speed with large genomes.

usage="USAGE:
bash 02-bsmap_qsub.sh <sample_list.txt> <genome.fa>"

```
#02-bsmap
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/02-bsmap_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa
```

**Methods summary**
Reads were aligned with bsmap v2.74 with the following parameters -v 5 to allow allow 5 mismatches, -r 0 to report only unique mapping pairs, -p 1, -q 20 to allow quality trimming to q20, -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG adapter sequence. Output file is in SAM format to allow custom QC and sorting (because some reads are not properly paired).

### Step 3 fix sam files

The output sam files produced by bsmap contain incorrect sam flags. Where PE reads map to different choromsomes the reads are still marked as correctly paird. This breaks picard tools, therefore, fix using samtools fixsam. Also make bams, sort and index (incase we want to vies these in IGV).

usage="USAGE:
bash 03-fix-sort-bsmap_qsub.sh <sample_list.txt>
"

```
#03-fix-sort-bsmap
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/03-fix-sort-bsmap_qsub.sh \
samples.txt
```

### Step 4 filter

Use picard tools to filter reads, removing duplicates, removing improperly paired reads and trimming overlapping reads so as to only count 'Cs' once per sequenced molecule. Also collect On target metrics.

usage="USAGE:
bash 04-filter_qsub.sh <sample_list.txt> <genome.fa> <CalculateHsMetrics_reference.bed>
for example:
"

```
#04-filter
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/04-filter_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/seqcapv2_onTarget-for-picard.bed
```

### Step 5 Extract methylation data

Use methylratio.py script from bsmap to extract methylation data. Then use awk to convert output to mC context (CG, CHG, CHH) and parse coordinate to zero-based format "BED" type for bedtools. Then use awk to convert to bedgraph format (chr, start, stop, ratio) and split into a separate file per context. Use bedGraphToBigWig to make a bigWig file for IGV. Use awk to get conversion rates using the chloroplast reads. Use bedtools to intersect bam with target regions then use awk to sum reads per region. Then use bedtools to intersect C and CT counts file with target regions and use awk to sum counts. Also count methylaiton per region using eff_CT incase we want this metric later.

usage="USAGE:
bash 05-summarise_methylation_qsub.sh <sample_list.txt> <genome.fa> <intersect_regions.bed>
for example:
"

```
#05-summarise_methylation
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/05-summarise_methylation_qsub.sh \
samples.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/BSseqcapv2_specific_regions.bed
```

output - *BSMAP_out.txt

chr     start end     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper
1       489     490     +       CHG     0.250   8.00    2       8       0       0       0.071   0.591
1       490     491     +       CG      0.875   8.00    7       8       0       0       0.529   0.978

### Step 7 100 bp tiles

Summarise methylation to 100 bp tiles across the genome. This is a 4-parter. First it runs the perl script `met_context_window.pl` to get C, CT, ratios and C_site_counts per 100bp window per context. Per sample any window with data is reported for all contexts even if some windows do not have data for a particular context. Next the R script `07-tiles_bed_to_bigWig.R` is called to fix the ends of the chromosomes which have windows the extend beyond the chromosome ends, this is necessary if we want to make bigWigs etc. This script also pulls in the C_sites per window information from the file `maize_v4_100pb_tiles_zBased_sites.txt`. This file was created by me, tiles windows across the genome, then merges with the file ` AGPv4_cg_sites_012017.bed` from Jawon c/o Qing c/o Jackie which gives C sites per 100 bp window. Note this C sites file does not include the contigs - might consider remaking this file to include the contigs (although they only account for 1.31% of the genome). Output is a bedgraph (fixed.sorted.bg; "chr", "start_zBased", "end", "ratio") and a text file with one-based coordinates (fixed.sorted.txt; "chr", "start", "end", "C", "CT", "ratio", "sites_with_data", "c_sites"). Then `bedGraphToBigWig` is used to make bigWigs from the tiles.

 usage="USAGE:
 bash 07-summarise_methylation_qsub.sh <sample_list.txt> <chrom.sizes file>
 for example:
 "

 Requires:
  - /home/springer/pcrisp/ws/refseqs/maize/maize_v4_100pb_tiles_zBased_sites.txt

```
#07-tiles_bed_to_bigWig
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/07-tiles_bed_to_bigWig_qsub.sh \
single_sample.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes
```

### 08-tiles_CHH_cov

This step takes the CHH output file from the 100 bp tiles script `*_BSMAP_out.txt.100.CHH.fixed.sorted.txt` and

1. Calculates CHH_cov (CHH_cov= CT/chh_sites).
2. Then using this field it filters on the arg $coverage_filter ($3); for example 2 or 5; which corrosponds to an average read coverage of 2 (or 5) accross all CHH sites in each 100 bp window, the output is significantly smaller file `*_BSMAP_out.txt.100.CHH_cov.txt`.

The idea (still being developed) is that:
1. These files can be used to explore the read coverage over the whole experiment; and
1. The tile index in this file could be combined with other samples to make a master list of tiles of interest (this is less memory intensive than combining the whole files).

usage="USAGE:
bash 08-tiles_analysis_qsub.sh <sample_list.txt> <data_folder> <coverage_filter>

For example:
```
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/08-tiles_analysis_CHH_cov_qsub.sh \
samples.txt \
tiles \
2
```
### 09-clean-up

Purges unnecesay files with the option to copy important files to home files space from scratch and also option to copy files to s3. Each of these three options (purge, copy to home and copy to s3) can be run independently.

See the script for more info and or the notebook.

### 10-tiles_filter_list

***Untested***

This step filters the `*_BSMAP_out.txt.100.CG.fixed.sorted.txt`, `*_BSMAP_out.txt.100.CHG.fixed.sorted.txt` and `*_BSMAP_out.txt.100.CHH.fixed.sorted.txt` files based on a list of tiles of interest. The idea is that once these smaller files are generated, they can be read into R in a stack for analysis all together.

For example:
```
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/10-tiles_analysis_filter_list_qsub.sh \
../samples.txt \
tiles \
../analysis_02_tiles_SeqCap_meta_140_samples/chh_2x_cov_62_sample_tile_list.tsv
```

--------

## Example pipeline execution - Laptop Steps

*It is possible these steps will be too intensive for the laptop with 100s of samples*

### Step 6 summarise-output

Use R to calculate ratios per region; add v4 annotations per region; and also merge with the read depth count file.

usage="USAGE:
bash 06-summarise-output-runner.sh <sample_list> <OnTargetCoverage_folder> <threads>
"

```
#06-summarise-output
bash ~/gitrepos/springerlab_methylation/SeqCap/06-summarise-output-runner.sh \
samples.txt \
OnTargetCoverage \
8
```

### Step 6b aggregate summarise-output

Use R to aggregate the summary output files into a single meta-table with all samples. Use this file to read into R and do further data exploration and calculate DMRs etc.

usage="USAGE:
Rscript 06-summarise-output-runner.sh <data_path> <outPrefix> <SeqCapEpi2_regions_annotation_v2_v4.csv_file>
"

```
Rscript \
~/gitrepos/springerlab_methylation/SeqCap/06b-aggregate-samples.R \
OnTargetCoverage_annotated \
Mei_final \
~/umn/refseqs/maize/SeqCap/Seqcap_ultimate_annotation_files/SeqCapEpi2_regions_annotation_v2_v4.csv
```

--------

# Scraping the log files

*Scrapers written so far*

1. Bsmap - to get initial the mapping rates
2. MarkDuplicates - to get duplication rates
3. Methratio - get final valid mapping used to extract methylation data

## Bsmap Scraper

```
cd logs/..._02-bsmap

for i in $(ls 02-bsmap_e*); do
SAMPLE=$(grep 'echo sample being mapped is' $i | cut -d " " -f 7)
MAPPED_PAIRS=$(grep 'pairs:' $i | cut -d " " -f 8)
MAPPED_PAIRS_PERCENTAGE=$(grep 'pairs:' $i | cut -d " " -f 9 | cut -d "(" -f2 | cut -d ")" -f1)
MAPPED_A=$(grep 'single a:' $i | cut -d " " -f 6)
MAPPED_A_PERCENTAGE=$(grep 'single a:' $i | cut -d " " -f 7 | cut -d "(" -f2 | cut -d ")" -f1)
MAPPED_B=$(grep 'single b:' $i | cut -d " " -f 6)
MAPPED_B_PERCENTAGE=$(grep 'single b:' $i | cut -d " " -f 7 | cut -d "(" -f2 | cut -d ")" -f1)
echo -e "$SAMPLE\t$MAPPED_PAIRS\t$MAPPED_PAIRS_PERCENTAGE\t$MAPPED_A\t$MAPPED_A_PERCENTAGE\t$MAPPED_B\t$MAPPED_B_PERCENTAGE"
done > bsmap_summary.tsv

```


## MarkDuplicates

```
cd /.../analysis/bsmapped_filtered

for i in $(cat ../../samples.txt); do
SAMPLE=$i
UNPAIRED_READS_EXAMINED=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt | cut -f 2)
STATS=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt)
UNPAIRED_READS_EXAMINED=$(grep 'Unknown Library' ${i}_MarkDupMetrics.txt | cut -f 2)
echo -e "$SAMPLE\t$STATS"
done > MarkDuplicates_scraped.tsv
```

## OnTargetMetrics

```
cd ../analysis/bsmapped_filtered

#Get header from a random sample
HEADER=$(grep 'BAIT_SET' US_1_Index9_S14_HsMetrics_noDuplicate.txt | head -1)
echo -e "Sample\t$HEADER" > OnTargetMetrics_scraped.tsv

#grep 'BAIT_SET' *_HsMetrics_noDuplicate.txt | head -1 > OnTargetMetrics_scraped.tsv

#scrape
for i in $(cat ../../samples.txt); do
SAMPLE=$i
STATS=$(grep '^SeqCapEpi2_v4_capture_space_sorted' ${i}_HsMetrics_noDuplicate.txt)
echo -e "$SAMPLE\t$STATS"
done >> OnTargetMetrics_scraped.tsv

```

## ConversionRate

```
cd ../analysis/ConversionRate

echo -e "Sample\tC_counts\tCT_counts\tConversionRate" > ConversionRate_scraped.tsv

for i in $(cat ../../samples.txt); do
SAMPLE=$i
File_data=${i}_conversion_rate.txt
#STATS=$(cat ${i}_conversion_rate.txt)
C_counts=$(cat $File_data | cut -d " " -f 1)
CT_counts=$(cat $File_data | cut -d " " -f 2)
ConversionRate=$(cat $File_data | cut -d " " -f 3)
echo -e "$SAMPLE\t$C_counts\t$CT_counts\t$ConversionRate"
done >> ConversionRate_scraped.tsv
```




## Methratio

```
cd logs/..._05-summarise_methylation

for i in $(ls 05-summarise_methylation_o*); do
SAMPLE=$(grep 'sample being mapped is' $i | cut -d " " -f 5)
VALID_MAPPINGS=$(grep 'total' $i | cut -d " " -f 2)
COVERED_CYTOSINES=$(grep 'total' $i | cut -d " " -f 5)
AVERAGE_COVERAGE=$(grep 'total' $i | cut -d " " -f 10)
echo -e "$SAMPLE\t$VALID_MAPPINGS\t$COVERED_CYTOSINES\t$AVERAGE_COVERAGE"
done > methratio_summary.tsv
```
