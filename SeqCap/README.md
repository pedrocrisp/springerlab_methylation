# SeqCap pipeline

Pipeline for processing SeqCap data.

--------

## Installation Notes
This pipeline requires a number of pieces of software, some are loaded as modules, some are expected to be in your path, depending on the step.

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
bash 04-filter_qsub.sh <sample_list.txt> <genome.fa> <CalculateHsMetrics_reference.bed> <intersect_regions.bed>
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

Use methylratio.py script from bsmap to extract methylation data. Then use awk to convert output to mC context and add a column for chromosome end column to bedtools. Then use awk to convert to bedgraph format and split into a separate file per context. Use bedGraphToBigWig to make a bigWig file for IGV. Use awk to get conversion rates using the chloroplast reads. Use bedtools to intersect bam with target regions then use awk to sum reads per region. Then use bedtools to intersect C and CT counts file with target regions and use awk to sum counts. Also count methylaiton per region using eff_CT incase we want this metric later.

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

--------

## Example pipeline execution - Laptop Steps

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
