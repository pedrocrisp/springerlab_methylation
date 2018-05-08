#!/bin/bash -l
#PBS -l walltime=24:00:00,nodes=1:ppn=8,mem=42gb
#PBS -N bismark
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

########## QC #################
set -xeuo pipefail

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: array_ID is ${PBS_ARRAYID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################

module load bowtie2/2.3.0
module load samtools/1.7

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p bismark_mapping
mkdir -p bismark_mC_summaries

########## Run #################

########## bismark alignment #######
# https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#appendix-ii-bismark
# --multicore 2 run to instances of bismark to speed things up
#(please note that a typical Bismark run will use several cores already
#(Bismark itself, 2 or 4 threads for Bowtie/Bowtie2, Samtools, gzip etc...) and ~10-16GB of memory per thread depending on the choice of aligner and genome.
# WARNING: Bismark Parallel is resource hungry!)
# # --basename ${ID} Specifying --basename in conjuction with --multicore is currently not supported (but we are aiming to fix this soon). Please lose either --basename or --multicore to proceed

bismark \
--bowtie2 \
--output_dir bismark_mapping \
--multicore 2 \
--genome ${genome_reference} \
-1 ${read_folder}/${ID}_R1_001_val_1.fq \
-2 ${read_folder}/${ID}_R2_001_val_2.fq

### rename
# bam
mv bismark_mapping/${ID}_R1_001_val_1_bismark_bt2_pe.bam \
bismark_mapping/${ID}_bismark_bt2_pe.bam
# report
mv bismark_mapping/${ID}_R1_001_val_1_bismark_bt2_PE_report.txt \
bismark_mapping/${ID}_bismark_bt2_PE_report.txt

########## extract methylation #######
# This step extract methylation values
# makes a bedGraph
# makes a genome coverage files
# --include_overlap this includes the overlap portion of PE reads, this will be useful for call read specific methylation but should be disabled if calling DMRs etc
# Please note that a typical process of extracting a BAM file and writing out .gz output streams will in fact use ~3 cores per value of --multicore <int> specified
# the main memory sort buffer when sorting the methylation information, default 2G
# --cytosine_report genome-wide methylation report for all cytosines in the genome
#

bismark_methylation_extractor \
--gzip \
--multicore 2 \
--paired-end \
--bedGraph \
--buffer_size 10G \
--cytosine_report \
--include_overlap \
--genome_folder ${genome_reference} \
--output bismark_mC_summaries \
bismark_mapping/${ID}_bismark_bt2_pe.bam
