#!/bin/bash -l
#PBS -N 21-call-umrs
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

module load R/3.3.2
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID
echo genome is $genome_prefix

mkdir -p sites
cd sites

annotation_suffix=_mC_domains_II_cov_5_sites_2_MR_${MR_percent}_UMR_${UMR_percent}

########## Run MODULE 1 #################

# Run R module to creat 100bp tile bed file
R -f ~/gitrepos/springerlab_methylation/SeqCap/21-call-umrs.R \
--args $reference_100bp_tiles $sample_list $annotation_suffix $chrom_sizes_path $coverage_filter_min $site_filter_min $MR_percent $UMR_percent

echo finished summarising
