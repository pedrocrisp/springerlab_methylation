#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=16gb
#PBS -N trim_galore
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
echo ------------------------------------------------------

########## Modules #################

module load fastqc/0.11.5
module load cutadapt/1.8.1

########## Set up dirs #################
#cd into work dir
cd "$PBS_O_WORKDIR"

#make trimmed folder
trimmedfolder=analysis/trimmed
mkdir -p $trimmedfolder

fastqcfolder=analysis/fastqc
mkdir -p $fastqcfolder

########## Run #################
trim_galore \
--phred33 \
--fastqc \
--fastqc_args "--noextract --outdir $fastqcfolder" \
-o $trimmedfolder --paired reads/${ID}_R1_001.fastq reads/${ID}_R2_001.fastq


