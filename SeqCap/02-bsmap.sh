#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=8gb
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
echo PBS: array_ID is ${PBS_ARRAYID}
echo ------------------------------------------------------

echo working dir is $PWD

#cd into work dir
echo changing to PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
echo working dir is now $PWD

########## Modules #################

module load samtools/1.3

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
mapfolder=analysis/bsmapped
mkdir -p $mapfolder
trimmed_folder=analysis/trimmed

########## Run #################

# align adapter trimmed datasets to B73 genome
        # -r 0: dont report repeat hits
        # -v 5: allow 5 mismatches (could also use -v 0.05 = 5% of read length)
        # -p 8: 8 threads/cores
        # -q 20: trim to q20

bsmap \
-a analysis/trimmed/F1-16_Index5_S1_R1_001_val_1.fq \
-b analysis/trimmed/F1-16_Index5_S1_R2_001_val_2.fq \
-d $genome_reference \
-o analysis/bsmapped/F1-16_Index5_S1.bam \
-v 5 \
-r 0 \
-p 1 \
-q 20 \
-A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG


