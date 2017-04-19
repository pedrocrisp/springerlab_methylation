#!/bin/bash -l
#PBS -l walltime=03:00:00,nodes=1:ppn=1,mem=24gb
#PBS -N fix-sort
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

#bsmap requires samtools < 1.0.0
module load samtools/0.1.18

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

cd analysis/bsmapped

########## Run #################

# make bam
samtools view -bS ${ID}.sam > ${ID}.bam

# sort by read name (needed for fixsam)
samtools sort -n ${ID}.bam ${ID}_nameSrt

# fix mate pairs
samtools fixmate ${ID}_nameSrt.bam ${ID}_nameSrt_fixed.bam

# co-ordinate sort
samtools sort ${ID}_nameSrt_fixed.bam ${ID}_sorted

# index
samtools index ${ID}_sorted.bam

# remove intermediate files
rm ${ID}.sam ${ID}_nameSrt.bam ${ID}_nameSrt_fixed.bam ${ID}.bam
