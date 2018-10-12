#!/bin/bash -l
#PBS -l walltime=30:00:00,nodes=1:ppn=8,mem=42gb
#PBS -N bsmapping
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

########## QC #################
set -xeuo pipefail

#add bsmap and associated scripts and samtools v0.1.18 (included in bsmap folder) to my path
PATH=~/software/bsmap-2.74:$PATH

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

# bsmap requires samtools < 1.0.0 (note: now output is sam and running samtools manually on my own)
# module load samtools/0.1.18 # no longer an available module
PATH=~/software/bsmap-2.74/samtools:$PATH

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p bsmapped

########## Run #################

# align adapter trimmed datasets to B73 genome
        # -r 0: [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default:1
        # -v 5: allow 5 mismatches (could also use -v 0.05 = 5% of read length)
        # -p 8: 8 threads/cores
        # -q 20: trim to q20

bsmap \
-a trimmed/${ID}_R1_001_val_1.fq \
-b trimmed/${ID}_R2_001_val_2.fq \
-d ${genome_reference} \
-o bsmapped/${ID}.sam \
-v 5 \
-r 0 \
-p 8 \
-q 20 \
-A $adapter_seq

echo total reads in sam

# consider calculating how many reads were used for mapping...
# samtools view -c bsmapped/${ID}.sam
