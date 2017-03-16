#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=16gb
#PBS -N 20170313_testpipeline_SeqCap_1_Mei
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

#set -xe
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

cd "$PBS_O_WORKDIR"

#reference file dir
refdir=/home/springer/pcrisp/ws/refseqs/maize

#get job ID - CHECK
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

echo $ID

#load modules
#module load python2/2.7.8
#module load java
#module load bedtools
#module load cutadapt/1.8.1
#module load bamtools
#module load fastqc/0.11.5
module load samtools/1.3 #bsmap is not working, maybe this is the issue?

which samtools

samtools --version

###################
        
        mkdir "${alignfolder}/${ID}"
        #cd "${alignfolder}/${ID}"
        
        #bsmap \
        #-a trimmed/F1-16_Index5_S1_R1_001_val_1.fq \
        #-b trimmed/F1-16_Index5_S1_R2_001_val_2.fq \
        #-d "${refdir}/Zea_mays.AGPv4.dna.toplevel.fa" \
        #-o bsmaped/F1-16_Index5_S1/F1-16_Index5_S1.bam \
        #-v 5 \
        #-r 0 \
        #-p 1 \
        #-q 20 \
        #-A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
        
        bsmap \
        -a ${trimmedfolder}/${ID}_R1_001_val_1.fq \
        -b ${trimmedfolder}/${ID}_R2_001_val_2.fq \
        -d ${refdir}/Zea_mays.AGPv4.dna.toplevel.fa \
        -o ${alignfolder}/${ID}/${ID}.bam \
        -v 5 \
        -r 0 \
        -p 1 \
        -q 20 \
        -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
        
         #-o "bsmaped/${ID}/${ID}.bam" \
