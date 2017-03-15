#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=64gb
#PBS -N 20170313_testpipeline_SeqCap_1_Mei
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

set -x
set -e

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

#cd to working dir
cd $workingdir

###################
        
        mkdir "${alignfolder}/${ID}"
        #cd "${alignfolder}/${ID}"
        
        bsmap \
        -a trimmed/F1-16_Index5_S1_R1_001_val_1.fq \
        -b trimmed/F1-16_Index5_S1_R2_001_val_2.fq \
        -d "${refdir}/Zea_mays.AGPv4.dna.toplevel.fa" \
        -o bsmaped/F1-16_Index5_S1/F1-16_Index5_S1.bam \
        -v 5 \
        -r 0 \
        -p 1 \
        -q 20 \
        -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
        
        #bsmap \
        #-a "${trimmedfolder}/${ID}_R1_001_val_1.fq" \
        #-b "${trimmedfolder}/${ID}_R2_001_val_2.fq" \
        #-d "${refdir}/Zea_mays.AGPv4.dna.toplevel.fa" \
        #-o "${alignfolder}/${ID}/${ID}.bam" \
        #-v 5 \
        #-r 0 \
        #-p 1 \
        #-q 20 \
        #-A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

cd $workingdir
