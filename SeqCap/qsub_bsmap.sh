#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=64gb
#PBS -N 20170313_testpipeline_SeqCap_1_Mei
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

set -x
set -e

module load samtools

cd /scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/analysis

bsmap -a trimmed/F1-16_Index5_S1_R1_001_val_1.fq -b trimmed/F1-16_Index5_S1_R2_001_val_2.fq -d /home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa -o bsmaped/F1-16_Index5_S1/F1-16_Index5_S1.bam -v 5 -r 0 -p 8 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

