#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=16gb
#PBS -N bsmapping_filter_summarise
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

#module load python2/2.7.8
module load java
#module load bedtools
module load bamtools
#bsmap requires samtools < 1.0.0
#module load samtools/0.1.18

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

cd analysis

########## Run #################

# remove PCR duplicates, must be sorted by coordinate using pickard
        java -jar /home/springer/pcrisp/software/picard.jar SortSam \
        INPUT=bsmapped/${ID}.bam \
        OUTPUT=bsmapped/${ID}_sorted.bam \
        SORT_ORDER=coordinate
        
        java -jar /home/springer/pcrisp/software/picard.jar MarkDuplicates \
        I=bsmapped/${ID}_sorted.bam \
        O=bsmapped/${ID}_sorted_MarkDup.bam \
        METRICS_FILE=bsmapped/${ID}_MarkDupMetrics.txt \
        REMOVE_DUPLICATES=true

        #DEPRECIATED
        # on-target HsMetrics
        # args to be passed from qsub script
        # ${CalculateHsMetrics_reference} eg ${refdir}/seqcapv2_onTarget-for-picard.bed
        #java -jar /home/springer/pcrisp/software/picard.jar CalculateHsMetrics \
        #I=bsmapped/$(ID}_sorted_MarkDup.bam \
        #O=bsmapped/$(ID}_HsMetrics_noDuplicate.txt \
        #BI=${CalculateHsMetrics_reference} \
        #TARGET_INTERVALS=${CalculateHsMetrics_reference} 
        
        # on-target CollectHsMetrics using pickard
        java -jar /home/springer/pcrisp/software/picard.jar CollectHsMetrics \
        I=bsmapped/${ID}_sorted_MarkDup.bam \
        O=bsmapped/${ID}_HsMetrics_noDuplicate.txt \
        CLIP_OVERLAPPING_READS=false \
        R=${genome_reference} \
        BAIT_INTERVALS=${CalculateHsMetrics_reference} \
        TARGET_INTERVALS=${CalculateHsMetrics_reference} 
        
        # keep properly paired reads using mabtools package
        bamtools filter \
        -isMapped true \
        -isPaired true \
        -isProperPair true \
        -in bsmapped/${ID}_sorted_MarkDup.bam \
        -out bsmapped/${ID}_sorted_MarkDup_pairs.bam
        
        # clip overlapping reads using bamUtils package
        bam clipOverlap \
        --in  bsmapped/${ID}_sorted_MarkDup_pairs.bam \
        --out bsmapped/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
        --stats
       
