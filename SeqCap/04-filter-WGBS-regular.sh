#!/bin/bash -l
#PBS -N filter_WGBS_regular
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
PATH=~/software/bsmap-2.74/samtools:$PATH

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

cd analysis
mkdir -p bsmapped_filtered

if [ "$paired_end" == "yes" ]
then

########## Run #################

        # fix improperly paird reads - specifically discordant read pairs that are incorrectly mraked as concordant by bsmap
        # picard FixMateInformation
        # didnt work...
        #java -jar /home/springer/pcrisp/software/picard.jar FixMateInformation \
        #I=bsmapped/${ID}.bam \
        #O=bsmapped_filtered/${ID}_sorted.bam

        #remove PCR duplicates, must be sorted by coordinate using pickard

        # bams already co-ordinate sorted by samtools, this step seems unnecessary
        # Also causesing issues: bsmap is reporting PE reads as properly mapped where they hit different chromosomes, solution: skip step
        # uncomment to re-instate

        #java -jar /home/springer/pcrisp/software/picard.jar SortSam \
        #INPUT=bsmapped/${ID}.bam \
        #OUTPUT=bsmapped/${ID}_sorted.bam \
        #SORT_ORDER=coordinate

        #mark duplicates
        #requires sorted input - using samtools sort in bsmap step (co-ordinate sorted)
        # if co-ordinate sorted then pairs where the mate is unmapped or has secondary alignment are not marked as duplicate
        # ASSUME_SORTED=true because sorting performed with samtools but samtools doesnt seem to add this flag to the headder
        java -jar /home/springer/pcrisp/software/picard.jar MarkDuplicates \
        I=bsmapped/${ID}_sorted.bam \
        O=bsmapped_filtered/${ID}_sorted_MarkDup.bam \
        METRICS_FILE=bsmapped_filtered/${ID}_MarkDupMetrics.txt \
        ASSUME_SORTED=true \
        CREATE_INDEX=False \
        REMOVE_DUPLICATES=true

        # keep properly paired reads using bamtools package
        # note that some reads marked as properly paired by bsmap actually map to different chromosomes
        bamtools filter \
        -isMapped true \
        -isPaired true \
        -isProperPair true \
        -in bsmapped_filtered/${ID}_sorted_MarkDup.bam \
        -out bsmapped_filtered/${ID}_sorted_MarkDup_pairs.bam

        # clip overlapping reads using bamUtils package
        bam clipOverlap \
        --in bsmapped_filtered/${ID}_sorted_MarkDup_pairs.bam \
        --out bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
        --stats

        #index bam
        # index
        samtools index bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

elif [ "$paired_end" == "no" ]
then

# remove duplicate
java -jar /home/springer/pcrisp/software/picard.jar MarkDuplicates \
I=bsmapped/${ID}_sorted.bam \
O=bsmapped_filtered/${ID}_sorted_MarkDup.bam \
METRICS_FILE=bsmapped_filtered/${ID}_MarkDupMetrics.txt \
ASSUME_SORTED=true \
CREATE_INDEX=False \
REMOVE_DUPLICATES=true

# rename to keep script naming convention consistent
mv bsmapped_filtered/${ID}_sorted_MarkDup.bam bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

# index
samtools index bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam


else

echo "library type not specified correctly, please indicate yes for PE or no for SE"

fi
