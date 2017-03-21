#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=32gb
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

module load python2/2.7.8
module load java
module load bedtools
module load bamtools
#bsmap requires samtools < 1.0.0
module load samtools/0.1.18

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p BSMAPratio
mkdir -p TempOut
mkdir -p OnTargetCoverage

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
        
        # extract methylation information using bsmap tool methratio.py
        python /home/springer/pcrisp/software/bsmap-2.90/methratio.py \
        -o BSMAPratio/${ID}_methratio.txt \
        -d ${genome_reference} \
        -u \
        -z \
        -r bsmapped/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
        
        #awk funciton for extracting methylation info from methratio.py output. Check with Qing what this is meant to do. Also try to figure out how to split this over multiple lines
        #awk '(NR>1){if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' BSMAPratio/${ID}.txt > BSMAPratio/${ID}_BSMAP_out.txt
        awk_make_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        printf $1, $2-1, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/))
                        print $1, $2-1, $2, $3, "CHG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/))
                        print $1, $2-1, $2, $3, "CHH", $5, $6, $7, $8, $9, $10, $11, $12;
                else
                        print $1, $2-1, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                ' 
        awk -F $'\t' "$awk_make_bed" BSMAPratio/${ID}_methratio.txt > BSMAPratio/${ID}_BSMAP_out.txt
        
        # conversion rate
        # awk -F"\t" '{if($1=="Pt") print}' "./BSMAPratio/"${ID}"_BSMAP_out.txt" | awk '{sum1 += $8; sum2 +=$9} END {print sum1"\t"sum2"\t"100-sum1/sum2*100}' > "./ConversionRate/"$i"_conversion_rate.txt"
        
        ## summarize read count, only use uniquley mapped and properly paried reads
        ## args to be passed from qsub script ${intersect_regions_ref} eg ${refdir}/BSseqcapv2_specific_regions.bed
        
        bedtools intersect \
        -abam bsmapped/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
        -b ${intersect_regions_ref} \
        -bed \
        -wa \
        -wb > TempOut/${ID}_specific_region_pairs_clipOverlap.txt
        
        awk -F"\t" '{sum[$16]++} END {for (i in sum) print i"\t"sum[i]}' \
        TempOut/${ID}_specific_region_pairs_clipOverlap.txt > \
        OnTargetCoverage/${ID}_specific_region_count_pairs_clipOverlap.txt
        
        # count methylated and unmethylated Cytosine per region
        bedtools intersect \
        -a BSMAPratio/${ID}_BSMAP_out.txt \
        -b ${intersect_regions_ref} \
        -wa \
        -wb > TempOut/${ID}_BSMAP_out_ontarget.txt
        
        awk \
        -F"\t" '{mC[$14"\t"$15"\t"$16"\t"$5] += $8; CT[$14"\t"$15"\t"$16"\t"$5] += $9; n[$14"\t"$15"\t"$16"\t"$5]++} END { for (j in mC) print j"\t"n[j]"\t"mC[j]"\t"CT[j]}' \
        TempOut/${ID}_BSMAP_out_ontarget.txt > OnTargetCoverage/${ID}_BSMAP_out_ontarget_mC.txt