#!/bin/bash

#SBATCH -J seqcapv2                        # job name
#SBATCH -o seqcapv2.%j.out        # output and error file name (%j expands to jobID)
#SBATCH -N 1                            # number of nodes requested
#SBATCH -n 24                                   # total number of tasks requested
#SBATCH -p normal           # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00         # run time (hh:mm:ss)
#SBATCH --mail-user=cauliqing@163.com
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH -A Springer_Vaughn           # <-- Allocation name to charge job against

module load python/2.7.11
#bsmap requires samtools < 1.0.0
module load samtools/0.1.18
module load java/1.8.77
module load bedtools


for i in F4-16_F5-16_Index2_S35	F4-16_F5-16_Index4_S37	F4-16_F5-16_Index5_S38	F4-16_F5-16_Index6_S39	F4-16_F5-16_Index7_S40	F4-16_F5-16_Index12_S45	F4-16_F5-16_Index13_S46	F4-16_F5-16_Index14_S47	F4-16_F5-16_Index15_S48	F4-16_F5-16_Index16_S49	F4-16_F5-16_Index18_S50	F4-16_F5-16_Index19_S51
do
        # trim adapter
        # cd /corral-tacc/projects/iplant/vaughn/springer_vaughn/eichten/seqcap_jan_2017/
        /work/02297/qingli/bin/trim_galore_zip/trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir /scratch/02297/qingli/seqCap_V2/B73_V4/rawdata/fastqc" -o /scratch/02297/qingli/seqCap_V2/B73_V4/rawdata/fastqc --paired "/corral-tacc/projects/iplant/vaughn/springer_vaughn/eichten/seqcap_jan_2017/"$i"_R1_001.fastq" "/corral-tacc/projects/iplant/vaughn/springer_vaughn/eichten/seqcap_jan_2017/"$i"_R2_001.fastq"

        # align adapter trimmed datasets to B73 genome
        cd /scratch/02297/qingli/seqCap_V2/B73_V4
        /scratch/02297/qingli/bsmap-2.90_2/bsmap -a "./rawdata/fastqc/"$i"_R1_001_val_1.fq" -b "./rawdata/fastqc/"$i"_R2_001_val_2.fq" -d /work/02297/qingli/db/AGPV4_20170119/Zea_mays.AGPv4.dna.toplevel.fa -o "./align_to_B73/BAM/"$i".bam" -v 5 -r 0 -p 8 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

        # reomve PCR duplicates, must be sorted by coordinate
        cd /scratch/02297/qingli/seqCap_V2/B73_V4/align_to_B73
        java -jar /work/02297/qingli/bin/picard-tools-1.102/SortSam.jar "INPUT=./BAM/"$i".bam" "OUTPUT=./BAM/"$i"_sorted.bam" SORT_ORDER=coordinate
        java -jar /work/02297/qingli/bin/picard-tools-1.102/MarkDuplicates.jar "I=./BAM/"$i"_sorted.bam" "O=./BAM/"$i"_sorted_MarkDup.bam" "METRICS_FILE=./DuplicateRate/"$i"_MarkDup.txt" REMOVE_DUPLICATES=true


        ## on-target HsMetrics
        cd /scratch/02297/qingli/seqCap_V2/B73_V4/align_to_B73
        java -jar /work/02297/qingli/bin/picard-tools-1.102/CalculateHsMetrics.jar "I=./BAM/"$i"_sorted_MarkDup.bam"  "O=./HsMetrics/"$i"_HsMetrics_noDuplicate.txt" BI=./REFfile/seqcapv2_onTarget-for-picard.bed TARGET_INTERVALS=./REFfile/seqcapv2_onTarget-for-picard.bed


        # keep properly paired reads
        cd /scratch/02297/qingli/seqCap_V2/B73_V4/align_to_B73
        /scratch/02297/qingli/bamtools/bin/bamtools filter -isMapped true -isPaired true -isProperPair true -in "./BAM/"$i"_sorted_MarkDup.bam" -out "./BAM/"$i"_sorted_MarkDup_pairs.bam"
        # clip overlapping reads
        /scratch/02297/qingli/bamUtil/bin/bam clipOverlap --in  "./BAM/"$i"_sorted_MarkDup_pairs.bam" --out "./BAM/"$i"_sorted_MarkDup_pairs_clipOverlap.bam" --stats


        # extract methylation information
        python /scratch/02297/qingli/bsmap-2.74/methratio.py -o "./BSMAPratio/"$i -d /work/02297/qingli/db/AGPV4_20170119/Zea_mays.AGPv4.dna.toplevel.fa -u -z -r "./BAM/"$i"_sorted_MarkDup_pairs_clipOverlap.bam"
        awk '(NR>1){if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' "./BSMAPratio/"$i > "./BSMAPratio/"$i"_BSMAP_out.txt"

        # conversion rate
        # awk -F"\t" '{if($1=="Pt") print}' "./BSMAPratio/"$i"_BSMAP_out.txt" | awk '{sum1 += $8; sum2 +=$9} END {print sum1"\t"sum2"\t"100-sum1/sum2*100}' > "./ConversionRate/"$i"_conversion_rate.txt"


        ## summarize read count, only use uniquley mapped and properly paried reads
        cd /scratch/02297/qingli/seqCap_V2/B73_V4/align_to_B73
        bedtools intersect -abam "./BAM/"$i"_sorted_MarkDup_pairs_clipOverlap.bam" -b ./REFfile/BSseqcapv2_specific_regions.bed -bed -wa -wb > ./TempOut/$i"_specific_region_pairs_clipOverlap.txt"
        awk -F"\t" '{sum[$16]++} END {for (i in sum) print i"\t"sum[i]}' ./TempOut/$i"_specific_region_pairs_clipOverlap.txt" > ./OnTargetCoverage/$i"_specific_region_count_pairs_clipOverlap.txt"

        # count methylated and unmethylated Cytosine per region
        bedtools intersect -a "./BSMAPratio/"$i"_BSMAP_out.txt" -b ./REFfile/BSseqcapv2_specific_regions.bed -wa -wb > "./TempOut/"$i"_BSMAP_out_ontarget.txt"
        awk -F"\t" '{mC[$14"\t"$15"\t"$16"\t"$5] += $8; CT[$14"\t"$15"\t"$16"\t"$5] += $9; n[$14"\t"$15"\t"$16"\t"$5]++} END { for (j in mC) print j"\t"n[j]"\t"mC[j]"\t"CT[j]}' "./TempOut/"$i"_BSMAP_out_ontarget.txt" > "./OnTargetCoverage/"$i"_BSMAP_out_ontarget_mC.txt"
done


