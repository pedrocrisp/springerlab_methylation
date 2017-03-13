#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/logs/20170313_testpipeline_SeqCap_1_Mei_o
#PBS -e /scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/logs/20170313_testpipeline_SeqCap_1_Mei_e
#PBS -N 20170313_testpipeline_SeqCap_1_Mei
#PBS -r n
#PBS -M pcrisp@umn.edu

#args passed
workingdir=$1
log_folder=$2 
fastqcfolder=$3
trimmedfolder=$4
alignedfolder=$5
$BSMAPratio_folder=$6
TempOut=$7

#get job ID - CHECK
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load python2/2.7.8
module load samtools
module load java
module load bedtools
module load cutadapt/1.8.1
module load bamtools

#cd to working dir
cd $workingdir

###################

# trim adapter
        trim_galore \
        --phred33 \
        --fastqc \
        --fastqc_args "--noextract --outdir $fastqcfolder" \
        -o $trimmedfolder --paired "${workingdir}/${ID}_R2_001.fastq" "${workingdir}/${ID}_R2_001.fastq"

        # align adapter trimmed datasets to B73 genome
        bsmap \
        -a "${trimmedfolder}/${ID}_R1_001_val_1.fq" \
        -b "${trimmedfolder}/${ID}_R2_001_val_2.fq" \
        -d /work/02297/qingli/db/AGPV4_20170119/Zea_mays.AGPv4.dna.toplevel.fa \
        -o "${alignedfolder}/${ID}.bam" \
        -v 5 \
        -r 0 \
        -p 8 \
        -q 20 \
        -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

        # reomve PCR duplicates, must be sorted by coordinate
        java -jar /home/springer/pcrisp/software/picard.jar SortSam \
        "INPUT=${alignedfolder}/${ID}.bam" \
        "OUTPUT=${alignedfolder}/${ID}_sorted.bam" \
        SORT_ORDER=coordinate
        
        java -jar /home/springer/pcrisp/software/picard.jar MarkDuplicates \
        "I=${alignedfolder}/${ID}_sorted.bam" \
        "O=${alignedfolder}/${ID}_sorted_MarkDup.bam" \
        "METRICS_FILE=.${alignedfolder}/${ID}_MarkDup.txt" \
        REMOVE_DUPLICATES=true

        ## on-target HsMetrics
        java -jar /home/springer/pcrisp/software/picard.jar CalculateHsMetrics \
        "I=./BAM/"$i"_sorted_MarkDup.bam" \
        "O=./HsMetrics/"$i"_HsMetrics_noDuplicate.txt" \
        BI=./REFfile/seqcapv2_onTarget-for-picard.bed \
        TARGET_INTERVALS=./REFfile/seqcapv2_onTarget-for-picard.bed

        # keep properly paired reads
        bamtools filter \
        -isMapped true \
        -isPaired true \
        -isProperPair true \
        -in "${alignedfolder}/${ID}_sorted_MarkDup.bam" \
        -out "${alignedfolder}/${ID}_sorted_MarkDup_pairs.bam"
        
        # clip overlapping reads
        bamtools clipOverlap \
        --in  "${alignedfolder}/${ID}_sorted_MarkDup_pairs.bam" \
        --out "${alignedfolder}/${ID}_sorted_MarkDup_pairs_clipOverlap.bam" \
        --stats

        # extract methylation information
        python /home/springer/pcrisp/software/bsmap-2.90/methratio.py \
        -o "$BSMAPratio_folder/${ID}" \
        -d /work/02297/qingli/db/AGPV4_20170119/Zea_mays.AGPv4.dna.toplevel.fa \
        -u \
        -z \
        -r "${alignedfolder}/${ID}_sorted_MarkDup_pairs_clipOverlap.bam"
        
        awk '(NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; 
                else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; 
                else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; 
                else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12
                }' \
                "$BSMAPratio_folder/${ID}" > "$BSMAPratio_folder/${ID}_BSMAP_out.txt"

        # conversion rate
        # awk -F"\t" '{if($1=="Pt") print}' "./BSMAPratio/"$i"_BSMAP_out.txt" | awk '{sum1 += $8; sum2 +=$9} END {print sum1"\t"sum2"\t"100-sum1/sum2*100}' > "./ConversionRate/"$i"_conversion_rate.txt"


        ## summarize read count, only use uniquley mapped and properly paried reads
        bedtools intersect \
        -abam "${alignedfolder}/${ID}_sorted_MarkDup_pairs_clipOverlap.bam" \
        -b ./REFfile/BSseqcapv2_specific_regions.bed \
        -bed \
        -wa \
        -wb > "${TempOut}/${ID}_specific_region_pairs_clipOverlap.txt"
        
        awk -F"\t" '{sum[$16]++} END {for (i in sum) print i"\t"sum[i]}' \
        "${TempOut}/${ID}_specific_region_pairs_clipOverlap.txt" > \
        "${OnTargetCoverage}/${ID}_specific_region_count_pairs_clipOverlap.txt"

        # count methylated and unmethylated Cytosine per region
        bedtools intersect \
        -a "$BSMAPratio_folder/${ID}_BSMAP_out.txt" \
        -b ./REFfile/BSseqcapv2_specific_regions.bed \
        -wa \
        -wb > "${TempOut}/${ID}_BSMAP_out_ontarget.txt"
        
        awk \
        -F"\t" '{mC[$14"\t"$15"\t"$16"\t"$5] += $8; CT[$14"\t"$15"\t"$16"\t"$5] += $9; n[$14"\t"$15"\t"$16"\t"$5]++} END { for (j in mC) print j"\t"n[j]"\t"mC[j]"\t"CT[j]}' "${TempOut}/${ID}_BSMAP_out_ontarget.txt" > "${OnTargetCoverage}/${ID}_BSMAP_out_ontarget_mC.txt"
