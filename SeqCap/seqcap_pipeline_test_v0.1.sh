#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=16gb
#PBS -N 20170313_testpipeline_SeqCap_1_Mei
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

set -x
set -e

#dir structure
#readsdir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/reads
#workingdir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/analysis
#log_folder=${workingdir}/logs_${timestamp}
#fastqcfolder=${workingdir}/fastqc
#trimmed=${workingdir}/trimmed
#alignfolder=${workingdir}/bsmaped
#BSMAPratio_folder=${workingdir}/BSMAPratio
#TempOut=${workingdir}/TempOut
#OnTargetCoverage=${workingdir}/OnTargetCoverage

#reference file dir
refdir=/home/springer/pcrisp/ws/refseqs/maize

#get job ID - CHECK
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

echo $ID

#load modules
module load python2/2.7.8
module load java
module load bedtools
module load cutadapt/1.8.1
module load bamtools
module load fastqc/0.11.5
module load samtools/1.3 #bsmap is not working, maybe this is the issue?

#cd to working dir
cd $workingdir

###################

# trim adapter
        trim_galore \
        --phred33 \
        --fastqc \
        --fastqc_args "--noextract --outdir $fastqcfolder" \
        -o $trimmedfolder --paired "${readsdir}/${ID}_R1_001.fastq" "${readsdir}/${ID}_R2_001.fastq"

        # align adapter trimmed datasets to B73 genome
        # -r 0: dont report repeat hits
        # -v 5: allow 5 mismatches (could also use -v 0.05 = 5% of read length)
        # -p 8: 8 threads/cores
        # -q 20: trim to q20
        
        #-o "${alignfolder}/${ID}.bam" \
        #-o "./${ID}.bam" \
        
        which samtools
        
        mkdir "${alignfolder}/${ID}"
        cd "${alignfolder}/${ID}"
        
        bsmap \
        -a "${trimmedfolder}/${ID}_R1_001_val_1.fq" \
        -b "${trimmedfolder}/${ID}_R2_001_val_2.fq" \
        -d "${refdir}/Zea_mays.AGPv4.dna.toplevel.fa" \
        -o "${alignfolder}/${ID}/${ID}.bam" \
        -v 5 \
        -r 0 \
        -p 8 \
        -q 20 \
        -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

cd $workingdir

        # remove PCR duplicates, must be sorted by coordinate
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
        BI=${refdir}/seqcapv2_onTarget-for-picard.bed \
        TARGET_INTERVALS="${refdir}/seqcapv2_onTarget-for-picard.bed

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
        -d ${refdir}/Zea_mays.AGPv4.dna.toplevel.fa \
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
        -b .${refdir}/BSseqcapv2_specific_regions.bed \
        -bed \
        -wa \
        -wb > "${TempOut}/${ID}_specific_region_pairs_clipOverlap.txt"
        
        awk -F"\t" '{sum[$16]++} END {for (i in sum) print i"\t"sum[i]}' \
        "${TempOut}/${ID}_specific_region_pairs_clipOverlap.txt" > \
        "${OnTargetCoverage}/${ID}_specific_region_count_pairs_clipOverlap.txt"

        # count methylated and unmethylated Cytosine per region
        bedtools intersect \
        -a "$BSMAPratio_folder/${ID}_BSMAP_out.txt" \
        -b ${refdir}/BSseqcapv2_specific_regions.bed \
        -wa \
        -wb > "${TempOut}/${ID}_BSMAP_out_ontarget.txt"
        
        awk \
        -F"\t" '{mC[$14"\t"$15"\t"$16"\t"$5] += $8; CT[$14"\t"$15"\t"$16"\t"$5] += $9; n[$14"\t"$15"\t"$16"\t"$5]++} END { for (j in mC) print j"\t"n[j]"\t"mC[j]"\t"CT[j]}' "${TempOut}/${ID}_BSMAP_out_ontarget.txt" > "${OnTargetCoverage}/${ID}_BSMAP_out_ontarget_mC.txt"
