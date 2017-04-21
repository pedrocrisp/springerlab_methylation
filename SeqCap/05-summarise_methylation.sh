#!/bin/bash -l
#PBS -l walltime=7:00:00,nodes=1:ppn=1,mem=60gb
#PBS -N summarise_methylation
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
#module load java
module load bedtools
#module load bamtools
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
mkdir -p ConversionRate

########## Run #################
        # The required input is all in the folder bsmapped_filtered from the prior step
        ########################
        # extract methylation information using bsmap tool methratio.py
        python /home/springer/pcrisp/software/bsmap-2.74/methratio.py \
        -o BSMAPratio/${ID}_methratio.txt \
        -d ${genome_reference} \
        -u \
        -z \
        -r bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam

        #awk funciton for extracting methylation info from methratio.py output. Check with Qing what this is meant to do. Also try to figure out how to split this over multiple lines
        #awk '(NR>1){if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' BSMAPratio/${ID}.txt > BSMAPratio/${ID}_BSMAP_out.txt

        #This awk makes the bed file at 1nt resolution rather than 2 as Qing did...
        awk_make_bed='BEGIN {OFS = FS} (NR>1){
                if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./))
                        print $1, $2, $2, $3, "CG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/))
                        print $1, $2, $2, $3, "CHG", $5, $6, $7, $8, $9, $10, $11, $12;
                else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/))
                        print $1, $2, $2, $3, "CHH", $5, $6, $7, $8, $9, $10, $11, $12;
                else
                        print $1, $2, $2, $3, "CNN", $5, $6, $7, $8, $9, $10, $11, $12
                }
                '
        #awk -F$"\t" "$awk_make_bed" "F1-16_Index5_S1_methratio.txt" > F1-16_Index5_S1.bed
        #BSMAPratio/$F1-16_Index5_S1_BSMAP_out.txt

        awk -F$"\\t" "$awk_make_bed" \
        "BSMAPratio/${ID}_methratio.txt" > "BSMAPratio/${ID}_BSMAP_out.txt"

        ########################
        #For genome browser

        # begGraph ratio files for tdfs
        awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $8/$9*100, $5
        }
        '

        # split by bedGraph by contex
        awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
          print $1, $2, $3, $4 > "BSMAPratio/"ID"_BSMAP_out_"$5".bedGraph"
        }
        '

        #pipe bedGraph to split by context (use dash to read from sdtin)
        awk -F$"\\t" "$awk_make_bedGraph" \
        "BSMAPratio/${ID}_BSMAP_out.txt" | \
        awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

        #Make bigWigs
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CG.bedGraph" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out_CG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHG.bedGraph" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out_CHG.bigWig"
        bedGraphToBigWig "BSMAPratio/${ID}_BSMAP_out_CHH.bedGraph" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out_CHH.bigWig"

        ########################

        # conversion rate
        #awk -F"\t" '{if($1=="Pt") print}' "./BSMAPratio/"${ID}"_BSMAP_out.txt" | awk '{sum1 += $8; sum2 +=$9} END {print sum1"\t"sum2"\t"100-sum1/sum2*100}' > "./ConversionRate/"$i"_conversion_rate.txt"

        # conversion rate pete - not sure if this is correct...
        # using plastid CHH unconverted C rate: 100-(sum(#C_counts)/sum(#CT_counts)*100)
        # Note that bsmap recomends using the eff_CT_counts (field #$7) which considers if there is a mismatch with the reverse strand.
        # However, Qing recommends just using the CT count (field #$9) becasue the reverse strand could equally be a sequencing error. Check this.
        awk -F$"\\t" \
        'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
        BSMAPratio/${ID}_BSMAP_out.txt | \
        awk '{sum1 += $9; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' > ConversionRate/${ID}_conversion_rate.txt
        # conversion rate pete - using eff_CT
        awk -F$"\\t" \
        'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
        BSMAPratio/${ID}_BSMAP_out.txt | \
        awk '{sum1 += $7; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}' > ConversionRate/${ID}_conversion_rate_eff_C.txt
        #debugging
        #awk -F$"\\t" \
        #'BEGIN {OFS = FS} {if($1=="Pt" && $5=="CHH") print}' \
        #BSMAPratio/F1-16_Index5_S1_BSMAP_out.txt | \
        #awk -F$"\\t" 'BEGIN {OFS = FS} {sum1 += $7; sum2 +=$8} END {print sum1, sum2 , 100-((sum2/sum1)*100)}'

        ########################
        #Target region analysis
        # count methylated and unmethylated Cytosine per specific target region

        ## count read depth over each specific target region, only use uniquley mapped and properly paried reads
        ## args to be passed from qsub script ${intersect_regions_ref} eg ${refdir}/BSseqcapv2_specific_regions.bed

        bedtools intersect \
        -abam bsmapped_filtered/${ID}_sorted_MarkDup_pairs_clipOverlap.bam \
        -b ${intersect_regions_ref} \
        -bed \
        -wa \
        -wb > TempOut/${ID}_specific_region_pairs_clipOverlap.txt

        awk -F$"\\t" 'BEGIN {OFS = FS} {sum[$16]++} END {for (i in sum) print i, sum[i]}' \
        TempOut/${ID}_specific_region_pairs_clipOverlap.txt | \
        sort -k1 -n > \
        OnTargetCoverage/${ID}_specific_region_count_pairs_clipOverlap.txt
        #debugging
        #awk -F$"\\t" 'BEGIN {OFS = FS} {sum[$16]++} END {for (i in sum) print i, sum[i]}' \
        #TempOut/F1-16_Index5_S1_specific_region_pairs_clipOverlap.txt | \
        #sort -k1 -n -t | head

        # subset the methylations counts file to only those sites that overlap target regions
        bedtools intersect \
        -a BSMAPratio/${ID}_BSMAP_out.txt \
        -b ${intersect_regions_ref} \
        -wa \
        -wb > TempOut/${ID}_BSMAP_out_ontarget.txt

        # use awk to total CTs in each region and count mC in each context to calculate ratios
        # *** Think about this method: it will be swayed by high count sites - should we average the ratios at each site?
        # Note that bsmap recomends using the eff_CT_counts (field #$7) which considers if there is a mismatch with the reverse strand. Check this.
        # However, Qing recommends just using the CT count (field #$9) becasue the reverse strand could equally be a sequencing error.
        awk \
        -F$"\\t" 'BEGIN {OFS = FS}
        {mC[$14"\\t"$15"\\t"$16"\\t"$5] += $8;
        CT[$14"\\t"$15"\\t"$16"\\t"$5] += $9;
        n[$14"\\t"$15"\\t"$16"\\t"$5]++} END {
        for (j in mC) print j, n[j], mC[j], CT[j]}' \
        TempOut/${ID}_BSMAP_out_ontarget.txt > OnTargetCoverage/${ID}_BSMAP_out_ontarget_mC.txt

        # count Cs again, this time use eff_CT_counts
        awk \
        -F$"\\t" 'BEGIN {OFS = FS}
        {mC[$14"\\t"$15"\\t"$16"\\t"$5] += $8;
        CT[$14"\\t"$15"\\t"$16"\\t"$5] += $7;
        n[$14"\\t"$15"\\t"$16"\\t"$5]++} END {
        for (j in mC) print j, n[j], mC[j], CT[j]}' \
        TempOut/${ID}_BSMAP_out_ontarget.txt > OnTargetCoverage/${ID}_BSMAP_out_ontarget_mCeff.txt

        # awk \
        # -F$"\\t" 'BEGIN {OFS = FS}
        # {mC[$14"\t"$15"\t"$16"\t"$5] += $8;
        # CT[$14"\t"$15"\t"$16"\t"$5] += $9;
        # n[$14"\t"$15"\t"$16"\t"$5]++} END {
        # for (j in mC) print j, n[j], mC[j], CT[j]}' \
        # TempOut/F1-16_Index5_S1_BSMAP_out_ontarget.txt |head

        ########################

        # 100 bp tiles using Qing's perl script
        #cd BSMAPratio
        #perl ~/gitrepos/springerlab_methylation/SeqCap/met_context_window.pl ${ID}_BSMAP_out.txt 100

        #Make bigWigs
        #bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CG.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        #"BSMAPratio/${ID}_BSMAP_out.txt.100.CG.bigWig"
        #bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CHG.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        #"BSMAPratio/${ID}_BSMAP_out.txt.100.CHG.bigWig"
        #bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CHH.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        #"BSMAPratio/${ID}_BSMAP_out.txt.100.CHH.bigWig"

echo finished summarising
