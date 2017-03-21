

# extract methylation information using bsmap tool methratio.py
        python /home/springer/pcrisp/software/bsmap-2.90/methratio.py \
        -o BSMAPratio/${ID}_methratio.txt \
        -d ${genome_reference} \
        -u \
        -z \
        -r bsmapped/${ID}_sorted_MarkDup_pairs_clipOverlap.bam
        
        #awk funciton for extracting methylation info from methratio.py output. Check with Qing what this is meant to do. Also try to figure out how to split this over multiple lines
        #awk '(NR>1){if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2-1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $1"\t"$2-1"\t"$2"\t"$3"\t""CNN""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' BSMAPratio/${ID}.txt > BSMAPratio/${ID}_BSMAP_out.txt
        
        cd BSMAPratio
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
        #awk -F $'\t' "$awk_make_bed" F1-16_Index5_S1_methratio.txt > F1-16_Index5_S1.bed
         
        awk -F $'\t' "$awk_make_bed" ${ID}_methratio.txt > ${ID}_BSMAP_out.txt
        
        cd ..
        
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
