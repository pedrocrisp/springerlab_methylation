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
