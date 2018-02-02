#!/bin/bash -l
#PBS -l walltime=8:00:00,nodes=1:ppn=1,mem=50gb
#PBS -N tiles_bigWigs
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

module load R/3.3.2
########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis
mkdir -p tiles

########## Run #################

        # 100 bp tiles using Qing's perl script
        # the output file format is actually a bed it is a zero-based coordinate
        perl ~/gitrepos/springerlab_methylation/SeqCap/met_context_window.pl ./BSMAPratio/${ID}_BSMAP_out.txt 100

        #move output into tiles folder
        mv -v ./BSMAPratio/${ID}_BSMAP_out.txt.100.*.bed ./tiles/

        #Run R moudle to:
        # 1. fix chr ends in BED file
        # 2. output one-based coordinate .txt file for analysis in R
        # 3. output zero-based coordinate bedGraph file (omit excess columns and sort)
        R -f ~/gitrepos/springerlab_methylation/SeqCap/07-tiles_bed_to_bigWig.R \
        --args ${ID} tiles ${reference_tile_file}

        #make bedGraph by sorting and removing cols 4 and 5 with awk

        # make bedGraph by removing cols 4 and 5 and sorting
        #cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.sorted.bg
        #cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.sorted.bg
        #cut -f1-3,6-6 ./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.bed | sort -k1,1 -k2,2n > ./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.sorted.bg

        #Make bigWigs
        # At some point soft code the reference, make variable in script call
        # To mkae the chrom_sizes file use seqtk compressseqtk comp Zea_mays.AGPv4.dna.toplevel.fa
        # seqtk comp Zea_mays.AGPv4.dna.toplevel.fa | cut -f 1-2 > Zea_mays.AGPv4.dna.toplevel.chrom.sizes

        bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CG.fixed.sorted.bg" ${chrom_sizes} \
        "./tiles/${ID}_BSMAP_out.txt.100.CG.bigWig"
        bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CHG.fixed.sorted.bg" ${chrom_sizes} \
        "./tiles/${ID}_BSMAP_out.txt.100.CHG.bigWig"
        bedGraphToBigWig "./tiles/${ID}_BSMAP_out.txt.100.CHH.fixed.sorted.bg" ${chrom_sizes} \
        "./tiles/${ID}_BSMAP_out.txt.100.CHH.bigWig"

        # remove .bg - not required
        # meh purge it in step 09 instead
        # rm -rv tiles/${ID}*.bg

echo finished summarising
