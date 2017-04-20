#!/bin/bash -l
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=32gb
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


########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID

#make adaligned folder bsmaped
cd analysis

########## Run #################

        cd BSMAPratio

        #Make bigWigs
        bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CG.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out.txt.100.CG.bigWig"
        bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CHG.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out.txt.100.CHG.bigWig"
        bedGraphToBigWig "${ID}_BSMAP_out.txt.100.CHH.bed" ~/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.chrom.sizes \
        "BSMAPratio/${ID}_BSMAP_out.txt.100.CHH.bigWig"

echo finished summarising
