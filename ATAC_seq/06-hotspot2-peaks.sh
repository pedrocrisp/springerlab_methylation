#!/bin/bash -l
#PBS -l walltime=24:00:00,nodes=1:ppn=6,mem=40gb
#PBS -N hotspot2-peaks
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

PATH=$PATH:~/bin/bedops_v2.4.35
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID
#sample_dir="${reads_folder}/${ID}"

#make out folder
cd analysis
outdir=${bam_folder}_hotspot_peaks
mkdir -p $outdir

########## Run #################

# run on default setting
# requires a centers file generated with the script ~/software/hotspot2-2.1.1/scripts/extractCenterSites.sh
# note that the chrom_sizes files is a bed file and requires start column
# to see usage: bash ~/software/hotspot2-2.1.1/scripts/hotspot2.sh -h

bash ~/software/hotspot2-2.1.1/scripts/hotspot2.sh \
-c $chrom_sizes \
-C $centers_starch \
${bam_folder}/${ID}*.bam \
$outdir
