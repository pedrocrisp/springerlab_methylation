#!/bin/bash -l
#PBS -l walltime=48:00:00,nodes=1:ppn=8,mem=80gb
#PBS -N bowtie2_batch
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

module load bowtie2/2.3.4.1
module load samtools/1.9

########## Set up dirs #################

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID
#sample_dir="${reads_folder}/${ID}"

# I just changed this bit to deal with nested and non-nested folder structures - it may break old analyses?...
fastqs="$(find $reads_folder -type f -name ${ID}*.fq*)"
# convert to array to count elements
fastqs_count=($fastqs)

#make adaligned folder bowtie2 (caution this will not fail if dir already exists)
outdir="align_bowtie2"
mkdir -p ${outdir}

# output structure
outsam="${outdir}/${ID}.sam"
tmpbam="${outdir}/${ID}.bam"
outbam="${outdir}/${ID}_sorted.bam"

########## Run #################

#-x bowtie index
#--phred33
#--end-to-end dont trim reads to enable alignment
#--mm memory-mapped I/O to load index allow multiple processes to use index
#-p number of threds to use
#-U fastq file path, and specifies reads are not paired
#-S file to write SAM alignemnts too (although default is stdout anyhow)
#-D and -R tell bowtie to try a little harder than normal to find alignments
#-L reduce substring length to 10 (default 22) as these are short reads
#-i reduce substring interval? more sensitive?
#-N max # mismatches in seed alignment; can be 0 or 1 (0)
#-D give up extending after <int> failed extends in a row (15)
# -k report N mapping locations
# --score-min L,0,0 sets formular for min aln score for alignment to be reported, eg readlength 36 * -0.6 + -0.6) = min 0; this means only report exact matches
# Bowtie 2 does not "find" alignments in any specific order,
# so for reads that have more than N distinct, valid alignments,
# Bowtie 2 does not guarantee that the N alignments reported are the best possible in terms of alignment score.
# Still, this mode can be effective and fast in situations where the user cares more about whether a read aligns
# (or aligns a certain number of times) than where exactly it originated.

if (( "${#fastqs_count[@]}" == 2 )); then

echo "paired reads"

# this is pretty non generic but I need a quick fix...
# this will only work if the suffix is consistent with trimmgalore output
fq_1="${reads_folder}/${ID}_R1_001_val_1.fq"
fq_2="${reads_folder}/${ID}_R2_001_val_2.fq"

bowtie2 \
-x $bt2_genome \
--phred33 \
--end-to-end \
--mm \
-k $multimapping_rate \
-D 20 \
-R 3 \
-N 0 \
-L 10 \
-i S,1,0.50 \
-p $bt2_threads \
-1 $fq_1 \
-2 $fq_2 \
-S "$outsam"

else
echo "assuming single end"

bowtie2 \
-x $bt2_genome \
--phred33 \
--end-to-end \
--mm \
-k $multimapping_rate \
-D 20 \
-R 3 \
-N 0 \
-L 10 \
-i S,1,0.50 \
-p $bt2_threads \
-U $fastqs \
-S "$outsam"

fi


###### sort and index
# filtering: -q 10 is MAPQ >= 10; 1 in 10 chance mapping location is wrong for example
samtools view -q $MAPQ_threshold -b -@ 8 $outsam | samtools sort -m 8 -@ 8 -o $outbam

#Make an index of the sorted bam file
samtools index ${outbam}

#Delete the temporary sam.
rm -v ${outsam}
