#pcrisp
#2017-04-04

#!/bin/bash -l
#PBS -l walltime=00:10:00,nodes=1:ppn=1,mem=4gb
#-o copy_logs_o
#-e copy_logs_e
#PBS -N copy_s3
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

#to run
qsub <~/gitrepos/springerlab_methylation/SeqCap/Misc/qsub_copy_to_s3.sh> <-v destination_s3_location=target_s3_location>

destination_s3_location=$1

# start copy
s3cmd \
sync \
-p \
./ \
$destination_s3_location
