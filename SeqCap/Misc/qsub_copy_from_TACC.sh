#THIS DOESNT WORK BECAUSE OF 2-STEP LOGIN ON TACC...
#USE THE SCREEN-INTERACTIVE METHOD ON LAB SERV

#!/bin/bash -l
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=4gb
#-o copy_TACC_o
#-e copy_TACC_e
#PBS -N copy_TACC
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

# cd to dir
cd ~/ws/rawdata/SeqCap_from_TACC
# start copy
rsync -rhivPt pcrisp@tacc:/corral-tacc/projects/iplant/vaughn/springer_vaughn/eichten/seqcap_jan_2017 ./
