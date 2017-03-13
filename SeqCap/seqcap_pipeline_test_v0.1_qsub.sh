#!/bin/bash

#script path and cat a record of what was run
script_dir=/home/springer/pcrisp/gitrepos/springerlab_methylation/
qsub_script=${script_dir}/seqcap_pipeline_test_v0.1.sh
cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"

timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

#make log folder
workingdir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei
log_folder=${workingdir}/logs
fastqcfolder=${workingdir}/fastqc

#run script and pass args
bash $script $workingdir $log_folder $fastqcfolder
