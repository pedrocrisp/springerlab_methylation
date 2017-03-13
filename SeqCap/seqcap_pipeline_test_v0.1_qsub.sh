#!/bin/bash

#make folder structure
timestamp=$(date +%Y%m%d-%H%M%S)

mkdir analysis

workingdir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/analysis

log_folder=${workingdir}/logs_${timestamp}
mkdir $log_folder

fastqcfolder=${workingdir}/fastqc
mkdir $fastqcfolder

trimmed=${workingdir}/trimmed
mkdir $trimmed

alignfolder=${workingdir}/bsmaped
mkdir $alignfolder

BSMAPratio_folder=${workingdir}/BSMAPratio
mkdir $BSMAPratio_folder

TempOut=${workingdir}/TempOut
mkdir $TempOut

OnTargetCoverage=${workingdir}/OnTargetCoverage
mkdir $OnTargetCoverage

#script path and cat a record of what was run
script_dir=/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap
script_to_qsub=${script_dir}/seqcap_pipeline_test_v0.1.sh
cat $script_to_qsub > "$log_folder/script.log"
cat $0 > "$log_folder/runner.log"

#run script and pass args
#bash $script $workingdir $log_folder $fastqcfolder $trimmed $alignfolder $BSMAPratio_folder $TempOut $OnTargetCoverage
qsub -t 1-8 -v LIST=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/samples.txt $script_to_qsub

#to run
#bash /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/seqcap_pipeline_test_v0.1_qsub.sh
