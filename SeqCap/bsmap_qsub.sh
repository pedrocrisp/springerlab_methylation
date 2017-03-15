#!/bin/bash
#set -x
set -e

usage="USAGE:
seqcap_pipeline_test_v0.1_qsub.sh <sample_list.txt>"

######### Setup ################
sample_list=$1
if [ "$#" -lt "1" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for methylome analysis"
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t=1
else
qsub_t="1-${number_of_samples}"
fi

echo "argument to be passed to qsub -t is '$qsub_t'"
########## Run #################



#make folder structure
timestamp=$(date +%Y%m%d-%H%M%S)

#get workingdir - not working try to fix later when I make script generic
#pbs_pwd=$(pwd)
#echo $pbs_pwd

workingdir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei

analysis_dir=/scratch.global/pcrisp/SeqCap_1_Mei/testpipeline_SeqCap_1_Mei/analysis

log_folder=${analysis_dir}/logs_bsmap_${timestamp}
mkdir $log_folder

trimmed=${analysis_dir}/trimmed

alignfolder=${analysis_dir}/bsmaped
mkdir $alignfolder

#script path and cat a record of what was run
script_dir=/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap
script_to_qsub=${script_dir}/bsmap.sh
cat $script_to_qsub > "$log_folder/script.log"
cat $0 > "$log_folder/runner.log"

#run script and pass args
#bash $script $workingdir $log_folder $fastqcfolder $trimmed $alignfolder $BSMAPratio_folder $TempOut $OnTargetCoverage
qsub -t $qsub_t \
-o ${log_folder}/testpipeline_SeqCap_1_Mei_o \
-e ${log_folder}/testpipeline_SeqCap_1_Mei_e \
-v LIST=${workingdir}/${sample_list},readsdir=${workingdir}/reads,workingdir=$analysis_dir,log_folder=$log_folder,trimmedfolder=$trimmed,alignfolder=$alignfolder \
$script_to_qsub

    #-wd $pbs_pwd \
    #-wd Execute the job from the directory specified in working_dir. - not working...
    #-o and -e pass the file locations for std out/error
    #-v additional variables to pass to the qsub script including the PBS_array list and the dir structures

#to run
#bash /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/bsmap_qsub.sh <sample_list.txt>
