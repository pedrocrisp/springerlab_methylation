#!/bin/bash -l
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=20gb
#PBS -N site_cov_tiles
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
output_folder=${data_folder}_filtered_${filter_suffix}
mkdir -p $output_folder

#get job ID
#use sed, -n supression pattern space, then 'p' to print item number {PBS_ARRAYID} eg 2 from {list}
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

echo sample being mapped is $ID
echo data folder is $data_folder
echo tile list is $tile_list

########## Run #################

        #Run R moudle to:
        # 1. make CHH coverage files for summary stats analysis
        R -f ~/gitrepos/springerlab_methylation/SeqCap/08b-alternate-tiles_analysis_filter_list_CT_GA.R \
        --args ${ID} $data_folder $tile_list $output_folder $minCHHs $minCHH_cov $minCHGs $minCHG_cov $minCGs $minCG_cov $minCHH_A_percent $minCHH_GA_cov $minCHG_A_percent $minCHG_GA_cov $minCG_A_percent $minCG_GA_cov

echo finished summarising
