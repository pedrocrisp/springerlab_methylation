#!/bin/bash -l
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=50gb
#PBS -N 05.2-analysis-contextMeans
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
echo data folder is $data_folder
echo coverage filter is $coverage_filter

mkdir -p analysis/contextMeans

########## Run #################

########################
#Make bedGraphs per context including CT

# begGraph ratio files for tdfs
awk_make_bedGraph='BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $8/$9*100, $5, $9
}
'

# split bedGraph by contex
awk_make_bedGraph_context='BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $4, $6 > "analysis/BSMAPratio/"ID"_BSMAP_out_"$5".bedGraph"
}
'

# split bedGraph by sub-contex
# only difference in the output file suffix (so that we make CG files for context and sub-contex: sainty check - they should be the same)
awk_make_bedGraph_subcontext='BEGIN {OFS = FS} (NR>1){
  print $1, $2, $3, $4, $6 > "analysis/BSMAPratio/"ID"_BSMAP_out_subcontext_"$5".bedGraph"
}
'

#pipe bedGraph to split by context (use dash to read from sdtin)
# per context
awk -F$"\\t" "$awk_make_bedGraph" \
"analysis/BSMAPratio/${ID}_BSMAP_out.txt" | \
awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_context" -

# per sub-context
awk -F$"\\t" "$awk_make_bedGraph" \
"analysis/BSMAPratio/${ID}_BSMAP_out_subcontext.txt" | \
awk -F$"\\t" -v ID=$ID "$awk_make_bedGraph_subcontext" -

#########################
        #Run R moudle to:
        # 1. make CHH coverage files for summary stats analysis
        R -f ~/gitrepos/springerlab_methylation/SeqCap/05.2-analysis-contextMeans-WGBS.R \
        --args ${ID} $data_folder $coverage_filter $loci_of_interst_file $loci_of_interest_name

#########################
# remove unnecesay bedGraphs

rm -rv analysis/BSMAPratio/${ID}*.bedGraph

echo finished summarising
