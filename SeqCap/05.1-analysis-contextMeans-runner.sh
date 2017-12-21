#!/bin/bash
#set -xe
set -xeuo pipefail

#Peter Crisp
#2017-04-08
#Bash parallel runner script for 05.1-analysis-contextMeans.R
#To be run in an interactive session
#eg qsub -I -l walltime=24:00:00,nodes=1:ppn=8,mem=40gb

###
#code to make it work on osx and linux
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

usage="USAGE:
05.1-analysis-contextMeans-runner.sh <sample_list> <threads>
eg bash ~/gitrepos/springerlab_methylation/SeqCap/05.1-analysis-contextMeans-runner.sh samples.txt 2"

######### Setup ################
sample_list=$1
threads=$2
#minLength=$3
#trimLength=$4
#minCoverage=$5

if [ "$#" -lt "2" ]
then
  echo $usage
  exit -1
else
  cat $sample_list
  echo Summarising methylation data Iniating $2 parallel jobs
fi
########## Run #################
module load parallel

#user defined variables that could be changed:
workingdir=./
script=$scriptdir/05.1-analysis-contextMeans.R
###

# function findSamples () {
#     find $BSMAPratio_folder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
# }

outdir="analysis_context"
mkdir -p ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="${outdir}/.logs_${timestamp}"
mkdir $logdir

cat $script > "${logdir}/script.log"
cat $0 > "${logdir}/runner.log"
cat $script

cat $sample_list | parallel -j $threads R -f $script --args {} $outdir \>${logdir}/{}.log 2\>\&1
