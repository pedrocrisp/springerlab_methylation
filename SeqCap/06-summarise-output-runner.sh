#!/bin/bash
#set -xe
set -xeuo pipefail

#Peter Crisp
#2017-04-08
#Bash parallel runner script for 06-summarise-output.R

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
06-summarise-output-runner.sh <sample_list> <OnTargetCoverage_folder> <threads>
eg bash ~/gitrepos/springerlab_methylation/SeqCap/06-summarise-output-runner.sh samples.txt OnTargetCoverage 6"

######### Setup ################
sample_list=$1
OnTargetCoverage_folder=$2
threads=$3
#minLength=$3
#trimLength=$4
#minCoverage=$5

if [ "$#" -lt "3" ]
then
  echo $usage
  exit -1
else
  cat $sample_list
  echo Summarising methylation data using folder = $1\n Iniating $2 parallel jobs
fi
########## Run #################


#user defined variables that could be changed:
workingdir=./
script=$scriptdir/06-summarise-output.R
###

# function findSamples () {
#     find $BSMAPratio_folder/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
# }

outdir="${OnTargetCoverage_folder}_annotated"
mkdir ${outdir}
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="${outdir}/.logs_${timestamp}"
mkdir $logdir

cat $script > "${logdir}/script.log"
cat $0 > "${logdir}/runner.log"
cat $script

cat $sample_list | parallel -j $threads R -f $script --args {} $OnTargetCoverage_folder $outdir \>${logdir}/{}.log 2\>\&1
