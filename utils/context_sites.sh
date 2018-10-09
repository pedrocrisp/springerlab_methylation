#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced.
set -x

###
#code to make script work on both osx and linux.

if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

#The fasta file
fa=$1

#basename
faname="$(basename $fa)"

echo $faname

#CG
outputFile="${faname%%.*}_CG.txt"
cat $fa | seqkit locate -i -p "CG" | cut -f 1-2,4-6 > sites_files/$outputFile

#CHG
outputFile="${faname%%.*}_CHG.txt"
cat $fa | seqkit locate -i -p "C[ATC]G" | cut -f 1-2,4-6 > sites_files/$outputFile

#CG
outputFile="${faname%%.*}_CHH.txt"
cat $fa | seqkit locate -i -p "C[ACT][ACT]" | cut -f 1-2,4-6 > sites_files/$outputFile

# to run
# find . -name "*.fastq." | parallel -j 20 bash ~/gitrepos/springerlab_methylation/utils/context_sites.sh {} >> $(date +%Y%m%d-%H%M%S)_make_sites.log 2>&1
