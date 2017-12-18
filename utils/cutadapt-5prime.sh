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

module load cutadapt/1.8.1

#cutadapt
#note the program has changed since it publication
#https://cutadapt.readthedocs.org/en/stable/guide.html#cutadapt-s-output


#The fastq file
fq=$1

#basename
fqname="$(basename $fq)"

echo $fqname

#new suffix for output file (using the same, but modify if desired)
outputFile="${fqname%%.*}.fastq"

# cutadapt remove 20 bp 5' end
cutadapt -u 20 -o trimmed/$outputFile $fq

gzip trimmed/$outputFile

# to run
# find . -name "*R2_001.fastq.gz" | parallel -j 12 bash ~/gitrepos/springerlab_methylation/utils/cutadapt-5ptime.sh {}
