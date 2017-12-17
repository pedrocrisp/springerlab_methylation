#!/bin/bash
#Set -e as an option, it tells the command line to exit the script immediately if it registers an error.
set -e

#Set -x as an option, it tells the computer to echo back each step before it is completed and the output is produced.
set -x

###
#code to make script work on both osx and linux. Essentially, it creates a file path to the script directory and saves this path as $0. In detail: 'if the operating system type is darwin (a mac), then use the greadlink function when the readlink function is called. Then use the greadlink function to find the script file named. In doing so, find the path of the script files directory and save as 'scriptdir'. This change is made for Macs because readlink doesn't run properly, but greadlink does. If the OS is not mac (eg. Linux), find the script file using the readlink function and save the path to the script file directory as 'scriptdir.' By using readlink to find where the scripts are, it means if this pipeline is copied onto another computer, the files can still be found.
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

#new suffix for output file
outputFile="${fqname%%.*}.fastq.gz"

# cutadapt remove 20 bp 5' end
cutadapt -u -20 -o trimmed/$outputFile fq

# to run
# find . -name "*R2_001.fastq.gz" | parallel -j 12 bash ~/gitrepos/springerlab_methylation/cutadapt-5ptime.sh {}
