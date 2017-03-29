#bin/bash -l
#PBS -l walltime=01:00:00,nodes=1:ppn=8,mem=16gb
#PBS -N extra_awk_ratio
#PBS -r n
#PBS -m abe
#PBS -M pcrisp@umn.edu

set -ex

# this script work but only if there are no zeros in the total CT column (which happens...); therefore this doent work
for ID in $(cat last_sample.txt);
  do
    cd analysis
    #count Cs and calculate methylation ratio
        awk \
        -F$"\\t" 'BEGIN {OFS = FS}
        {mC[$14"\\t"$15"\\t"$16"\\t"$5] += $8;
        CT[$14"\\t"$15"\\t"$16"\\t"$5] += $9;
        n[$14"\\t"$15"\\t"$16"\\t"$5]++} END {
        for (j in mC) print j, n[j], mC[j], CT[j], mC[j]/CT[j]}' \
        TempOut/${ID}_BSMAP_out_ontarget.txt > OnTargetCoverage/${ID}_BSMAP_out_ontarget_mC_ratio.txt

        #count Cs again, this time use eff_CT_counts
        awk \
        -F$"\\t" 'BEGIN {OFS = FS}
        {mC[$14"\\t"$15"\\t"$16"\\t"$5] += $8;
        CT[$14"\\t"$15"\\t"$16"\\t"$5] += $7;
        n[$14"\\t"$15"\\t"$16"\\t"$5]++} END {
        for (j in mC) print j, n[j], mC[j], CT[j], mC[j]/CT[j]}' \
        TempOut/${ID}_BSMAP_out_ontarget.txt > OnTargetCoverage/${ID}_BSMAP_out_ontarget_mCeff_ratio.txt
        cd -
      done
