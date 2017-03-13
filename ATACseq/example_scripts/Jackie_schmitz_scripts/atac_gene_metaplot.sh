#Jan_11_2017
#atac_gene_metaplot

#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atac_gene_metaplot_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atac_gene_metaplot_e
#PBS -N atac_gene_metaplot
#PBS -V
#PBS -r n

ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load bedtools/2.25.0

cd /scratch.global/nosha003/schmitz

## Bins - all genes
#closest gene to bin --> distances based on from bin (need to flip sign to make dist from gene)
#sort -k1,1n -k2,2n atacseq/counts/${ID}.100bpcounts.bed | sed '1d' > atacseq/counts/${ID}.100bpcounts_sort.bed
#bedtools closest -D a -t first -a atacseq/counts/${ID}.100bpcounts_sort.bed -b /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 > atacseq/metaplot/${ID}_gene.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data.pl -i atacseq/metaplot/${ID}_gene.txt -o atacseq/metaplot/${ID}_gene_metaplot.txt


## Bins - tissue specific genes
#root specific B73
#perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i rnaseq/rpm_matrix_rootspecific_B73.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o rnaseq/rootspecific_B73_gene.gff3
#closest gene to bin
#bedtools closest -D a -t first -a atacseq/counts/${ID}.100bpcounts_sort.bed -b rnaseq/rootspecific_B73_gene.gff3 > atacseq/metaplot/${ID}_gene_rootB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data.pl -i atacseq/metaplot/${ID}_gene_rootB73.txt -o atacseq/metaplot/${ID}_gene_rootB73_metaplot.txt

#leaf specific B73
#perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i rnaseq/rpm_matrix_leafspecific_B73.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o rnaseq/leafspecific_B73_gene.gff3
#closest gene to bin
#bedtools closest -D a -t first -a atacseq/counts/${ID}.100bpcounts_sort.bed -b rnaseq/leafspecific_B73_gene.gff3 > atacseq/metaplot/${ID}_gene_leafB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data.pl -i atacseq/metaplot/${ID}_gene_leafB73.txt -o atacseq/metaplot/${ID}_gene_leafB73_metaplot.txt


## Bins - expression specific genes
#active B73leaf genes
#perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i /scratch.global/nosha003/schmitz/rnaseq/rpm_matrix_active_B73leaf.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o /scratch.global/nosha003/schmitz/rnaseq/active_B73leaf_gene.gff3
#closest gene to bin
#bedtools closest -D a -t first -a atacseq/counts/${ID}.100bpcounts_sort.bed -b rnaseq/active_B73leaf_gene.gff3 > atacseq/metaplot/${ID}_gene_activeB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data.pl -i atacseq/metaplot/${ID}_gene_activeB73.txt -o atacseq/metaplot/${ID}_gene_activeB73_metaplot.txt

#active B73leaf genes
#perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i /scratch.global/nosha003/schmitz/rnaseq/rpm_matrix_silent_B73leaf.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o /scratch.global/nosha003/schmitz/rnaseq/silent_B73leaf_gene.gff3
#closest gene to bin
#bedtools closest -D a -t first -a atacseq/counts/${ID}.100bpcounts_sort.bed -b rnaseq/silent_B73leaf_gene.gff3 > atacseq/metaplot/${ID}_gene_silentB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data.pl -i atacseq/metaplot/${ID}_gene_silentB73.txt -o atacseq/metaplot/${ID}_gene_silentB73_metaplot.txt


#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atac_gene_metaplot.sh