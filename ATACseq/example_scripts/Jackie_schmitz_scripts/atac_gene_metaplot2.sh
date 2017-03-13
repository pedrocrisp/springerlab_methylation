#Feb_27_2017
#atac_gene_metaplot2

#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atac_gene_metaplot_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atac_gene_metaplot_e
#PBS -N atac_gene_metaplot
#PBS -V
#PBS -r n

#load modules
module load bedtools/2.25.0

## Bins - all genes
#closest gene to bin --> distances based on from bin (need to flip sign to make dist from gene)
sed '1d' /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_repcombine.txt | sort -k 1,1n -k 2,2n > /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_repcombine_noheader.txt
bedtools closest -D a -t first -a /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_repcombine_noheader.txt -b /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 > /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data_fpkm.pl -i /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene.txt -o /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_metaplot.txt

## Bins - expression specific genes
#active B73leaf genes
perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i /scratch.global/nosha003/schmitz/rnaseq/rpm_matrix_active_B73leaf.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o /scratch.global/nosha003/schmitz/rnaseq/active_B73leaf_gene.gff3
#closest gene to bin
bedtools closest -D a -t first -a /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_repcombine_noheader.txt -b /scratch.global/nosha003/schmitz/rnaseq/active_B73leaf_gene.gff3 > /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_activeB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data_fpkm.pl -i /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_activeB73.txt -o /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_activeB73_metaplot.txt

#active B73leaf genes
perl /home/springer/nosha003/schmitz/scripts/gene_subset_gff.pl -i /scratch.global/nosha003/schmitz/rnaseq/rpm_matrix_silent_B73leaf.txt -I /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 -o /scratch.global/nosha003/schmitz/rnaseq/silent_B73leaf_gene.gff3
#closest gene to bin
bedtools closest -D a -t first -a /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_repcombine_noheader.txt -b /scratch.global/nosha003/schmitz/rnaseq/silent_B73leaf_gene.gff3 > /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_silentB73.txt
#make file for metaplot
perl /home/springer/nosha003/wenli_data/scripts/metaplot_data_fpkm.pl -i /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_silentB73.txt -o /scratch.global/nosha003/schmitz/atacseq/metaplot/atac_fpkm_gene_silentB73_metaplot.txt


#qsub /home/springer/nosha003/schmitz/scripts/atac_gene_metaplot2.sh