#! /usr/bin/perl -w
# bedtools closest between atacseq fpkm and genes... overlay with RNAseq data
# 27_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I))|| (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "atac_chr\tatac_start\tatac_end\tatac_fpkm\tgene_chr\tensembl\ttype\tgene_start\tgene_end\tdot\tstarnd\tdot\tgene_id\tgene_atac_dist\ttranscript\t207_Leaf_rep1\t207_Leaf_rep2\t207_Root_rep1\t207_Root_rep2\tB73_Leaf_rep1\tB73_Leaf_rep2\tB73_Root_rep1\tB73_Root_rep2\tPHJ89_Leaf_rep1\tW22_Leaf_rep1\tW22_Leaf_rep2\tW22_Root_rep1\n";

my %gene;
while (my $line = <$in_fh>) {
  chomp $line;

  # 1	1	100	N/A	1	ENSEMBLPEP	gene	44289	49837	.	+	.	ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene

  my ($chr, $start, $end, $fpkm, $gene_chr, undef, undef, $gene_start, $gene_end, undef, $strand, undef, $id, $dist) = split ("\t", $line);

  $id =~ s/ID=//;
  my ($gene_id, undef) = split(";", $id);

  $gene{$gene_id} = $line;
}

my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  # transcript	207_Leaf_rep1	207_Leaf_rep2	207_Root_rep1	207_Root_rep2	B73_Leaf_rep1	B73_Leaf_rep2	B73_Root_rep1	B73_Root_rep2	PHJ89_Leaf_rep1	W22_Leaf_rep1	W22_Leaf_rep2	W22_Root_rep1
  # gene:Zm00001d001763	0.156249291995396	0.721215240468036	6.13192190117661	6.78783344935149	0.691360612239848	0.584624930952423	0.129278688266033	0.576504300595972	10.016420397003	5.61185924782014	0.715600724775479	0.663561776044532

  my ($gene_id, undef) = split ("\t", $line2);

  if (exists ($gene{$gene_id})) {
    print $out_fh "$gene{$gene_id}\t$line2\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
