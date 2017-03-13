#! /usr/bin/perl -w
# subset gene gff file to only include genes with tissue specific expression
# gene_subset_gff.pl
# 11_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation_atacseq_nonzero data files
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my %genes;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($gene, undef) = split ("\t", $line);

  $genes{$gene} = 'NA';
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, undef, undef, $start, $end, undef, $strand, undef, $id) = split ("\t", $line2);

# 1	ENSEMBLPEP	gene	44289	49837	.	+	.	ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene

  my ($gene, undef) = split (";", $id);
  $gene =~ s/ID=//;
  #print "$gene\t";

  if (exists ($genes{$gene})) {
    print $out_fh "$line2\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
