#! /usr/bin/perl -w
# take metaplot created file for all genes and subset only for specific genes
# atac_metaplot_subsetgenes.pl
# 17_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in files 
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my $header = <$in_fh>;
my $header2 = <$in2_fh>;
print $out_fh "chr\tbinstart\tbinstop\tbinid\tchr2\tgene\tgenestart\tgenestop\tgenesize\tstrand\tdistance\tstrand_distance\trelative_distance\treal_distance\tcount\tlog2fc\n";

my %subset;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($gene, $basemean, $log2fc, $lfcse, $stat, $pval, $padj) = split ("\t", $line);

# "gene:Zm00001d001766"	245.91193709757	-1.41978090277316	0.772298793483918	-1.83838290924732	0.0660060017226724	0.187474082556142

  $gene =~ s/"//g;
  $gene =~ s/gene://;

  $subset{$gene} = $log2fc;
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, $binstart, $binstop, $binid, $chr2, $gene, $genestart, $genestop, $genesize, $strand, $distance, $strand_distance, $relative_distance, $real_distance, $count) = split ("\t", $line2);

# 1	101	200	ID=bin1-2	1	Zm00001d027230	44289	49837	5548	+	-44089	-44089	-44089	-44089	78

  if (exists ($subset{$gene})) {
    print $out_fh "$line2\t$subset{$gene}\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;

