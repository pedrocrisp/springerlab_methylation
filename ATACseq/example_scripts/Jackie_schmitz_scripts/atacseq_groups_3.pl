#! /usr/bin/perl -w
# subset data and create new column to define groups
# atacseq_groups.pl
# 10_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out\n";
our ($opt_i, $opt_o, $opt_h);
getopts("i:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation_atacseq_nonzero data files
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tbinstart\tbinend\tCG_ratio\tCHG_ratio\tCHH_ratio\tatacseq_count\tgroup\n";

my $group;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $binstart, $binend, $CG, $CHG, $CHH, $atac) = split ("\t", $line);

  if ($atac < 20 && $atac > 0) {
    $group = "1";
  }

  if ($atac < 30 && $atac > 20) {
    $group = "2";
  }

  if ($atac > 30) {
    $group = "3";
  } 

  print $out_fh "$chr\t$binstart\t$binend\t$CG\t$CHG\t$CHH\t$atac\t$group\n";
}

close $in_fh;
close $out_fh;
exit;
