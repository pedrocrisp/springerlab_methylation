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
  my $log2_atac = log($atac) / log(2);

  if ($log2_atac <= -3) {
    $group = "1";
  }

  if ($log2_atac > -3 && $log2_atac < -1) {
    $group = "2";
  }

  if ($log2_atac >= -1 && $log2_atac <= 1) {
    $group = "3";
  }

  if ($log2_atac > 1 && $log2_atac < 3 ) {
    $group = "4";
  }

  if ($log2_atac >= 3) {
    $group = "5";
  }
  
  print $out_fh "$chr\t$binstart\t$binend\t$CG\t$CHG\t$CHH\t$atac\t$group\n";
}

close $in_fh;
close $out_fh;
exit;
