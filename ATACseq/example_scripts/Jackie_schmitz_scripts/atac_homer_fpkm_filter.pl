#! /usr/bin/perl -w
# script to normalize homer peak counts (# reads/ total reads / peak length * 1mil)
# atac_homer_fpkm_filter.pl
# 13_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out -u out2 -v out3\n";
our ($opt_i, $opt_o, $opt_n, $opt_h);
getopts("i:o:n:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_n)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $count) = split ("\t", $line);
  my ($chr, $start, $end) = split ("_", $id);
  
  my $peak_length = $end - $start;

  my $fpkm = $count / $opt_n / $peak_length * 1000000 * 1000;
  
  print $out_fh "$chr\t$start\t$end\t$id\t$fpkm\n";
}

close $in_fh;
close $out_fh;
exit;

