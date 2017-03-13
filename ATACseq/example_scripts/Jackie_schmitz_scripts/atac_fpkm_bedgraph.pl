#! /usr/bin/perl -w
# script to convert fpkm matrix file into bedgraph viewable in IGV
# atac_fpkm_bedgraph.pl
# 22_Jan_2017
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

#read in files
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $leaf_count, $root_count) = split ("\t", $line);
  
  # ID=bin1-10000	N/A	2.84601385875935
  
  $id =~ s/ID=bin//;
  my ($chr, $bin) = split ("-", $id);
  my $binstart = $bin - 1;
  my $binend = $binstart + 100;
  
  print $out_fh "$chr\t$binstart\t$binend\t$leaf_count";
  
}

close $in_fh;
close $out_fh;
exit;

