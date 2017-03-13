#! /usr/bin/perl -w
# create binomial file indicating 100bp bins of macs and homer peaks to make venn diagram 
# 28_Jan_2017
# Jaclyn_Noshay

# Input peaks overlapping 100bp bins files for macs and homer peaks and determine which bins have peaks in both, one, or none

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
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my %bins;
while (my $line = <$in_fh>) {
  chomp $line;

  # 1	B73v4	bin	1	100	.	.	.	ID=bin1-1

  my ($chr, $genome, $feature, $start, $end, undef, undef, undef, $bin) = split ("\t", $line);
  $bin =~ s/ID=//;

  my $binid = $chr . "_" . $start . "_" . $end;

  $bins{$binid} = 'NA';
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($bin, $count) = split("\t", $line2);
  
  if (exists ($bins{$bin})) {
    
			  


