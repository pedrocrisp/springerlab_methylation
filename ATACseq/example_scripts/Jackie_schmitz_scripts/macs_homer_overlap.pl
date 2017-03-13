#! /usr/bin/perl -w
# make file of all overlapping peaks (in MACs and HOMER peak calling)
# 28_Jan_2017
# Jaclyn_Noshay

# Input peaks overlapping 100bp bins files for macs and homer peaks and determine which bins have peaks in both

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
  my ($bin, $count) = split("\t", $line);

  $bins{$bin} = 'NA';
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($bin, $count) = split("\t", $line2);
  
  if (exists ($bins{$bin})) {
    print $out_fh "$bin\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
    
			  


