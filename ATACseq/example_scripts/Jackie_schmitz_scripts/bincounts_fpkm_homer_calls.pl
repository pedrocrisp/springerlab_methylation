#! /usr/bin/perl -w
# create file with counts for every 100bp bin and add column to declare if peak is in HOMER peak (input 100bp bin counts file and homer bins file)
# 23_Feb_2017
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

#read in gff 100bp bin file and bedtools intersect file
#create hash of gff file and then read through and count from bedtools file
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "id\tcount\thomer\n";
my %bins;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;

  # ID=bin1-1000	N/A	0.632447524168745

  my ($bin, $leaf_count, $root_count) = split ("\t", $line);

  $bins{$bin} = $leaf_count;
}

my $keys = keys %bins;
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  # ID=bin1-25	4

  my ($id, $count) = split ("\t", $line2);

  if (exists ($bins{$id})) {
    print $out_fh "$id\t$bins{$id}\tYES\n";
    delete $bins{$id};
  }
}

while (my ($key, $value) = each %bins) {
  print $out_fh "$key\t$value\tNO\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
