#! /usr/bin/perl -w
# create file with counts for every 100bp bin 
# 28_Nov_2016
# Jaclyn_Noshay

# Input bedtools intersect file and input gff 100bp bins file and count how many reads for each 100bp bin

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

my %bins;
my $count = 0;
while (my $line = <$in_fh>) {
  chomp $line;

  # 1	B73v4	bin	1	100	.	.	.	ID=bin1-1

  my ($chr, $genome, $feature, $start, $end, undef, undef, undef, $bin) = split ("\t", $line);
  #$bin =~ s/ID=//;

  my $binid = $chr . "_" . $start . "_" . $end;

  $bins{$binid} = 0;
}

my $keys = keys %bins;
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  # 1	1257	2608	256	1169	1351	0.8652850	1	2401	2500	.	.	.	ID=bin1-25	4

  my ($peakchr, $peakstart, $peakend, undef, undef, undef, undef, $chr, $start, $end, undef, undef, undef, $id, $count) = split ("\t", $line2);

  my $binid = $chr . "_" . $start . "_" . $end;

  if ($start >= $peakstart && $end <= $peakend) {
    if (exists ($bins{$binid})) {
      print $out_fh "$id\t$count\n";
    }
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
