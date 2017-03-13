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
  $bin =~ s/ID=//;

  my $binid = $chr . "_" . $start . "_" . $end;

  $bins{$binid} = 0;
}

my $keys = keys %bins;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  
  # 10	23620951	23621051	HWI-D00238R:193:H9GR0ADXX:1:1101:1340:2217	1	-	23620951	23621051	0,0,0	1	100,	0,	10	B73v4	bin	23620901	23621000	.	.	.	ID=bin10-236210

  my @bedtools = split ("\t", $line2);
  my $chr = $bedtools[12];
  my $start = $bedtools[15];
  my $end = $bedtools[16];
  my $id = $bedtools[20];
  my (undef, $bin) = split ("=", $id);

  my $binid = $chr . "_" . $start . "_" . $end;

  if (exists ($bins{$binid})) {
    $bins{$binid} += 1;
  }

  #print $out_fh "$binid\t$bins{$binid}\n";
}

while (my ($key, $value) = each %bins) {
  print $out_fh "$key\t$value\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
