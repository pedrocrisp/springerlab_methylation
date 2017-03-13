#! /usr/bin/perl -w
# take bedtools intersect between te and atac peaks file and separate into complete and partial overlap
# atac_te_overlap.pl
# 13_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out \n";
our ($opt_i, $opt_o, $opt_h);
getopts("i:o:u:t:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

my ($class, $overlap);
while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $start, $end, $width, undef, $te_chr, undef, undef, $te_start, $te_end, undef, undef, undef, $te_id, $dist) = split ("\t", $line);

  # 1	261551	262252	702	*	1	HelitronScanner	helitron	251639	262224	.	-	.	ID=DHH00002B73v406099;Name=DHH00002B73v406099_Hip1_LCV5p16_LCV3p28	0

  my $atac_length = $end - $start;
  my $te_length = $te_end - $te_start;

  if ($dist == 0) {
    if ($te_start > $start && $te_start < $end) {
      $overlap = $end - $te_start;
    }
    if ($te_end > $start && $te_end < $end) {
      $overlap = $te_end - $start;
    }
    if ($start < $te_start && $end > $te_end) {
      $overlap = $te_length;
    }
    if ($start > $te_start && $end < $end) {
      $overlap = $atac_length;
    }
  }
  else {
    $overlap = 'NA';
  }
  
  if ($dist == 0) {  
    if ($te_length == $atac_length) {
      $class = "complete_overlap";
    }
    elsif ($te_length < $atac_length) {
      $class = "te_within";
    }
    else {
      $class = "partial_overlap";
    }
  }
  else {
    $class = "no_overlap";
  }
  
  print $out_fh "$line\t$class\t$overlap\n";
}

close $in_fh;
close $out_fh;
exit;

