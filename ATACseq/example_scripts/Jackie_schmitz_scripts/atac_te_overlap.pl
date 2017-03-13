#! /usr/bin/perl -w
# take bedtools intersect between te and atac peaks file and separate into complete and partial overlap
# atac_te_overlap.pl
# 13_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -u out2\n";
our ($opt_i, $opt_o, $opt_u, $opt_t, $opt_h);
getopts("i:o:u:t:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_u)) || (!(defined $opt_t)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;
open (my $out2_fh, '>', $opt_u) || die;
open (my $out3_fh, '>', $opt_t) || die;

my $class;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $start, $end, undef, $family, $te, $te_start, $te_end, undef, undef, undef, $te_id, $overlap) = split ("\t", $line);

  # 1	262093	262904	1	HelitronScanner	helitron	251639	262224	.	-	.	ID=DHH00002B73v406099;Name=DHH00002B73v406099_Hip1_LCV5p16_LCV3p28	131

  my $atac_length = $end - $start;
  
  if ($overlap == $atac_length) {
    $class = "complete_overlap";
    print $out_fh "$line\n";
  }
  elsif ($overlap > $atac_length) {
    $class = "te_within";
  }
  else {
    $class = "partial_overlap";
    print $out2_fh "$line\n";
  }

  print $out3_fh "$line\t$class\n";
}

close $in_fh;
close $out_fh;
close $out2_fh;
exit;

