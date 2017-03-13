#! /usr/bin/perl -w
# use bedtools intersect from bins peaks and gene gff output... pull out data to use for metaplot
# TEfamily_macs.pl
# 2_Jan_2017
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

#read in bedtools intersect file 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tbinstart\tbinstop\tCG\tCHG\tCHH\tbin\thomer_chr\thomer_start\thomer_end\thomer_size\tpeak\tdistance\trelative_distance\treal_distance\n";

my $strand_dist;
my $relative_dist;
my $real_dist;
my $homersize;
my $diff;
my $decimal;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $binstart, $binend, $CG, $CHG, $CHH, $homer_chr, $homer_start, $homer_end, $dist) = split ("\t", $line);

  # 1	6354	6812	1	6300	6400	0.88	0.78	0.08	0

  my $bin = $chr . "_" . $binstart  . "_" . $binend;
  my $peak = $homer_chr . "_" . $homer_start . "_" . $homer_end;
  $homersize = ($homer_end - $homer_start);

  if ($dist > 0) {
    $relative_dist = $dist + 1000;
    $real_dist = $dist;
  }
  elsif ($dist == 0) {
    if ($binstart < $homer_start) {
      $diff = $binend - $homer_start;
    }
    if ($binstart > $homer_start) {
      $diff = $homer_end - $binstart;
    }
    $decimal = $diff / $homersize;
    $relative_dist = $decimal * 1000;
    $real_dist = $diff;
  }
  else {
    $real_dist = $dist;
    $relative_dist = $dist;
  }

  print $out_fh "$chr\t$binstart\t$binend\t$CG\t$CHG\t$CHH\t$bin\t$homer_chr\t$homer_start\t$homer_end\t$homersize\t$peak\t$dist\t$relative_dist\t$real_dist\n";
}

close $in_fh;
close $out_fh;
exit;

