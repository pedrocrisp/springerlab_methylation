#! /usr/bin/perl -w
# script to combine methylation data file and atacseq data file by bin
# methylation_atacseq.pl
# 9_Jan_2017
# Jaclyn_Noshay

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
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tbinstart\tbinend\tCG_ratio\tCHG_ratio\tCHH_ratio\tatacseq_count\n";

my %bins;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $leaf, $root) = split ("\t", $line);
  
  # binid	B73_leaf_fpkm	B73_root_fpkm	leaf_root_ratio
  # ID=bin1-1000006	3.37938477286395	0.948671286253117	3.56222942744605
  
  $id =~ s/ID=bin//;
  my ($chr, $value) = split ("-", $id);
  my $binend = $value * 100;
  my $binstart = $binend - 99;
  my $binid = $chr . "_" . $binstart . "_" . $binend;
  #print "$value\t";
  #print "$binid\t";
  
  $bins{$binid} = $leaf;
}

my $binid2;
my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, $binstart, $binend, $CG, undef, undef, undef, $CHG, undef, undef, undef, $CHH, undef, undef, undef) = split ("\t", $line2);
  
  # 1	1200	1300	0.70	28	40	6	0.56	14	25	4	0.06	11	193	28
  
  my $binstart2 = $binstart + 1;
  $binid2 = $chr . "_" . $binstart2 . "_" . $binend;

  #print "$binid2\t";
  
  if (exists ($bins{$binid2})) {
    print $out_fh "$chr\t$binstart\t$binend\t$CG\t$CHG\t$CHH\t$bins{$binid2}\n";
  }
}

close $in_fh;
close $out_fh;
exit;

