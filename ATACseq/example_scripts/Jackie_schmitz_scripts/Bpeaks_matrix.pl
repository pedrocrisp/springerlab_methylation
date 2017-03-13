#! /usr/bin/perl -w
# make file of bins for Bpeaks with matrix file
# 29_Jan_2017
# Jaclyn_Noshay

# Input homer Bpeaks 100bp bins file and fpkm matrix file and output fpkm matrix file with only Bpeaks 100bp bins
# Also calculate B minus P

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

my $header = <$in2_fh>;
my $bp;
print $out_fh "binid\tchr\tstart\tend\t207_Leaf\t207_Root\tB73_Leaf\tB73_Root\tPHJ89_Leaf\tW22_Leaf\tW22_Root\tBminusPleaf\n";
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  # binid	207_Leaf	207_Root	B73_Leaf	B73_Root	PHJ89_Leaf	W22_Leaf	W22_Root

  my ($binid, $x207leaf, $x207root, $Bleaf, $Broot, $PHJleaf, $Wleaf, $Wroot) = split("\t", $line2);

  $binid =~ s/ID=bin//;
  my ($chr, $binnumber) = split ("-", $binid);
  my $binend = $binnumber * 100;
  my $binstart = $binend - 99;
  my $bin = $chr . "_" . $binstart . "_" . $binend;
  
  if ($Bleaf eq 'N/A' | $x207leaf eq 'N/A') {
      $bp = "N/A";
  }
  else {
    $bp = $Bleaf - $x207leaf;
  }
  
  if (exists ($bins{$bin})) {
    print $out_fh "$chr\t$binstart\t$binend\t$line2\t$bp\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;		      
