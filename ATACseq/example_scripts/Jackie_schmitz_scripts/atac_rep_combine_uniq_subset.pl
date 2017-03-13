#! /usr/bin/perl -w
# combine 100bp atacseq counts for rep 1 and rep 2 for each genotype and tissue
# atac_rep_combine.pl
# 18_Jan_2017
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

#read in files 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

my ($x207_Leaf, $x207_Root, $B73_Leaf, $B73_Root, $PHJ89_Leaf, $W22_Leaf, $W22_Root);
my ($x207_Leaf2, $x207_Root2, $B73_Leaf2, $B73_Root2, $PHJ89_Leaf2, $W22_Leaf2, $W22_Root2);

print $out_fh "binid\tB73_Leaf\tB73_Root\n";

my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($binid, $B73_Leaf_rep1, $B73_Leaf_rep2, $B73_Root_rep1, $B73_Root_rep2) = split ("\t", $line);
  
  if (!($B73_Leaf_rep1 eq "N/A")) {
    if (!($B73_Leaf_rep2 eq "N/A")) {
      $B73_Leaf = ($B73_Leaf_rep1 + $B73_Leaf_rep2) / 2;
      $B73_Leaf2 = ($B73_Leaf / ((12197384 + 11475564)/2)) * 10000000;

    }
  }
  else {
    $B73_Leaf = 'N/A';
    $B73_Leaf2 = 'N/A';
  }
  
  if (!($B73_Root_rep1 eq "N/A")) {
    if (!($B73_Root_rep2 eq "N/A")) {
      $B73_Root = ($B73_Root_rep1 + $B73_Root_rep2) / 2;
      $B73_Root2 = ($B73_Root / ((22785968 + 8837209)/2)) * 10000000;
    }
  }
  else {
    $B73_Root = 'N/A';
    $B73_Root2 = 'N/A';
  }
  
  #print $out_fh "$binid\t$B73_Leaf\t$B73_Root\n";
  print $out_fh "$binid\t$B73_Leaf2\t$B73_Root2\n";
  
}

close $in_fh;
close $out_fh;
exit;

  
    
