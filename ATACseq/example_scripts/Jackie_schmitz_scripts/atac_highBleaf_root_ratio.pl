#! /usr/bin/perl -w
# take fpkm matrix file for atacseq uniq data and pull out high B73 leaf and high B73 root bins, replace 0 with 0.001, calculate ratio, print two files (one for high leaf and one for high root)
# atac_highBleaf_root_ratio.pl
# 10_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -u out2\n";
our ($opt_i, $opt_o, $opt_u, $opt_h);
getopts("i:o:u:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_u)) || (defined $opt_h) ) {
  print "$usage";
}

#read in files 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;
open (my $out2_fh, '>', $opt_u) || die;

print $out_fh "binid\tB73_leaf_fpkm\tB73_root_fpkm\tleaf_root_ratio\n";
print $out2_fh "binid\tB73_root_fpkm\tB73_leaf_fpkm\troot_leaf_ratio\n";

my ($leaf_fpkm, $root_fpkm, $ratio1, $ratio2);
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($binid, $B73_Leaf, $B73_Root) = split ("\t", $line);

  if ($B73_Leaf eq "N/A") {
    $B73_Leaf = 0.001;
  }
  if ($B73_Root eq "N/A") {
    $B73_Root = 0.001;
  }


  if ($B73_Leaf >= 3) {
    $leaf_fpkm = $B73_Leaf;
    $root_fpkm = $B73_Root;
    $ratio1 = $leaf_fpkm / $root_fpkm;
    print $out_fh "$binid\t$leaf_fpkm\t$root_fpkm\t$ratio1\n";
  }

  if ($B73_Root >= 3) {
    $root_fpkm = $B73_Root;
    $leaf_fpkm = $B73_Leaf;
    $ratio2 = $root_fpkm / $leaf_fpkm;
    print $out2_fh "$binid\t$root_fpkm\t$leaf_fpkm\t$ratio2\n";
  } 
}

close $in_fh;
close $out_fh;
close $out2_fh;
exit;
