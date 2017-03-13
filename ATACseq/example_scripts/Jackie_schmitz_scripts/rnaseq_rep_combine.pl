#! /usr/bin/perl -w
# combine reps for rnaseq rpm data set
# rnaseq_rep_combine.pl
# 23_Jan_2017
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

print $out_fh "transcript\t207_Leaf\t207_Root\tB73_Leaf\tB73_Root\tPHJ89_Leaf\tW22_Leaf\tW22_Root\n";

my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($transcript, $x207_Leaf_rep1, $x207_Leaf_rep2, $x207_Root_rep1, $x207_Root_rep2, $B73_Leaf_rep1, $B73_Leaf_rep2, $B73_Root_rep1, $B73_Root_rep2, $PHJ89_Leaf_rep1, $W22_Leaf_rep1, $W22_Leaf_rep2, $W22_Root_rep1) = split ("\t", $line);

# transcript	207_Leaf_rep1	207_Leaf_rep2	207_Root_rep1	207_Root_rep2	B73_Leaf_rep1	B73_Leaf_rep2	B73_Root_rep1	B73_Root_rep2	PHJ89_Leaf_rep1	W22_Leaf_rep1	W22_Leaf_rep2	W22_Root_rep1

  my $x207_Leaf = ($x207_Leaf_rep1 + $x207_Leaf_rep2) / 2;
  my $x207_Root = ($x207_Root_rep1 + $x207_Root_rep2) / 2;
  my $B73_Leaf = ($B73_Leaf_rep1 + $B73_Leaf_rep2) / 2;
  my $B73_Root = ($B73_Root_rep1 + $B73_Root_rep2) / 2;
  my $PHJ89_Leaf = $PHJ89_Leaf_rep1;
  my $W22_Leaf =  ($W22_Leaf_rep1 + $W22_Leaf_rep2) / 2;
  my $W22_Root = $W22_Root_rep1;

    print $out_fh "$transcript\t$x207_Leaf\t$x207_Root\t$B73_Leaf\t$B73_Root\t$PHJ89_Leaf\t$W22_Leaf\t$W22_Root\n";
}

close $in_fh;
close $out_fh;
exit;

  
    
