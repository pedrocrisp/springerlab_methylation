#! /usr/bin/perl -w
# deseq_compile.pl
# This script will take multiple deseq output files and compile into single matrix file with log2fc
# Jaclyn Noshay
# Jan_5_2017

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in bedtools intersect file 
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "transcript\tv4leaf_log2fc\tv4leaf_pval\tv4root_log2fc\tv4root_pval\tv4B_207leaf_log2fc\tv4B_207leaf_pval\tv4B_207root_log2fc\tv4B_207root_pval\tv4B_PHJleaf_log2fc\tv4B_PHJleaf_pval\tv4B_W22leaf_log2fc\tv4B_W22leaf_pval\tv4B_W22root_log2fc\tv4B_W22root_pval\tv4Bleaf_Broot_log2fc\tv4Bleaf_Broot_pval\n";

my %log;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($gene, $v4leaf, $v4root, $x207leaf, $x207root, $PHJleaf, $W22leaf, $W22root, $Bleafroot) = split ("\t", $line);

# transcript	v4leaf	v4root	v4B_207leaf	v4B_207root	v4B_PHJleaf	v4B_W22leaf	v4B_W22root v4Bleaf_Broot

  my $log2fc = $v4leaf . "_" . $v4root . "_" . $x207leaf . "_" . $x207root . "_" . $PHJleaf . "_" . $W22leaf . "_" . $W22root . "_" . $Bleafroot;

  $log{$gene} = $log2fc;
}


my ($v4leaf_log, $v4root_log, $x207leaf_log, $x207root_log, $PHJleaf_log, $W22leaf_log, $W22root_log, $Bleafroot_log);
my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($gene2, $v4leaf_p, $v4root_p, $x207leaf_p, $x207root_p, $PHJleaf_p, $W22leaf_p, $W22root_p, $Bleafroot_p) = split ("\t", $line2);

  if (exists ($log{$gene2})) {
    ($v4leaf_log, $v4root_log, $x207leaf_log, $x207root_log, $PHJleaf_log, $W22leaf_log, $W22root_log, $Bleafroot_log) = split ("_", $log{$gene2});
  }

  print $out_fh "$gene2\t$v4leaf_log\t$v4leaf_p\t$v4root_log\t$v4root_p\t$x207leaf_log\t$x207leaf_p\t$x207root_log\t$x207root_p\t$PHJleaf_log\t$PHJleaf_p\t$W22leaf_log\t$W22leaf_p\t$W22root_log\t$W22root_p\t$Bleafroot_log\t$Bleafroot_p\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;

