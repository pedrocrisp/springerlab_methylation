#! /usr/bin/perl -w
# take rnaseq rpm data and deseq data and compile into single matrix by gene
# rnaseq_rpm_deseq_compile.pl
# 23_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out -\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (defined $opt_h) ) {
  print "$usage";
}

#read in files (rpm data and deseq data)
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my $header = <$in_fh>;
my $header2 = <$in2_fh>;
print $out_fh "transcript\tx207leaf\tx207root\tbleaf\tbroot\tphjleaf\twleaf\twroot\ttranscript\tv4leaf_log2fc\tv4leaf_pval\tv4root_log2fc\tv4root_pval\tv4B_207leaf_log2fc\tv4B_207leaf_pval\tv4B_207root_log2fc\tv4B_207root_pval\tv4B_PHJleaf_log2fc\tv4B_PHJleaf_pval\tv4B_W22leaf_log2fc\tv4B_W22leaf_pval\tv4B_W22root_log2fc\tv4B_W22root_pval\tv4Bleaf_Broot_log2fc\tv4Bleaf_Broot_pval\n";

my %genes;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($transcript, $x207leaf, $x207root, $bleaf, $broot, $phjleaf, $wleaf, $wroot) = split ("\t", $line);

# transcript	207_Leaf	207_Root	B73_Leaf	B73_Root	PHJ89_Leaf	W22_Leaf	W22_Root

  $genes{$transcript} = $line;
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($transcript, undef) = split ("\t", $line2);

    # transcript	v4leaf_log2fc	v4leaf_pval	v4root_log2fc	v4root_pval	v4B_207leaf_log2fc	v4B_207leaf_pval	v4B_207root_log2fc	v4B_207root_pval	v4B_PHJleaf_log2fc	v4B_PHJleaf_pval	v4B_W22leaf_log2fc	v4B_W22leaf_pval	v4B_W22root_log2fc	v4B_W22root_pval	v4Bleaf_Broot_log2fc	v4Bleaf_Broot_pval

  if (exists ($genes{$transcript})) {
    print $out_fh "$genes{$transcript}\t$line2\n";
  }
}
			    
close $in_fh;
close $in2_fh;
close $out_fh;
exit;

  
    
