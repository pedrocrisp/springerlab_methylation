#! /usr/bin/perl -w
# take bedtools intersect (between syntenic and atac) output files and compare genes across two files
# 1_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out\n";
our ($opt_i, $opt_o, $opt_I, $opt_h);
getopts("i:o:I:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "Bchr\tBstart\tBend\tBgene\tPchr\tPstart\tPend\tPgene\tBpeak_chr\tBpeak_start\tBpeak_end\ttranscript\tv4B_207leaf_log2fc\tv4B_207leaf_pval\n";


my %genes;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($Bchr, $Bstart, $Bend, $Bgene, $Pchr, $Pstart, $Pend, $Pgene, $chr2, $peakstart, $peakend) = split("\t", $line);

  # Bchr	Bstart	Bend	Bgene	Pchr	Pstart	Pend	Pgene	Bpeak_chr	Bpeak_start	Bpeak_end

  $genes{$Bgene} = $line;
}

while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($transcript, $log2fc, $pval) = split("\t", $line2);

  # transcript	v4B_207leaf_log2fc	v4B_207leaf_pval

  $transcript =~ s/gene://;

  if (exists ($genes{$transcript})) {
    print $out_fh "$genes{$transcript}\t$line2\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;

