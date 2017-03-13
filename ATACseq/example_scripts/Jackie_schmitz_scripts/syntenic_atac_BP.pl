#! /usr/bin/perl -w
# take bedtools intersect (between syntenic and atac) output files and compare genes across two files
# 1_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out -b out2 -p out3\n";
our ($opt_i, $opt_o, $opt_I, $opt_b, $opt_p, $opt_h);
getopts("i:o:I:b:p:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (!(defined $opt_b)) || (!(defined $opt_p)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;
open (my $out2_fh, '>', $opt_b) || die;
open (my $out3_fh, '>', $opt_p) || die;

print $out_fh "Bchr\tBstart\tBend\tBgene\tBpeak_chr\tBpeak_start\tBpeak_end\tPchr\tPstart\tPend\tPgene\tPpeak_chr\tPpeak_start\tPpeak_end\n";
print $out2_fh "Bchr\tBstart\tBend\tBgene\tPchr\tPstart\tPend\tPgene\tBpeak_chr\tBpeak_start\tBpeak_end\n";
print $out3_fh "Pchr\tPstart\tPend\tPgene\tBchr\tBstart\tBend\tBgene\tPpeak_chr\tPpeak_start\tPpeak_end\n";

my %genes;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($Bchr, $Bstart, $Bend, $Bgene, $Pchr, $Pstart, $Pend, $Pgene, $chr2, $peakstart, $peakend) = split("\t", $line);

  # 2	1246628	1251628	Zm00001d001804	2	1247319	1247520

  $genes{$Bgene} = $line;
}

my $keys = keys %genes;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($Pchr, $Pstart, $Pend, $Pgene, $Bchr, $Bstart, $Bend, $Bgene, $chr2, $peakstart, $peakend) = split("\t", $line2);

  if (exists ($genes{$Bgene})) {
    print $out_fh "$genes{$Bgene}\t$line2\n";
    delete $genes{$Bgene};
  }
  else {
    print $out3_fh "$line2\n";
    delete $genes{$Bgene};
  }
}

while (my ($key, $value) = each %genes) {
  print $out2_fh "$value\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
close $out2_fh;
close $out3_fh;
exit;

