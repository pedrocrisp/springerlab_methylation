#! /usr/bin/perl -w
# homer peak call filtering (minimum peak density of 0.1 (# of reads divided by peak length))
# atac_rep_combine.pl
# 9_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -p param\n";
our ($opt_i, $opt_o, $opt_p, $opt_h);
getopts("i:o:p:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_p)) || (defined $opt_h) ) {
  print "$usage";
}

#read in files 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $start, $end) = split ("\t", $line);

  my $length = $end - $start;
  my $density = $opt_p / $length;
  
  if ($density > 0.1) {
    print $out_fh "$line\n";
  }
}

close $in_fh;
close $out_fh;
exit;

  
