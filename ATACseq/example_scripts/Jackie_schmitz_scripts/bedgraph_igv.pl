#! /usr/bin/perl -w
# this script will take a bedgraph file and convert it to a .igv file in order to create a .tdf file for IGV viewing
# sam_igv.pl
# 9_Jan_2017
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

#read in bedgraph file
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tstart\tend\tfeature\tcount\n";

while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $start, $end, $value) = split ("\t", $line);

# 1	0	3	1

  my $feature = "rnaseq_bedgraph";

  print $out_fh "$chr\t$start\t$end\t$feature\t$value\n";
}

close $in_fh;
close $out_fh;
exit;

