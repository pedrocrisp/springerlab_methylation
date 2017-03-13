#! /usr/bin/perl -w
# this script will take a 100bp count file and convert it to a .igv file in order to create a .tdf file for IGV viewing
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

#read in counts file
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tstart\tend\tfeature\tcount\n";

while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $count) = split ("\t", $line);

# 1_33982501_33982600	0

  my ($chr, $start, $end) = split ("_", $id);

  my $feature = "atacseq_count";

  print $out_fh "$chr\t$start\t$end\t$feature\t$count\n";
}

close $in_fh;
close $out_fh;
exit;

