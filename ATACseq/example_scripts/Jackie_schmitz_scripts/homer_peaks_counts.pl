#! /usr/bin/perl -w
# script to make counts file for homer peaks
# homer_peaks_counts.pl
# 13_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out \n";
our ($opt_i, $opt_o, $opt_I, $opt_h);
getopts("i:o:I:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (defined $opt_h) ) {
  print "$usage";
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

# create hash of homer peaks
my %peaks;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($chr, $start, $end, undef) = split ("\t", $line);
  
  my $peak_id = $chr . "_" . $start . "_" . $end;
  $peaks{$peak_id} = 0;
}
close $in_fh;


while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, $start, $end, undef) = split ("\t", $line2);
  
  my $peak_id2 = $chr . "_" . $start . "_" . $end;
  
  if (exists ($peaks{$peak_id2})) {
    $peaks{$peak_id2} += 1;
  }
}

while (my ($key, $value) = each %peaks) {
  print $out_fh "$key\t$value\n";
}

close $in2_fh;
close $out_fh;
exit;

