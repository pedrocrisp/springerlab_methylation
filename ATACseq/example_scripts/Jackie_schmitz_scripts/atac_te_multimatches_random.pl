#! /usr/bin/perl -w
# take output of bedtools closest (te and randomized atac data) and determine most overlap for multi te matches
# atac_te_file.pl
# 18_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -I in2\n";
our ($opt_i, $opt_o, $opt_I, $opt_h);
getopts("i:o:I:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my %multi;
my ($chr, $start, $end, $width, $te_chr, $te_start, $te_end, $te_id, $dist, $class, $overlap, $name);
while (my $line = <$in_fh>) {
  chomp $line;
  ($chr, $start, $end, $width, undef, $te_chr, undef, undef, $te_start, $te_end, undef, undef, undef, $te_id, $dist, $class, $overlap) = split ("\t", $line);

  #1	261551	262252	702	*	1	rtracklayer	sequence_feature	243770	251243	.	*	.	ID=RLG00001B73v400962;sup=RLG	10309	no_overlap	NA

  $name = $chr . "_" . $start . "_" . $end;

  $multi{$name}{$overlap} = $line;
}

my ($chr2, $start2, $end2, $width2, $te_chr2, $te_start2, $te_end2, $te_id2, $dist2, $class2, $overlap2, $name2);
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  ($chr2, $start2, $end2, $width2, undef, $te_chr2, undef, undef, $te_start2, $te_end2, undef, undef, undef, $te_id2, $dist2, $class2, $overlap2) = split ("\t", $line2);

  $name2 = $chr2 . "_" . $start2 . "_" . $end2;
  #print "$overlap2\t";

  if ($dist2 == 0) {
    if (exists ($multi{$name2})) {
      
      foreach my $key (keys %{$multi{$name2}})
	{	  
	  if ($overlap2 ne 'NA' && $overlap2 > $key) {
	    print $out_fh "$line2\n";
	    delete $multi{$name2};
	  }
	}
    }
  }
}

for $name ( sort keys %multi ) {
    for $overlap ( sort keys %{ $multi{$name} } ) {
         print $out_fh "$multi{$name}{$overlap}\n";
    }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
