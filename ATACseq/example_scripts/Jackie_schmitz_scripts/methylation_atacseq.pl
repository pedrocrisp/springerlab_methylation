#! /usr/bin/perl -w
# script to combine methylation data file and atacseq data file by bin
# methylation_atacseq.pl
# 9_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in methylation and atacseq files
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tbinstart\tbinend\tCG_ratio\tCHG_ratio\tCHH_ratio\tatacseq_count\n";

my %bins;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $leaf_count, $root_count) = split ("\t", $line);
  
  # ID=bin1-1
  
  $id =~ s/ID=bin//;
  my ($chr, $bin) = split ("-", $id);
  my $binstart = $bin - 1;
  my $binend = $binstart + 100;
  my $binid = $chr . "_" . $binstart . "_" . $binend;
  #print "$binid\t";
  
  $bins{$binid} = $leaf_count;
}

my $binid2;
my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, $binstart, $binend, $CG, undef, undef, undef, undef, $CHG, undef, undef, undef, undef, $CHH, undef, undef, undef, undef) = split ("\t", $line2);
  
  # chr	start	end	CG_ratio	CG_C	CG_CT	CG_sites	CG_cov	CHG_ratio	CHG_C	CHG_CT	CHG_sites	CHG_cov	CHH_ratio	CHH_C	CHH_CT	CHH_sites	CHH_cov

  $binid2 = $chr . "_" . $binstart . "_" . $binend;
  
  if (exists ($bins{$binid2})) {
    print $out_fh "$chr\t$binstart\t$binend\t$CG\t$CHG\t$CHH\t$bins{$binid2}\n";
  }
}

close $in_fh;
close $out_fh;
exit;

