#! /usr/bin/perl -w
# take subset of methylation file (CHG < 0.20) and atac file (fpkm > 10) and find overlapping bins for venn
# 24_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -I in2 -o out\n";
our ($opt_i, $opt_I, $opt_o, $opt_h);
getopts("i:I:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_I))|| (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

# create hash of methylation bins and read through with atac bins... print file with declaration of overlap or uniq
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "bin\tCHG_ratio\tatac_fpkm\thomer_peak\tvenn_call\n";

my %bins;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;

  # chr	start	end	CG_ratio	CG_C	CG_CT	CG_sites	CG_cov	CHG_ratio	CHG_C	CHG_CT	CHG_sites	CHG_cov	CHH_ratio	CHH_C	CHH_CT	CHH_sites	CHH_cov

  my ($chr, $start, $end, $CG_ratio, undef, undef, undef, $CG_cov, $CHG_ratio, undef, undef, undef, $CHG_cov, $CHH_ratio, undef, undef, undef, $CHH_cov) = split ("\t", $line);

  my $bin = $chr . "_" . $start . "_" . $end;

  if ($CHG_ratio ne "NA") {
    if ($CHG_ratio < 0.20) {
      $bins{$bin} = $CHG_ratio;
    }
  }
}

my $keys = keys %bins;
my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  # ID=bin1-25	2.95696167625595	YES

  my ($id, $fpkm, $homer) = split ("\t", $line2);

  $id =~ s/ID=bin//;
  my ($chr, $value) = split("-", $id);
  my $end = $value * 100;
  my $start = $end - 100;

  my $bin = $chr . "_" . $start ."_" . $end;

  if ($fpkm ne "N/A") {
    if ($fpkm > 10) {
      if (exists ($bins{$bin})) {
	print $out_fh "$bin\t$bins{$bin}\t$fpkm\t$homer\toverlap\n";
	delete $bins{$bin};
      }
      else {
	print $out_fh "$bin\tNA\t$fpkm\t$homer\tuniq_atac\n";
      }
    }
  }
}

while (my ($key, $value) = each %bins) {
  print $out_fh "$key\t$value\tNA\tNA\tuniq_CHG\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
