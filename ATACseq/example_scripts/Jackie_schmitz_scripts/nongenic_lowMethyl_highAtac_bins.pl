#! /usr/bin/perl -w
# take lowMethyl highAtac venn file and remove genic bins (found from bedtools intersect of 100bp bins and v4 genes gff file)
# 27_Feb_2017
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

# create hash of genic and read through with methyl_atac_venn bins... print file of non-genic bins
open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "bin\tCHG_ratio\tatac_fpkm\thomer_peak\tvenn_call\n";

my %genic;
while (my $line = <$in_fh>) {
  chomp $line;

  # 1	B73v4	bin	44201	44300	.	.	.	ID=bin1-443	1	44289	49837	+	ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene	11

  my ($chr, undef, undef, $start, $end, undef) = split ("\t", $line);

  my $binstart = $start - 1;
  my $bin = $chr . "_" . $binstart . "_" . $end;

  $genic{$bin} = 'NA';
}

my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;

  #1_1500_1600	NA	14.7848083812798	YES	uniq_atac

  my ($id, undef) = split ("\t", $line2);

  if (! (exists ($genic{$id}))) {
    print $out_fh "$line2\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
