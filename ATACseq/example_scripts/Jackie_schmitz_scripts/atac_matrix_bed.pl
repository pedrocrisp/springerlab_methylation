#! /usr/bin/perl -w
# take atacseq_counts_matrix and split into bed files
# atac_rep_combine.pl
# 19_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -u out2 -b out3 -r out4 -p out5 -w out6 -x out7 -\n";
our ($opt_i, $opt_o, $opt_u, $opt_b, $opt_r, $opt_p, $opt_w, $opt_x, $opt_h);
getopts("i:o:u:b:r:p:w:x:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_u)) || (!(defined $opt_b)) || (!(defined $opt_r)) || (!(defined $opt_p)) || (!(defined $opt_w)) || (!(defined $opt_x)) || (defined $opt_h) ) {
  print "$usage";
}

#read in files 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;
open (my $out2_fh, '>', $opt_u) || die;
open (my $out3_fh, '>', $opt_b) || die;
open (my $out4_fh, '>', $opt_r) || die;
open (my $out5_fh, '>', $opt_p) || die;
open (my $out6_fh, '>', $opt_w) || die;
open (my $out7_fh, '>', $opt_x) || die;

print $out_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out2_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out3_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out4_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out5_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out6_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";
print $out7_fh "chr\tbinstart\tbinstop\tdot\tdot\tdot\tbinid\tfpkm\n";

my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($binid, $x207leaf, $x207root, $bleaf, $broot, $phjleaf, $wleaf, $wroot) = split ("\t", $line);

  my $bin = $binid;
  $bin =~ s/ID=bin//;
  my ($chr, $end) = split ("-", $bin);
  my $binstop = $end * 100;
  my $binstart = $binstop - 99;
  #print "$end\t$binstop\t$binstart\n";
  
  print $out_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$x207leaf\n";
  print $out2_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$x207root\n";
  print $out3_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$bleaf\n";
  print $out4_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$broot\n";
  print $out5_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$phjleaf\n";
  print $out6_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$wleaf\n";
  print $out7_fh "$chr\t$binstart\t$binstop\t.\t.\t.\t$binid\t$wroot\n";
}

close $in_fh;
close $out_fh;
close $out2_fh;
close $out3_fh;
close $out4_fh;
close $out5_fh;
close $out6_fh;
close $out7_fh;
exit;

  
    
