#! /usr/bin/perl -w
# take 1:1 syntenic genes file and only pull out genes/variants where alignment includes no Ns (no gaps)
# 1_Feb_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -u out2\n";
our ($opt_i, $opt_o, $opt_u, $opt_h);
getopts("i:o:u:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_u)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;
open (my $out2_fh, '>', $opt_u) || die;

my $header = <$in_fh>;
print $out_fh "$header\n";
while (my $line = <$in_fh>) {
  chomp $line;

  # B73	PH207	B73loc	PH207loc	QueryPos	HitPos	Strand	REF	Alt	HSPalnLen	B73_Non-N	PH207_Non-N	FracID	FracAlnQuery	FracAlnHit	HitalnLen

  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line);

 if ($b_nonN == 5000) {
   if ($p_nonN == 5000) {
     print $out_fh "$line\n";
   }
 }

  else {
    print $out2_fh "$line\n";
  }
}

close $in_fh;
close $out_fh;
close $out2_fh;
exit;
