#! /usr/bin/perl -w
# take 1:1 syntenic genes file and only pull out "indels" 
# 30_Jan_2017
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
print $out_fh "$header\tvariant_type\n";
while (my $line = <$in_fh>) {
  chomp $line;

  # B73	PH207	B73loc	PH207loc	QueryPos	HitPos	Strand	REF	Alt	HSPalnLen	B73_Non-N	PH207_Non-N	FracID	FracAlnQuery	FracAlnHit	HitalnLen

  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line);

  my $Blength = length($ref);
  my $Plength = length($alt);

# criteria is that indel has length > 1bp and no Ns

  # if ($Blength > 1) {
  #   if ($b_nonN == 5000) {
  #     if ($p_nonN == 5000) {
  # 	print $out_fh "$line\n";
  #     }
  #   }
  # }

  # if ($Plength > 1) {
  #   if ($b_nonN == 5000) {
  #     if ($p_nonN == 5000) {
  # 	print $out_fh "$line\n";
  #     }
  #   }
  # }

  # if (! ($Blength > 1)) {
  #   if (! ($Plength > 1)) {
  #     if ($b_nonN == 5000) {
  # 	if ($p_nonN == 5000) {
  # 	  print $out2_fh "$line\n";
  # 	}
  #     }
  #   }
  # }

# criteria is tha indel is when BorP is "-" and BorP is > 10bp

  if ($Blength = "-" && $Plength > 10) {
    print $out_fh "$line\tinsertion\n";
  }
  
  if ($Plength = "-" && $Blength > 10) {
    print $out_fh "$line\tdeletion\n";
  }
  
  if (!($Blength = "-" && $Plength > 10)) {
    if (!($Plength = "-" && $Blength > 10)) {
      print $out2_fh "$line\tnon_indel\n";
    }
  }
}

close $in_fh;
close $out_fh;
close $out2_fh;
exit;
