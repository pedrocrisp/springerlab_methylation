#! /usr/bin/perl -w
# take indel or nonindel syntenic gene files and 5kb upstream
# 30_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -u out2\n";
our ($opt_i, $opt_o, $opt_u,  $opt_h);
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
my $blank = <$in_fh>;

print $out_fh "Bchr\tB_5kb\tB_start\tB_end\tBgene\n";
print $out2_fh "Pchr\tP_5kb\tP_start\tP_end\tPgene\n";
while (my $line = <$in_fh>) {
  chomp $line;

  # B73	PH207	B73loc	PH207loc	QueryPos	HitPos	Strand	REF	Alt	HSPalnLen	B73_Non-N	PH207_Non-N	FracID	FracAlnQuery	FracAlnHit	HitalnLen

  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line);

  # 2:353456-358456
  my ($Bchr, $Bcoord) = split(":", $Bpos);
  my ($Bstart, $Bend) = split("-", $Bcoord);

  my ($Pchr, $Pcoord) = split(":", $Ppos);
  my ($Pstart, $Pend) = split("-", $Pcoord);

  my $Bup = $Bstart - 5000;
  my $Pup = $Pstart - 5000;

  if ($Bup >= 0) {
    print $out_fh "$Bchr\t$Bup\t$Bstart\t$Bend\t$B\n";
  }

  if ($Pup >= 0) {
    print $out2_fh "$Pchr\t$Pup\t$Pstart\t$Pend\t$P\n";
  }
}
  
close $in_fh;
close $out_fh;
close $out2_fh;
exit;
