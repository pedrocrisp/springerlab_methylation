#! /usr/bin/perl -w
# take 1:1 syntenic gene files (nogaps, indel, nonindel, etc) and make two bed files (B coordinates and P coordinates)
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
while (my $line = <$in_fh>) {
  chomp $line;
  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line);

  # 3:134437050-134442050
  my ($Bchr, $Bcoord) = split (":", $Bpos);
  my ($Pchr, $Pcoord) = split (":", $Ppos);
  my ($Bstart, $Bend) = split ("-", $Bcoord);
  my ($Pstart, $Pend) = split ("-", $Pcoord);

  print $out_fh "$Bchr\t$Bstart\t$Bend\t$B\t$Pchr\t$Pstart\t$Pend\t$P\n";
  print $out2_fh "$Pchr\t$Pstart\t$Pend\t$P\t$Bchr\t$Bstart\t$Bend\t$B\n";  
}

close $in_fh;
close $out_fh;
close $out2_fh;
exit;

