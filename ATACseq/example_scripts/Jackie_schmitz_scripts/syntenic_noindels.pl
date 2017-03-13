#! /usr/bin/perl -w
# take 1:1 syntenic genes file and only pull out "indels" (variants greater than 1 base)
# 30_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out -I in2\n";
our ($opt_i, $opt_o, $opt_I, $opt_u, $opt_h);
getopts("i:o:I:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (!(defined $opt_I)) || (defined $opt_h) ) {
  print "$usage";
  #exit;
}

open (my $in_fh, '<', $opt_i) || die;
open (my $in2_fh, '<', $opt_I) || die;
open (my $out_fh, '>', $opt_o) || die;

my $header = <$in_fh>;
my $blank = <$in_fh>;
print $out_fh "$header\n";

my %genes;
while (my $line = <$in_fh>) {
  chomp $line;

  # B73	PH207	B73loc	PH207loc	QueryPos	HitPos	Strand	REF	Alt	HSPalnLen	B73_Non-N	PH207_Non-N	FracID	FracAlnQuery	FracAlnHit	HitalnLen

  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line);

  $genes{$B} = $line;
}

my $header2 = <$in2_fh>;
my $blank2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($B, $P, $Bpos, $Ppos, $query_pos, $hit_pos, $strand, $ref, $alt, $HSPalnLet, $b_nonN, $p_nonN, $fracID, $fracAlnQuery, $fracAlnHit, $HitaltLen) = split("\t", $line2);

  if (! (exists ($genes{$B}))) {
    print $out_fh "$line2\n";
  }
}

close $in_fh;
close $out_fh;
close $in2_fh;
exit;
