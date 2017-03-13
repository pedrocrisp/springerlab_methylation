#! /usr/bin/perl -w
# combine atac peaks file with information about te overlap
# atac_te_file.pl
# 13_Jan_2017
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

print $out_fh "ID\tchr\tstart(0-based)\tstop\tdistance\tcategory\tHistone_profile\tclass\tte_start\tte_end\toverlap\tratio_overlap\tte_class\tte_family\tte_id\tatac_te_dist\n";

my %peaks;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $chr, $start, $end, $distance, $category, $histone) = split ("\t", $line);
  
  my $name = $chr . "_" . $start . "_" . $end;
  $peaks{$name} = $line;
}

my $ratio;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($chr, $start, $end, undef, undef, undef, $te_start, $te_end, undef, undef, undef, $te_id, $dist, $class, $overlap) = split ("\t", $line2);
  
  #1	132224	133416	1	rtracklayer	sequence_feature	122983	128278	.	*	.	ID=RLC00032B73v400024;sup=RLC	3947
  
  my $name2 = $chr . "_" . $start . "_" . $end;
  my ($id, $id_name) = split(";", $te_id);
  $id =~ s/ID=//;
  $id_name =~ s/sup=//;
  my $fam = substr($id, 0, 8);

  if (! ($overlap eq 'NA')) {
    my $atac_length = $end - $start;
    $ratio = $overlap / $atac_length;
  }
  else {
    $ratio = 'NA';
  }

  if (exists ($peaks{$name2})) {
    print $out_fh "$peaks{$name2}\t$class\t$te_start\t$te_end\t$overlap\t$ratio\t$id_name\t$fam\t$id\t$dist\n";
    delete $peaks{$name2};
  }
}

while (my ($key, $value) = each %peaks) {
  print $out_fh "$value\tNA\tNA\tNA\tNA\tNA\tNA\n";
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;
