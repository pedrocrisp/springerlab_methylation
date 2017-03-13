#! /usr/bin/perl -w
# add TE descriptors to atac te file (including complete TE coordinates, nested, etc.)
# atac_te_descriptors.pl
# 20_Jan_2017
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

print $out_fh "id\tchr\tstart(0-based)\tstop\tdistance\tcategory\tHistone_profile\tclass\tte_start\tte_end\toverlap\tratio_overlap\tte_class\tte_family\tte_id\tatac_te_dist\tte\tstart\tend\tstrand\tigv_coords\tdisjoined_length\tnested\tin_intron\tltr_coord\torder\n";

my %te;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  my ($id, $start, $end, $strand, $igv_coords, $disjoined_length, $nested, $in_intron, $ltr_coord, $order) = split ("\t", $line);

  # te	start	end	strand	igv_coords	disjoined_length	nested	in_intron	ltr_coord	order

  $te{$id} = $line;
  #print "$id\t";
}

my $header2 = <$in2_fh>;
while (my $line2 = <$in2_fh>) {
  chomp $line2;
  my ($id, $chr, $start, $end, $dist, $category, $histone, $class, $te_start, $te_end, $overlap, $ratio_overlap, $te_class, $te_family, $te_id, $atac_dist) = split ("\t", $line2);

  # ID	chr	start(0-based)	stop	distance	category	Histone_profile	class	overlap	te_class	te_family	te_id	atac_te_dist

  #print "$te_id\t";
  if (exists ($te{$te_id})) {
    print $out_fh "$line2\t$te{$te_id}\n";
  }
}

close $in_fh;
close $in2_fh;
close $out_fh;
exit;

