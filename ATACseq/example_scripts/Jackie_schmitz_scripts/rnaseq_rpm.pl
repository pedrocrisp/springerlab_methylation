#! /usr/bin/perl -w
# rnaseq_rpm.pl
# This script will calculated rpm values from matrix (gene counts/(total gene row / 1,000)
# Jaclyn Noshay
# Jan_3_2017

use strict;
use Getopt::Long();
my $usage = "\nUsage: $0 --rpm_out <output file for rpm values> --matrix_in <matrix containing genes and counts> --help <help on usage of script>\n\n";
my ($matrix_in, $rpm_out, $help);

Getopt::Long::GetOptions('rpm_out=s' => \$rpm_out,
			 'matrix_in=s' => \$matrix_in,
			 'h|help' => \$help
			);

if (defined($help)) {
  die $usage;
}

# Open files and calculate rpm
open (my $matrix_in_fh, '<', $matrix_in) or die "Can't open file $matrix_in\n\n";

my $header = <$matrix_in_fh>;
my ($sum_s1, $sum_s2, $sum_s3, $sum_s4, $sum_s5, $sum_s6, $sum_s7, $sum_s8, $sum_s9, $sum_s10, $sum_s11, $sum_s12);

while (my $line = <$matrix_in_fh>) {
  chomp $line;
  my ($genes, $s1, $s2, $s3, $s4, $s5, $s6, $s7, $s8, $s9, $s10, $s11, $s12) = split ("\t", $line);
  $sum_s1 += $s1;
  $sum_s2 += $s2;
  $sum_s3 += $s3;
  $sum_s4 += $s4;
  $sum_s5 += $s5;
  $sum_s6 += $s6;
  $sum_s7 += $s7;
  $sum_s8 += $s8;
  $sum_s9 += $s9;
  $sum_s10 += $s10;
  $sum_s11 += $s11;
  $sum_s12 += $s12;
}
close $matrix_in_fh;

my $s1_mil = $sum_s1 / 1000000;
my $s2_mil = $sum_s2 / 1000000;
my $s3_mil = $sum_s3 / 1000000;
my $s4_mil = $sum_s4 / 1000000;
my $s5_mil = $sum_s5 / 1000000;
my $s6_mil = $sum_s6 / 1000000;
my $s7_mil = $sum_s7 / 1000000;
my $s8_mil = $sum_s8 / 1000000;
my $s9_mil = $sum_s9 / 1000000;
my $s10_mil = $sum_s10 / 1000000;
my $s11_mil = $sum_s11 / 1000000;
my $s12_mil = $sum_s12 / 1000000;

open (my $matrix_in_fh2, '<', $matrix_in) or die "Can't open file $matrix_in\n\n";
open (my $rpm_out_fh, '>', $rpm_out) or die "Can't open file $rpm_out\n\n";
my $header2 = <$matrix_in_fh2>;
print $rpm_out_fh "$header2";

while (my $line2 = <$matrix_in_fh2>) {
  chomp $line2;
  my ($genes, $s1, $s2, $s3, $s4, $s5, $s6, $s7, $s8, $s9, $s10, $s11, $s12) = split ("\t", $line2);
  my $s1_rpm = $s1/$s1_mil;
  my $s2_rpm = $s2/$s2_mil;
  my $s3_rpm = $s3/$s3_mil;
  my $s4_rpm = $s4/$s4_mil;
  my $s5_rpm = $s5/$s5_mil;
  my $s6_rpm = $s6/$s6_mil;
  my $s7_rpm = $s1/$s7_mil;
  my $s8_rpm = $s2/$s8_mil;
  my $s9_rpm = $s3/$s9_mil;
  my $s10_rpm = $s4/$s10_mil;
  my $s11_rpm = $s5/$s11_mil;
  my $s12_rpm = $s6/$s12_mil;
  
  print $rpm_out_fh "$genes\t$s1_rpm\t$s2_rpm\t$s3_rpm\t$s4_rpm\t$s5_rpm\t$s6_rpm\t$s7_rpm\t$s8_rpm\t$s9_rpm\t$s10_rpm\t$s11_rpm\t$s12_rpm\n";
}

close $matrix_in_fh2;
close $rpm_out_fh;

exit;
















