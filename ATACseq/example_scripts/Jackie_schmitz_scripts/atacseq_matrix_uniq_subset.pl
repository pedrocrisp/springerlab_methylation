#! /usr/bin/perl -w
# atacseq_matrix.pl
# This script will make a matrix for the schmitz ATACseq dataset
# Jaclyn Noshay
# Jan_18_2017

use strict;
use Getopt::Long();
my $usage = "\nUsage: $0 --matrix_out <output file for matrix> --bins <list of 100bp bins> --help <help on usage of script>\n\n";
my ($matrix_out, $bins, $help);

Getopt::Long::GetOptions('matrix_out=s' => \$matrix_out,
			 'bins=s' => \$bins,
			 'h|help' => \$help
			);

if (defined($help)) {
    die $usage;
}

# Open genes and make hash of genes to be filled by files
open (my $bins_fh, '<', $bins) or die "Can't open file $bins\n\n";
my %hoh;

while (my $line = <$bins_fh>) {
    chomp $line;
    my ($chr, $genotype, $bin, $start, $end, undef, undef, undef, $id) = split ("\t", $line);
# 1	B73v4	bin	1	100	.	.	.	ID=bin1-1

    $hoh{$id}{'B73_Leaf_rep1'} = 'N/A';
    $hoh{$id}{'B73_Leaf_rep2'} = 'N/A';  
    $hoh{$id}{'B73_Root_rep1'} = 'N/A';  
    $hoh{$id}{'B73_Root_rep2'} = 'N/A';      
}

# Hard code in the files want to include in the matrix
# Could also do this by putting all transcript files into one directory and use open dir
my @files = qw | /home/springer/nosha003/schmitz/atacseq/counts/B73_Leaf_rep1_uniq.bowtie.100bpcounts.bed
		 /home/springer/nosha003/schmitz/atacseq/counts/B73_Leaf_rep2_uniq.bowtie.100bpcounts.bed
		 /home/springer/nosha003/schmitz/atacseq/counts/B73_Root_rep1_uniq.bowtie.100bpcounts.bed
		 /home/springer/nosha003/schmitz/atacseq/counts/B73_Root_rep2_uniq.bowtie.100bpcounts.bed
	       |;

my $file_fh;

# Open output file
open (my $matrix_out_fh, '>', $matrix_out) or die "\nCan't open file $matrix_out\n\n";

# Print header line (hard coded for this example)
print $matrix_out_fh "bins";

foreach my $file (@files) {
    if ($file =~ /\/home\/springer\/nosha003\/schmitz\/atacseq\/counts\/(.+)\_uniq.bowtie.100bpcounts.bed/) {
	my $name = $1;
	#print "$name\n";
	print $matrix_out_fh "\t$name";

	open ($file_fh, '<', $file) or die "Can't open file $file\n\n";
	my $header = <$file_fh>;
	while (my $line = <$file_fh>) {
	  # chr	binstart	binstop	dot	dot	dot	binid	count
	    chomp $line;
	    my ($chr, $binstart, $binstop, undef, undef, undef, $binid, $count) = split ("\t", $line);
	    if (exists $hoh{$binid}) {
		$hoh{$binid}{$name} = $count;
	      }
	  }
      }
  }
print $matrix_out_fh "\n";

# Fill in matrix
foreach my $key (sort {$a cmp $b} keys %hoh) {
    my $e = $hoh{$key}{'B73_Leaf_rep1'};
    my $f = $hoh{$key}{'B73_Leaf_rep2'};
    my $g = $hoh{$key}{'B73_Root_rep1'};
    my $h = $hoh{$key}{'B73_Root_rep2'};
    
    print $matrix_out_fh "$key\t$e\t$f\t$g\t$h\n";   
}
close $matrix_out_fh;
exit;
