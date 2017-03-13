#! /usr/bin/perl -w
# rnaseq_matrix.pl
# This script will make a matrix for the schmitz RNAseq dataset
# Jaclyn Noshay
# Jan_3_2017

use strict;
use Getopt::Long();
my $usage = "\nUsage: $0 --matrix_out <output file for matrix> --transcripts <list of transcripts> --help <help on usage of script>\n\n";
my ($matrix_out, $transcripts, $help);

Getopt::Long::GetOptions('matrix_out=s' => \$matrix_out,
			 'transcripts=s' => \$transcripts,
			 'h|help' => \$help
			);

if (defined($help)) {
    die $usage;
}

# Open genes and make hash of genes to be filled by files
open (my $transcripts_fh, '<', $transcripts) or die "Can't open file $transcripts\n\n";
my %hoh;

while (my $line = <$transcripts_fh>) {
    chomp $line;
    $line =~ s/>//g;
    $hoh{$line}{'207_Leaf_rep1'} = 'N/A';
    $hoh{$line}{'207_Leaf_rep2'} = 'N/A';
    $hoh{$line}{'207_Root_rep1'} = 'N/A';
    $hoh{$line}{'207_Root_rep2'} = 'N/A';
    $hoh{$line}{'B73_Leaf_rep1'} = 'N/A';
    $hoh{$line}{'B73_Leaf_rep2'} = 'N/A';  
    $hoh{$line}{'B73_Root_rep1'} = 'N/A';  
    $hoh{$line}{'B73_Root_rep2'} = 'N/A';  
    $hoh{$line}{'PHJ89_Leaf_rep1'} = 'N/A';  
    $hoh{$line}{'W22_Leaf_rep1'} = 'N/A';  
    $hoh{$line}{'W22_Leaf_rep2'} = 'N/A';  
    $hoh{$line}{'W22_Root_rep1'} = 'N/A';    
}

# Hard code in the files want to include in the matrix
# Could also do this by putting all transcript files into one directory and use open dir
my @files = qw | /scratch.global/nosha003/schmitz/rnaseq/align/htseq/207_Leaf_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/207_Leaf_rep2.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/207_Root_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/207_Root_rep2.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/B73_Leaf_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/B73_Leaf_rep2.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/B73_Root_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/B73_Root_rep2.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/PHJ89_Leaf_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/W22_Leaf_rep1.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/W22_Leaf_rep2.htseqcount
		 /scratch.global/nosha003/schmitz/rnaseq/align/htseq/W22_Root_rep1.htseqcount
	       |;

my $file_fh;

# Open output file
open (my $matrix_out_fh, '>', $matrix_out) or die "\nCan't open file $matrix_out\n\n";

# Print header line (hard coded for this example)
print $matrix_out_fh "transcript";

foreach my $file (@files) {
    if ($file =~ /\/scratch.global\/nosha003\/schmitz\/rnaseq\/align\/htseq\/(.+)\.htseqcount/) {
	my $name = $1;
	print "$name\n";
	print $matrix_out_fh "\t$name";

	open ($file_fh, '<', $file) or die "Can't open file $file\n\n";
	#my $header = <$file_fh>;
	while (my $line = <$file_fh>) {
	    chomp $line;
	    my ($gene, $counts) = split ("\t", $line);
	    if (exists $hoh{$gene}) {
		$hoh{$gene}{$name} = $counts;
	      }
	  }
      }
  }
print $matrix_out_fh "\n";

# Fill in matrix
foreach my $key (sort {$a cmp $b} keys %hoh) {
    my $a = $hoh{$key}{'207_Leaf_rep1'};
    my $b = $hoh{$key}{'207_Leaf_rep2'};
    my $c = $hoh{$key}{'207_Root_rep1'};
    my $d = $hoh{$key}{'207_Root_rep2'};
    my $e = $hoh{$key}{'B73_Leaf_rep1'};
    my $f = $hoh{$key}{'B73_Leaf_rep2'};
    my $g = $hoh{$key}{'B73_Root_rep1'};
    my $h = $hoh{$key}{'B73_Root_rep2'};
    my $i = $hoh{$key}{'PHJ89_Leaf_rep1'};
    my $j = $hoh{$key}{'W22_Leaf_rep1'};
    my $k = $hoh{$key}{'W22_Leaf_rep2'};
    my $l = $hoh{$key}{'W22_Root_rep1'};
    
    print $matrix_out_fh "$key\t$a\t$b\t$c\t$d\t$e\t$f\t$g\t$h\t$i\t$j\t$k\t$l\n";   
}
close $matrix_out_fh;
exit;
