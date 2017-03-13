#! /usr/bin/perl -w
# deseq_compile.pl
# This script will take multiple deseq output files and compile into single matrix file with log2fc
# Jaclyn Noshay
# Jan_5_2017

use strict;
use Getopt::Long();
my $usage = "\nUsage: $0 --out <output file for compiled deseq data> -in <transcript file> --help <help on usage of script>\n\n";
my ($in, $out, $help);

Getopt::Long::GetOptions('out=s' => \$out,
			 'in=s' => \$in,
			 'h|help' => \$help
			);

if (defined($help)) {
  die $usage;
}

# Open genes and make hash of genes to be filled by files
open (my $in_fh, '<', $in) or die "Can't open file $in\n\n";
my %hoh;

while (my $line = <$in_fh>) {
    chomp $line;
    $line =~ s/>//g;
    $hoh{$line}{'v4leaf'} = 'N/A';  
    $hoh{$line}{'v4root'} = 'N/A';
    $hoh{$line}{'v4B_207leaf'} = 'N/A';
    $hoh{$line}{'v4B_207root'} = 'N/A';
    $hoh{$line}{'v4B_PHJleaf'} = 'N/A';
    $hoh{$line}{'v4B_W22leaf'} = 'N/A';
    $hoh{$line}{'v4B_W22root'} = 'N/A';
    $hoh{$line}{'v4Bleaf_Broot'} = 'N/A';
}

# Hard code in the files want to include in the matrix
# Could also do this by putting all transcript files into one directory and use open dir
my @files = qw | /scratch.global/nosha003/schmitz/rnaseq/deseq/v4leafresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4rootresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4B_207leafresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4B_207rootresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4B_PHJleafresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4B_W22leafresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4B_W22rootresDESeq.txt
		 /scratch.global/nosha003/schmitz/rnaseq/deseq/v4Bleaf_BrootresDESeq.txt
	       |;

my $file_fh;

# Open output file
open (my $out_fh, '>', $out) or die "\nCan't open file $out\n\n";

# Print header line (hard coded for this example)
print $out_fh "transcript";

foreach my $file (@files) {
    if ($file =~ /\/scratch.global\/nosha003\/schmitz\/rnaseq\/deseq\/(.+)\.txt/) {
	my $name = $1;
	$name =~ s/resDESeq//;
	print "$name\n";
	print $out_fh "\t$name";

	open ($file_fh, '<', $file) or die "Can't open file $file\n\n";
	#my $header = <$file_fh>;
	while (my $line = <$file_fh>) {
	    chomp $line;
	    my ($gene, $log2fc) = split ("\t", $line);
	    $gene =~ s/"//g;
	    if (exists $hoh{$gene}) {
		$hoh{$gene}{$name} = $log2fc;
	      }
	  }
      }
  }
print $out_fh "\n";

# Fill in matrix
foreach my $key (sort {$a cmp $b} keys %hoh) {
    my $a = $hoh{$key}{'v4leaf'};
    my $b = $hoh{$key}{'v4root'};
    my $c = $hoh{$key}{'v4B_207leaf'};
    my $d = $hoh{$key}{'v4B_207root'};
    my $e = $hoh{$key}{'v4B_PHJleaf'};
    my $f = $hoh{$key}{'v4B_W22leaf'};
    my $g = $hoh{$key}{'v4B_W22root'};
    my $h = $hoh{$key}{'v4Bleaf_Broot'};
    
    print $out_fh "$key\t$a\t$b\t$c\t$d\t$e\t$f\t$g\t$h\n";   
}
close $out_fh;
exit;
