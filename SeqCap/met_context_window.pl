#!/bin/perl
#source Qing
#updated by pete 2017-04-17
#output is a proper zero based BED file per context
#each context has the same windows, so long as one of the context has data
#windows without data in any context are omited

my $INFILE = $ARGV[0];
my $WINDOW = $ARGV[1] || 100;

# perl ./met_context_window_lin.pl example_input.txt 100

my %dat;
my $CGOUTFILE = "$INFILE.$WINDOW.CG.bed";
my $CHGOUTFILE = "$INFILE.$WINDOW.CHG.bed";
my $CHHOUTFILE = "$INFILE.$WINDOW.CHH.bed";

# Read in data
open (IN, $INFILE) or die;
open (CGOUT, ">$CGOUTFILE") or die;
open (CHGOUT, ">$CHGOUTFILE") or die;
open (CHHOUT, ">$CHHOUTFILE") or die;

print STDERR "Reading in report file...\n";
my $ts = localtime();
print STDERR $ts, "\n";

while (my $l = <IN>) {

	chomp($l);
	my @f = split(/\t/, $l);

	#udated column numbers for my file which has a new col 2
	# note this is using CT rather than eff_CT
	my $chr = $f[0];
	my $pos = $f[2];
	my $ccnt = $f[7];
	my $ctcnt = $f[8];
	my $context = $f[4];

	my $bin = int(($pos-1 )/ $WINDOW);

	#$dat{$chr}->{$bin}->{'total'}++;
	$dat{$chr}->{$bin}->{$context}->{'C'}+=$ccnt;
	$dat{$chr}->{$bin}->{$context}->{'CT'}+=$ctcnt;
	#count sites
	$dat{$chr}->{$bin}->{$context}->{'sites'}++;
	#print STDERR "$chr\t$bin\t$met\n";

}

print STDERR "Printing out window file...\n";
$ts = localtime();
print STDERR $ts, "\n";

# Emit window BED file (base zero coordinates)
# to get one based changed ($W * $WINDOW) to (($W * $WINDOW) + 1)
print STDERR "Printout time\n";
for my $C (sort {$a <=> $b} keys %dat) {

	for my $W (sort {$a <=> $b} keys %{ $dat{$C} }) {
		# CG output
		if($dat{$C}->{$W}->{'CG'}->{'CT'}!=0){
			my $cgfrac = ($dat{$C}->{$W}->{'CG'}->{'C'} / $dat{$C}->{$W}->{'CG'}->{'CT'});
			print CGOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CG'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CG'}->{'C'}, "\t", $dat{$C}->{$W}->{'CG'}->{'CT'}, "\t", $cgfrac, "\n";
		}else{
			print CGOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CG'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CG'}->{'C'}, "\t", $dat{$C}->{$W}->{'CG'}->{'CT'}, "\tNA\n";
		}
		# CHG output
		if($dat{$C}->{$W}->{'CHG'}->{'CT'}!=0){
			my $chgfrac = ($dat{$C}->{$W}->{'CHG'}->{'C'} / $dat{$C}->{$W}->{'CHG'}->{'CT'});
			print CHGOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CHG'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CHG'}->{'C'}, "\t", $dat{$C}->{$W}->{'CHG'}->{'CT'}, "\t", $chgfrac, "\n";
		}else{
			print CHGOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CHG'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CHG'}->{'C'}, "\t", $dat{$C}->{$W}->{'CHG'}->{'CT'}, "\tNA\n";
		}

		# CHH output
		if($dat{$C}->{$W}->{'CHH'}->{'CT'}!=0){
			my $chhfrac = ($dat{$C}->{$W}->{'CHH'}->{'C'} / $dat{$C}->{$W}->{'CHH'}->{'CT'});
			print CHHOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CHH'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CHH'}->{'C'}, "\t", $dat{$C}->{$W}->{'CHH'}->{'CT'}, "\t", $chhfrac, "\n";
		}else{
			print CHHOUT $C, "\t", ($W * $WINDOW), "\t", (($W + 1) * $WINDOW) , "\t", $dat{$C}->{$W}->{'CHH'}->{'sites'}, "\t", $dat{$C}->{$W}->{'CHH'}->{'C'}, "\t", $dat{$C}->{$W}->{'CHH'}->{'CT'}, "\tNA\n";
		}

	}

}

close CGOUT;
close CHGOUT;
close CHHOUT;
close IN;
