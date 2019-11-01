#!/usr/bin/perl -w

use strict;
no warnings 'uninitialized';

## In case things get disrupted repeatedly:
my $start = 1;

## Goes through all results and results directories from command-line DFE-alpha to pull out relevant information

if (@ARGV != 2) {
    die "\nUsage: parse-command-line-DFEalpha.pl <tissue> <# of iterations>\n";
}
my $tissue = shift(@ARGV);
my $num = shift(@ARGV);

## Set the path
my $path = `pwd`;
chomp $path; 		# strip terminal white space from the path

## Print the OUT as:
my $out = $tissue."_dfe_results_".$num.".txt";

open(OUT, "> $out") or die "Cannot write outfile: $!";
print OUT "tissue\tresample\tN1\tN2\tt2\tf0\tNw\tbeta\tEs\talpha\tomega_a\tlambda\tsel_diverg\tNeS_0..1_propmut\tNeS_1..10_propmut\tNeS_10..100_propmut\tNeS_100..inf_propmut\t\n";


my $file = my $line = "";
my %dfe = ();## Make as a HoA
my @line = ();
my $i = 0;

## For each of the resampled sets:
for $i ($start..$num) {

    $file = $path."/".$tissue."_selected_results_".$i."/est_dfe.out";
    open (FILE, '<', $file) or die "Cannot open $file: $!";
    print "I am opening $file now...\n";
    
    while ($line = <FILE>) {
	chomp($line);
	## split the line on white space
	@line = split(/\s+/, $line);
	#N1 100 N2 4 t2 22.8475 Nw 57.43 b 0.0500 Es -267910321460.079773 f0 0.976947500 L -1303.1421

	$dfe{$i}[0] = $line[1];
	$dfe{$i}[1] = $line[3];
	$dfe{$i}[2] = $line[5];
	$dfe{$i}[3] = $line[13];
	$dfe{$i}[4] = $line[7];
	$dfe{$i}[5] = $line[9];
	$dfe{$i}[6] = $line[11];
    }
    close(FILE);
    
    $file = $path."/".$tissue."_alpha_results_".$i.".txt";
    open (FILE, '<', $file) or die "Cannot open file: $!";
    print "I am opening $file now...\n";

	while ($line = <FILE>) {
	chomp($line);
	## split the line on white space
	@line = split(/\s+/, $line);
	#lambda 0.053327 selected_divergence 0.012742 alpha 1.000000 omega_A 0.238933

	$dfe{$i}[7] = $line[5];
	$dfe{$i}[8] = $line[7];
	$dfe{$i}[9] = $line[1];
	$dfe{$i}[10] = $line[3];
    }
    close(FILE);
    
    $file = $path."/".$tissue."_prop_del_muts_".$i.".txt";
    open (FILE, '<', $file) or die "Cannot open file: $!";
    print "I am opening $file now...\n";

    while ($line = <FILE>) {
	chomp($line);
	## split the line on white space
	@line = split(/\s+/, $line);
	#0.000000 1.000000 0.193756 1.000000 10.000000 0.023642 10.000000 100.000000 0.026526 100.000000 -99.000000 0.756076 
	$dfe{$i}[11] = $line[2];
	$dfe{$i}[12] = $line[5];
	$dfe{$i}[13] = $line[8];
	$dfe{$i}[14] = $line[11];
    }
    close(FILE);
}

foreach $i (sort keys %dfe) {
    print OUT "$tissue\t$i\t";
    for my $j (0..$#{ $dfe{$i} } ) {
	print OUT "$dfe{$i}[$j]\t";
    }
    print OUT "\n";
}
