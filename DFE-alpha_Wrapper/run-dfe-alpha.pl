#!/usr/bin/perl -w

use strict;
no warnings 'uninitialized';

#### Runs DFE-alpha

if (@ARGV != 2) {
	die "\nUsage: run-dfe-alpha.pl <tissue> <# of iterations>\n";
}

my $tissue = shift(@ARGV);
my $iter = shift(@ARGV);

my $file = my $out = "";
for (my $i = 1; $i <= $iter; $i++) {
	
	## est_dfe run 1: neutral estimation
	$file = $tissue."-config-neutral-".$i.".txt";
	system("/opt/dfe-alpha/est_dfe -c $file");
	
	## est_dfe run 2: selected estimation
	$file = $tissue."-config-selected-".$i.".txt";
	system("/opt/dfe-alpha/est_dfe -c $file");
	
	## est_alpha_omega, run 3: estimates alpha and omega
	$file = $tissue."-config-alpha-".$i.".txt";
	system("/opt/dfe-alpha/est_alpha_omega -c $file");
	
	## prop_muts_in_s_range, run 4: estimates the proportion of deleterious mutations
	$file = $tissue."_selected_results_".$i."/est_dfe.out";
	$out = $tissue."_prop_del_muts_".$i.".txt";
	system("/opt/dfe-alpha/prop_muts_in_s_ranges -c $file -o $out"); 
}