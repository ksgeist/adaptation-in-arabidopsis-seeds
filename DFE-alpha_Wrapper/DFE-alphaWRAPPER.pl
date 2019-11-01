#!/usr/bin/perl -w

use strict;
no warnings 'uninitialized';
use IO::Handle;
use Getopt::Long;
use List::Util 'shuffle';

## Takes Polymorphorama output and creates Site Frequency Spectra (SFS) suitable for use
## with DFE-alpha
## Also generates the necessary divergence and config files

## Updating so that we now give counts from an alignment that is all population sequences
## against the reference genome OF THAT SPECIES
## (so we're not calculating interspecific divergence; we're calculating "divergence"
## between population sequences and it's reference genome)

#########################################################################################
## !!!!!!!!!!!!!
## Important: Assumes that Gene List 1 is the focal gene set, and Gene List 2 is background set if that's the resampling scheme you're using
#########################################################################################

## Initialize variables:
my $bs_w_repl;	
my $bs_wo_repl;				
my $labelshuffle;
my $bkgd;				## Indicates whether you are including a larger set as genelist2
my $summary;			## Summary file from Polymorphorama
my $rep_poly;			## Replacement polymorphism frequencies from Polymorphorama
my $syn_poly;			## Synonymous polymorphism from Polymorphorama
my $paml_div;			## PAML divergence file
my $glist1;				## List for gene set 1; assumes each gene ID is on its own line
my $glist2;				## List for gene set 2; assumes each gene ID is on its own line
my $num_iter;			## The number of times you want to resample (number of permutations)
my $set1;				## The name for the first gene set
my $set2;				## The name for the second gene set
my $run;				## Parameter to launch DFE-alpha

my %synsites = my %repsites = my %poly_syn_total = my %poly_rep_total = my %syn_div = my %rep_div = my %rep_poly = my %syn_poly = ();
my @line = my @glist1 = my @glist2 = ();
my $var = my $vector = my $sampsize = "";
my $gene = "";

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			SET OPTIONS					##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions(	"bs_w_repl" => \$bs_w_repl, # optional Boolean; set to false if not given
			"bs_wo_repl" => \$bs_wo_repl, # optional Boolean; set to false if not given
			"bkgd" => \$bkgd, # optional Boolean; set to false if not given
			"labelshuffle" => \$labelshuffle, # optional Boolean; set to false if not given
			"run" => \$run, #optional Boolean, false if not given
			"summary:s" => \$summary, # optional string parameter
			"rep_poly:s" => \$rep_poly, # optional string parameter
			"syn_poly:s" => \$syn_poly, # optional string parameter
			"paml_div:s" => \$paml_div, # optional string parameter
			"num_iter:i" => \$num_iter, # optional integer parameter
			"glist1:s" => \$glist1, # optional string parameter
			"glist2:s" => \$glist2, # optional string parameter
			"set1:s" => \$set1, # optional string parameter
			"set2:s" => \$set2) # optional string parameter
			or die "\n***WARNING: Error in command line arguments. I'm terminating.\n\n"; 

# N.B.: When an argument, is called only as a flag, its value is 0. When it is not called
# at all, it is null

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			SET ERRORS					##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

if($num_iter > 0) {
	if(! $bs_w_repl) {
		if(! $bs_wo_repl) {
			if(! $labelshuffle) {
				die "\n***WARNING: No resampling method was provided. I'm terminating.***\n\n"; }}}}
						
if (! $paml_div) {
	print "\nNo alternate divergence file given. Divergences from Polymorphorama will be used.\n\n" }

if (! $set1) {
	print "\nNo name for gene set 1 provided. A generic identifier will be used.\n\n"; 
	$set1 = "set1"; }

if (defined $glist2){
	if (! $set2) {
		print "\nNo name for gene set 2 provided. A generic identifier will be used.\n\n"; 
		$set2 = "set2"; }	
}

if (! defined $num_iter) {
	print "\nYou did not tell me how many times to resample. Only running point estimates.\n\n";
	$num_iter = 0; }



if (! $glist1) {
	die "\n***WARNING: You must give me at least one gene list. I'm terminating.***\n\n"; }


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			SUMMARY FILE				##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

open(SUMMARY, $summary) or die "\nCould not open the Polymorphorama Summary file: $!\n\n";
while (my $line = <SUMMARY>) {
	chomp $line;
	## Check for carriage returns and strip if needed:
	$line =~ s/\r|\n//g;

	@line = split(/\t/, $line);
	## Assigns to genes
	$poly_syn_total{$line[0]} = $line[4];
	$poly_rep_total{$line[0]} = $line[17];

	if(!defined $paml_div) {	## If no information from PAML was provided, use Polymorphorama divergence info:
		$synsites{$line[0]} = $line[3];		
		$repsites{$line[0]} = $line[16];
		$syn_div{$line[0]} = $line[6];
		$rep_div{$line[0]} = $line[19];
	}
}
close(SUMMARY);


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			REPLACEMENT POLY			##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

open(REP_POLY, $rep_poly) or die "\nCould not open the Polymorphorama Replacement Frequencies file: $!\n\n";
## Take the replacement polymorphism vectors from Polymorphorama 
## Store in a hash with gene as key and vector 1/xi --> xi/xi as value; SKIP the first
## number as it is divergence data

while (my $line = <REP_POLY>) {
	$vector = "";
	chomp $line;
	## Check for carriage returns and strip if needed:
	$line =~ s/\r|\n//g;

	@line = split(/\t/, $line);
	for (my $g = 2; $g <= $#line; $g++) {
		$vector = $vector . $line[$g] . "\t";
		$rep_poly{$line[0]} = $vector;
	}
	$sampsize = $#line;		
}
close(REP_POLY);


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			SYNONYMOUS POLY		    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

open(SYN_POLY, $syn_poly) or die "\nCould not open the Polymorphorama Neutral Frequencies file: $!\n\n";
## Do the same with the neutral polymorphism vectors from Polymorphorama
while (my $line = <SYN_POLY>) {
	$vector = "";
	chomp $line;
	## Check for carriage returns and strip if needed:
	$line =~ s/\r|\n//g;

	@line = split(/\t/, $line);
	
	for (my $b = 2; $b <= $#line; $b++) {		## Skips the first bin and second bins; bin 0 is the gene ID and bin 1 is the divergence calculated by Polymorphorama
		$vector = $vector . $line[$b] . "\t";
		$syn_poly{$line[0]} = $vector;
	}	
}
close(SYN_POLY);


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			PAML DIVERGENCE		    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##

if (defined $paml_div) {
	open(PAML, $paml_div) or die "\nCould not open the PAML divergence file: $!\n\n";
	while (my $line = <PAML>) {
		$vector = "";
		chomp $line;
		## Skip the header line or any line containing NAs:
		if ($line !~ /^gene/ & $line !~ /NA/) {
			@line = split(/\t/, $line);

		# 	gene	lnL	t	S	N	dN/dS	dN	dS	NG.dNdS	NG.dN	NG.dS
			$synsites{$line[0]} = $line[3];
			$repsites{$line[0]} = $line[4];
			$rep_div{$line[0]} = sprintf "%.0f", $line[6]*$line[4];
			$syn_div{$line[0]} = sprintf "%.0f", $line[7]*$line[3];

			$var = sprintf "%.0f", $line[6]*$line[4];
		}	
	}
	close(PAML);
}

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			GENE LIST 1			    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
my $gct1 = 0;
open (LIST1, $glist1) or die "\nCould not open the gene list for $set1: $!\n\n";
while (my $line = <LIST1>) {
	if ($line =~ /(^\S+)/) {
		push (@glist1, $1);
#   		print "$1\n";
		$gct1++;
	}
}
close(LIST1);


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			GENE LIST 2			    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
my $gct2 = 0;

if (defined $glist2) {
	open (LIST2, $glist2) or die "\nCould not open the gene list for $set2: $!\n\n";
	while (my $line = <LIST2>) {
		if ($line =~ /(^\S+)/) {
			push (@glist2, $1);
# 			print "$1\n";
			$gct2++;
		}
	}
	close(LIST2);
}

###### WARNING IF SECOND GENE LIST IS TOO SMALL:
if (defined $glist2) {
 	if ($gct2 < $gct1 && !defined $bkgd){		## If the second gene set is less than the first set
		print "Heads up! Are you sure you chose the correct resampling method? I have noticed that your $glist2 is much smaller than your $glist1.\n"
	}}

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##					RESAMPLING CODE	   				##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

my $i = 1;
my $rep_div_scalar = my $syn_div_scalar = my $rep_sites_scalar = my $syn_sites_scalar = my $g = my $roundrepsites = my $roundsynsites = my $rep0thbin = my $syn0thbin = my $poly_rep_total_scalar = my $poly_syn_total_scalar = 0;
my %resampledgenes = my %rep_sfs = my %syn_sfs = my %syngenes = my %repgenes = ();
my @rep_sfs = my @syn_sfs = my @resamplist = my @comblist = my @newglist1 = my @newglist2 = ();
my $dnds;
my $pnps;
my $alpha;
my %overlap = ();
my $newglist1ref = my $newglist2ref = my $resamp1_ref = my $resamp2_ref = my $temp_ref = "";
my $n1 = my $n2 = 0;


##########################################
##	   GENERATE A POINT ESTIMATE ONLY	##
##########################################
### Code for a -num_iter of 0 option. Generates the point estimate on one or two gene sets. 
if ($num_iter == 0) {
	$i = 0;
	####################################################
	## Proceed with making the SFSs:				  ##
	####################################################

		## Reset at the beginning of each iteration
		$rep_div_scalar = $syn_div_scalar = $rep_sites_scalar = $syn_sites_scalar = $roundrepsites = $roundsynsites = $rep0thbin = $syn0thbin = $poly_rep_total_scalar = $poly_rep_total_scalar = 0;
		%rep_sfs = %syn_sfs = %syngenes = %repgenes = ();
		@rep_sfs = @syn_sfs = @resamplist = ();
	
		## Create OUTFILES
		my $repout = $set1."-replacement-SFS-".$i.".txt";
		my $synout = $set1."-neutral-SFS-".$i.".txt";
		my $sfs = $set1."-SFS-".$i.".txt";
		my $divout = $set1."-divergence-".$i.".txt";
		my $config1 = $set1."-config-neutral-".$i.".txt";
		my $config2 = $set1."-config-selected-".$i.".txt";
		my $config3 = $set1."-config-alpha-".$i.".txt";
		my $dir1 = $set1."_neutral_results_".$i;
		my $dir2 = $set1."_selected_results_".$i;
		my $alphafile = $set1."_alpha_results_".$i.".txt";
		my $handalphafile = $set1."-handcalc-alpha-results-".$i.".txt";
		my $genesfile1 = $set1."-gene-list-".$i.".txt";
		
		my $repout2 = $set2."-replacement-SFS-".$i.".txt";
		my $synout2 = $set2."-neutral-SFS-".$i.".txt";
		my $sfs2 = $set2."-SFS-".$i.".txt";
		my $divout2 = $set2."-divergence-".$i.".txt";
		my $config1_2 = $set2."-config-neutral-".$i.".txt";
		my $config2_2 = $set2."-config-selected-".$i.".txt";
		my $config3_2 = $set2."-config-alpha-".$i.".txt";
		my $dir1_2 = $set2."_neutral_results_".$i;
		my $dir2_2 = $set2."_selected_results_".$i;
		my $alphafile2 = $set2."_alpha_results_".$i.".txt";
		my $handalphafile2 = $set2."-handcalc-alpha-results-".$i.".txt";
		my $genesfile2 = $set2."-gene-list-".$i.".txt";
		
		## Open OUTFILES:
		open(SFS, "> $sfs") or die "Cannot write outfile: $!";
		open(DIV, "> $divout") or die "Cannot write outfile: $!";
		open(CONFIG1, "> $config1") or die "Cannot write outfile: $!";
		open(CONFIG2, "> $config2") or die "Cannot write outfile: $!";
		open(CONFIG3, "> $config3") or die "Cannot write outfile: $!";
		open(HANDA, "> $handalphafile") or die "Cannot write outfile: $!";
		open(GENES1, "> $genesfile1") or die "Cannot write outfile: $!";

	#############################################
	## This samples across each of the resampled genes
		foreach my $k (0..$gct1) {
			$gene = $glist1[$k];
			print GENES1 "$gene\n";

			## We have to add the "divergences" to the 0th bin or it will be incorrect:
			## To create the 0th bin of the array, we take the total number of sites and subtract those contributing to both divergence and polymorphism (the unmutated sites):
			$rep0thbin = sprintf "%.0f", ($repsites{$gene}-($rep_div{$gene}+$poly_rep_total{$gene}));
			$syn0thbin = sprintf "%.0f", ($synsites{$gene}-($syn_div{$gene}+$poly_syn_total{$gene}));	
	# 		print "LIST 1: $gene\t$repsites{$gene}\t$roundrepsites\t$rep0thbin\t$rep_div{$gene}\t$poly_rep_total{$gene}\n";
	
			## Dual hashes because duplicates get overwritten:
			$g++;
			## The value is the full SFS, including that new 0th bin (unmutated sites) with the gene NUMBER as key
			$rep_sfs{$g} = $rep0thbin."\t$rep_poly{$gene}";		## Gene NUMBER as key, with gene IDs as values
			$repgenes{$g} = $gene;
	# 		print "LIST 1: $gene\t$rep_sfs{$repgenes{$g}}\n";
			$syn_sfs{$g} = $syn0thbin."\t$syn_poly{$gene}";
			$syngenes{$g} = $gene;

			## Sum the sites and divergences:
			## Rounds the number of total replacement and polymorphism sites for that gene	
			$roundrepsites = sprintf "%.0f", $repsites{$gene};
			$roundsynsites = sprintf "%.0f", $synsites{$gene};
			$rep_sites_scalar += $roundrepsites;
			$syn_sites_scalar += $roundsynsites;
			$rep_div_scalar += $rep_div{$gene};						## Dn
			$syn_div_scalar += $syn_div{$gene};						## Ds
			$poly_rep_total_scalar += $poly_rep_total{$gene};		## Pn
			$poly_syn_total_scalar += $poly_syn_total{$gene};		## Ps
		}

		#############################################
		## Now we just need to sum the SFS for DFE-alpha, plus print our various out files
		foreach $gene (keys %rep_sfs) {
			my @sumrep = split(/\t/, $rep_sfs{$gene});
			my @sumsyn = split(/\t/, $syn_sfs{$gene});
		
			foreach my $i (0..$#sumrep) {
				$rep_sfs[$i] += $sumrep[$i];
				$syn_sfs[$i] += $sumsyn[$i];
			}
		}

		## Series of print statements
		print HANDA "Dn\tDs\tDn/Ds\tPn\tPs\tPn/Ps\talpha\tN\n";
		if ($syn_div_scalar > 0) {
			$dnds = $rep_div_scalar/$syn_div_scalar;
		}
		else {
			$dnds = "NA";
		}

		if ($poly_rep_total_scalar > 0) {
			if ($poly_syn_total_scalar > 0) {
				$pnps = $poly_rep_total_scalar/$poly_syn_total_scalar;
			}
		}
		else {
			$pnps = "NA";
		}

		if ($dnds ne "NA" && $pnps ne "NA"){
			if ($dnds > 0) {
				$alpha = 1-(($pnps)/($dnds));
			}
			else {
				$alpha = "NA";
			}
		}
		else {
			$alpha = "NA";
		}
		print HANDA "$rep_div_scalar\t$syn_div_scalar\t$dnds\t$poly_rep_total_scalar\t$poly_syn_total_scalar\t$pnps\t$alpha\t$n1\n";
	
		print SFS "1\n";
		print SFS $sampsize."\n";
		print SFS join "\t", @rep_sfs;
		print SFS "\t$rep_div_scalar\n";
	
		print SFS join "\t", @syn_sfs;
		print SFS "\t$syn_div_scalar\n";
	
		print DIV "1\t".int($rep_sites_scalar)."\t".$rep_div_scalar."\n";
		print DIV "0\t".int($syn_sites_scalar)."\t".$syn_div_scalar."\n";

		print CONFIG1 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG1 "sfs_input_file\t$sfs\n";
		print CONFIG1 "est_dfe_results_dir\t".$dir1."\n";
		print CONFIG1 "site_class\t0\n";
		print CONFIG1 "fold\t1\n";
		print CONFIG1 "epochs\t2\n";
		print CONFIG1 "search_n2\t1\n";
		print CONFIG1 "t2_variable\t1\n";
		print CONFIG1 "t2\t50\n";
	
		print CONFIG2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG2 "sfs_input_file\t$sfs\n";
		print CONFIG2 "est_dfe_results_dir\t".$dir2."\n";
		print CONFIG2 "est_dfe_demography_results_file\t".$dir1."/est_dfe.out\n";
		print CONFIG2 "site_class\t1\n";
		print CONFIG2 "fold\t1\n";
		print CONFIG2 "epochs\t2\n";
		print CONFIG2 "mean_s_variable\t1\n";
		print CONFIG2 "mean_s\t-0.1\n";
		print CONFIG2 "beta_variable\t1\n";
		print CONFIG2 "beta\t0.5\n";
	
		print CONFIG3 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG3 "divergence_file\t".$divout."\n";
		print CONFIG3 "est_alpha_omega_results_file\t".$alphafile."\n";
		print CONFIG3 "est_dfe_results_file\t".$dir2."/est_dfe.out\n";
		print CONFIG3 "neut_egf_file\t".$dir1."/neut_egf.out\n";
		print CONFIG3 "sel_egf_file\t".$dir2."/sel_egf.out\n";
		print CONFIG3 "do_jukes_cantor\t1\n";
		print CONFIG3 "remove_poly\t0\n";

		## Reset at the beginning of each iteration
		$rep_div_scalar = $syn_div_scalar = $rep_sites_scalar = $syn_sites_scalar = $roundrepsites = $roundsynsites = $rep0thbin = $syn0thbin = $poly_rep_total_scalar = $poly_rep_total_scalar = 0;
		%rep_sfs = %syn_sfs = %syngenes = %repgenes = ();
		@rep_sfs = @syn_sfs = ();

	##############################################################
	## Now we need to do all this again for the second gene set:##
	##############################################################
	## Open OUTFILES		
	
	if(defined $glist2) {
		open(SFS2, "> $sfs2") or die "Cannot write outfile: $!";
		open(DIV2, "> $divout2") or die "Cannot write outfile: $!";
		open(CONFIG1_2, "> $config1_2") or die "Cannot write outfile: $!";
		open(CONFIG2_2, "> $config2_2") or die "Cannot write outfile: $!";
		open(CONFIG3_2, "> $config3_2") or die "Cannot write outfile: $!";
		open(HANDA2, "> $handalphafile2") or die "Cannot write outfile: $!";	
		open(GENES2, "> $genesfile2") or die "Cannot write outfile: $!";

		foreach my $k (0..$gct2) {
			$gene = $glist2[$k];
			print GENES2 "$gene\n";
		
			## To create the 0th bin of the array, we take the total number of sites and subtract those contributing to both divergence and polymorphism (the unmutated sites):
			$rep0thbin = sprintf "%.0f", ($repsites{$gene}-($rep_div{$gene}+$poly_rep_total{$gene}));
			$syn0thbin = sprintf "%.0f", ($synsites{$gene}-($syn_div{$gene}+$poly_syn_total{$gene}));	
#			print "LIST 2: $gene\t$repsites{$gene}\t$roundrepsites\t$rep0thbin\t$rep_div{$gene}\t$poly_rep_total{$gene}\n";
			## Sum all of the polymorphic sites for a given gene
	#		$rep0thbin = sprintf "%.0f", $repsites{$gene}-$poly_rep_total{$gene};
	#		$syn0thbin = sprintf "%.0f", $synsites{$gene}-$poly_syn_total{$gene};
			
			## Dual hashes because duplicates get overwritten:
			$g++;
			$rep_sfs{$g} = $rep0thbin."\t$rep_poly{$gene}";
			$repgenes{$g} = $gene;
			$syn_sfs{$g} = $syn0thbin."\t$syn_poly{$gene}";
			$syngenes{$g} = $gene;

			## Sum the sites and divergences:
			$roundrepsites = sprintf "%.0f", $repsites{$gene};
			$roundsynsites = sprintf "%.0f", $synsites{$gene};
			$rep_sites_scalar += $roundrepsites;		
			$syn_sites_scalar += $roundsynsites;
			$rep_div_scalar += $rep_div{$gene};						## Dn
			$syn_div_scalar += $syn_div{$gene};						## Ds
			$poly_rep_total_scalar += $poly_rep_total{$gene};		## Pn
			$poly_syn_total_scalar += $poly_syn_total{$gene};		## Ps
		}
		foreach $gene (keys %rep_sfs) {
			my @sumrep = split(/\t/, $rep_sfs{$gene});
			my @sumsyn = split(/\t/, $syn_sfs{$gene});

			foreach my $i (0..$#sumrep) {
				$rep_sfs[$i] += $sumrep[$i];
				$syn_sfs[$i] += $sumsyn[$i];
			}
		}
		
		print HANDA2 "Dn\tDs\tDn/Ds\tPn\tPs\tPn/Ps\talpha\tN\n";
		if ($syn_div_scalar > 0) {
			$dnds = $rep_div_scalar/$syn_div_scalar;
		}
		else {
			$dnds = "NA";
		}

		if ($poly_rep_total_scalar > 0) {
				if ($poly_syn_total_scalar > 0) {
				$pnps = $poly_rep_total_scalar/$poly_syn_total_scalar;
			}
		}
		else {
			$pnps = "NA";
		}

		if ($dnds ne "NA" && $pnps ne "NA"){
			if ($dnds > 0) {
				$alpha = 1-(($pnps)/($dnds));
			}
			else {
				$alpha = "NA";
			}
		}
		else {
			$alpha = "NA";
		}

		## Series of print statements
		print HANDA2 "$rep_div_scalar\t$syn_div_scalar\t$dnds\t$poly_rep_total_scalar\t$poly_syn_total_scalar\t$pnps\t$alpha\t$n2\n";
	
		print SFS2 "1\n";
		print SFS2 $sampsize."\n";
		print SFS2 join "\t", @rep_sfs;
		## Add an n+1-th bin of 0s -- this is xi/xi... but with Polymorphorama we do not have that information. So, we use 0 instead.
	# 	print SFS2 "\t0\n";
	## Newest update: rather than 0 in the last bin, we're doing divergence instead (with the rationale being that is every site differs [xi/xi bin] then isn't this divergence? the only other option would be it's the total differences seen relative to the REFERENCE genome... and i still can't suss out what they mean exactly!)
		print SFS2 "\t$rep_div_scalar\n";
		print SFS2 join "\t", @syn_sfs;
		print SFS2 "\t$syn_div_scalar\n";
		print DIV2 "1\t".int($rep_sites_scalar)."\t".$rep_div_scalar."\n";
		print DIV2 "0\t".int($syn_sites_scalar)."\t".$syn_div_scalar."\n";
		print CONFIG1_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG1_2 "sfs_input_file\t$sfs2\n";
		print CONFIG1_2 "est_dfe_results_dir\t".$dir1_2."\n";
		print CONFIG1_2 "site_class\t0\n";
		print CONFIG1_2 "fold\t1\n";
		print CONFIG1_2 "epochs\t2\n";
		print CONFIG1_2 "search_n2\t1\n";
		print CONFIG1_2 "t2_variable\t1\n";
		print CONFIG1_2 "t2\t50\n";
		print CONFIG2_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG2_2 "sfs_input_file\t$sfs2\n";
		print CONFIG2_2 "est_dfe_results_dir\t".$dir2_2."\n";
		print CONFIG2_2 "est_dfe_demography_results_file\t".$dir1_2."/est_dfe.out\n";
		print CONFIG2_2 "site_class\t1\n";
		print CONFIG2_2 "fold\t1\n";
		print CONFIG2_2 "epochs\t2\n";
		print CONFIG2_2 "mean_s_variable\t1\n";
		print CONFIG2_2 "mean_s\t-0.1\n";
		print CONFIG2_2 "beta_variable\t1\n";
		print CONFIG2_2 "beta\t0.5\n";
		print CONFIG3_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG3_2 "divergence_file\t".$divout2."\n";
		print CONFIG3_2 "est_alpha_omega_results_file\t".$alphafile2."\n";
		print CONFIG3_2 "est_dfe_results_file\t".$dir2_2."/est_dfe.out\n";
		print CONFIG3_2 "neut_egf_file\t".$dir1_2."/neut_egf.out\n";
		print CONFIG3_2 "sel_egf_file\t".$dir2_2."/sel_egf.out\n";
		print CONFIG3_2 "do_jukes_cantor\t1\n";
		print CONFIG3_2 "remove_poly\t0\n";		
}}


####################################################
## Resample according to specified scheme:		  ##
####################################################

## Iterative Resampling to create some number of SFS's:
for $i (1..$num_iter) {

	$n1 = $gct1; 
	$n2 = $gct2;
			
	##########################
	##	   BOOTSTRAPPING	##
	##########################

	if ($bs_w_repl == 1) {
			print "We are bootstraptrapping the genes in $glist1 WITH replacement...\n";
	 		($resamp1_ref) = resample_with(@glist1);
	 		
	}
	if ($bs_wo_repl == 1) {
			print "We are bootstraptrapping $n1 genes from $glist2 WITHOUT replacement...\n";
	 		($resamp1_ref) = resample_without($n1, @glist2);
	}
	
	##########################
	##	   LABEL SHUFFLING	##
	##########################

	if ($labelshuffle == 1) {
			if (! defined $glist2) {
				die "\n***WARNING: You must give me two gene lists if you want to shuffle 	
				labels. I'm terminating.***\n\n"; }

			##########################################
			##	   LABEL SHUFFLING	NO SUBSAMPLING	##
			##########################################			
			if(! defined $bkgd) {
			## Here we want to simply take the two lists, reshuffle them,
			## and pull out to originally-sized sets
			
				## We need to make sure they aren't overlapping
				($newglist1ref, $newglist2ref) = overlap(\@glist1, \@glist2);
				
				## Assigns the new gene lists, even if there was no overlap
				@newglist1 = @{$newglist1ref};
				@newglist2 = @{$newglist2ref};
				
				## Assigns new numbers of genes, even if there was no overlap
				$n1 = scalar(@newglist1);
				$n2 = scalar(@newglist2);
				
				print "We are shuffling the genes in $glist1 & $glist2 create reshuffled sets of sizes $n1 and $n2, respectively...\n\n";

				## Combines the two gene lists together into a bigger list
 				@comblist = (@newglist1, @newglist2);
 				## Performs the label shuffling on the two lists		
 				($resamp1_ref, $resamp2_ref) = label_shuffle($n1, $n2, @comblist);
 			}

			##################################################
			##	   LABEL SHUFFLING	SUBSAMPLING BACKGROUND	##
			##################################################	
			
			if($bkgd == 1) {
			## We sample a list1 sized set from list2, then shuffle their labels
			## Before we can combine them, make sure they're not overlapping sets:
			($newglist1ref, $newglist2ref) = overlap(\@glist1, \@glist2);
			
			## Assigns the new gene lists, even if there was no overlap
			@newglist1 = @{$newglist1ref};
			@newglist2 = @{$newglist2ref};
			
			## Assigns new numbers of genes, even if there was no overlap
			$n1 = scalar(@newglist1);
			$n2 = scalar(@newglist2);

			print "There are now $n1 and $n2 genes in the lists instead of $gct1 and $gct2, respectively.\n\n";

			print "We are shuffling the genes in $glist1 & $glist2 after sampling $n1 genes without replacement from $glist2...\n";
			
			#############
			## Step 1: ## Resample the larger list (list #2) without replacement to get something the same size as the focal list #1
			#############
			$temp_ref = resample_without($n1, @newglist2);

			#############
			## Step 2: 	## Turn the reference indices back into genes
			#############

			foreach my $k (0..$#$temp_ref) { 
# 				print "$$temp_ref[$k]\n";
				push(@resamplist, $$temp_ref[$k]); 
			}

			## Reset $n2 as the new size:
			$n2 = scalar(@resamplist);

			## ERROR FLAG:
				if ($n2 != $n1) {
					die "\n\n*** WARNING: There is an issue with your gene lists during resampling for label shuffling. Terminating.***\n";
				}

			#############
			## Step 3: ## Make a combined list:
			#############
			@comblist = (@newglist1, @resamplist);		## Combines the two gene lists together into a bigger list

			#############
			## Step 4: ## Send them through the label shuffling subroutine:
			#############
			($resamp1_ref, $resamp2_ref) = label_shuffle($n1, $n2, @comblist);
		}
	} 

	####################################################
	## Proceed with making the SFSs:				  ##
	####################################################

		## Reset at the beginning of each iteration
		$rep_div_scalar = $syn_div_scalar = $rep_sites_scalar = $syn_sites_scalar = $roundrepsites = $roundsynsites = $rep0thbin = $syn0thbin = $poly_rep_total_scalar = $poly_rep_total_scalar = 0;
		%rep_sfs = %syn_sfs = %syngenes = %repgenes = ();
		@rep_sfs = @syn_sfs = @resamplist = ();
	
		## Create OUTFILES
		my $repout = $set1."-replacement-SFS-".$i.".txt";
		my $synout = $set1."-neutral-SFS-".$i.".txt";
		my $sfs = $set1."-SFS-".$i.".txt";
		my $divout = $set1."-divergence-".$i.".txt";
		my $config1 = $set1."-config-neutral-".$i.".txt";
		my $config2 = $set1."-config-selected-".$i.".txt";
		my $config3 = $set1."-config-alpha-".$i.".txt";
		my $dir1 = $set1."_neutral_results_".$i;
		my $dir2 = $set1."_selected_results_".$i;
		my $alphafile = $set1."_alpha_results_".$i.".txt";
		my $handalphafile = $set1."-handcalc-alpha-results-".$i.".txt";
		my $genesfile1 = $set1."-gene-list-".$i.".txt";
	
		my $repout2 = $set2."-replacement-SFS-".$i.".txt";
		my $synout2 = $set2."-neutral-SFS-".$i.".txt";
		my $sfs2 = $set2."-SFS-".$i.".txt";
		my $divout2 = $set2."-divergence-".$i.".txt";
		my $config1_2 = $set2."-config-neutral-".$i.".txt";
		my $config2_2 = $set2."-config-selected-".$i.".txt";
		my $config3_2 = $set2."-config-alpha-".$i.".txt";
		my $dir1_2 = $set2."_neutral_results_".$i;
		my $dir2_2 = $set2."_selected_results_".$i;
		my $alphafile2 = $set2."_alpha_results_".$i.".txt";
		my $handalphafile2 = $set2."-handcalc-alpha-results-".$i.".txt";
		my $genesfile2 = $set2."-gene-list-".$i.".txt";

		## Open OUTFILES:
		open(SFS, "> $sfs") or die "Cannot write outfile: $!";
		open(DIV, "> $divout") or die "Cannot write outfile: $!";
		open(CONFIG1, "> $config1") or die "Cannot write outfile: $!";
		open(CONFIG2, "> $config2") or die "Cannot write outfile: $!";
		open(CONFIG3, "> $config3") or die "Cannot write outfile: $!";
		open(HANDA, "> $handalphafile") or die "Cannot write outfile: $!";
		open(GENES1, "> $genesfile1") or die "Cannot write outfile: $!";

	#############################################
	## This samples across each of the resampled genes
		foreach my $k (0..$#$resamp1_ref) {
			$gene = $$resamp1_ref[$k];
			print GENES1 "$gene\n";

			## We have to add the "divergences" to the 0th bin or it will be incorrect:
			## To create the 0th bin of the array, we take the total number of sites and subtract those contributing to both divergence and polymorphism (the unmutated sites):
			$rep0thbin = sprintf "%.0f", ($repsites{$gene}-($rep_div{$gene}+$poly_rep_total{$gene}));
			$syn0thbin = sprintf "%.0f", ($synsites{$gene}-($syn_div{$gene}+$poly_syn_total{$gene}));	
	# 		print "LIST 1: $gene\t$repsites{$gene}\t$roundrepsites\t$rep0thbin\t$rep_div{$gene}\t$poly_rep_total{$gene}\n";
	
			## Dual hashes because duplicates get overwritten:
			$g++;
			## The value is the full SFS, including that new 0th bin (unmutated sites) with the gene NUMBER as key
			$rep_sfs{$g} = $rep0thbin."\t$rep_poly{$gene}";		## Gene NUMBER as key, with gene IDs as values
			$repgenes{$g} = $gene;
	# 		print "LIST 1: $gene\t$rep_sfs{$repgenes{$g}}\n";
			$syn_sfs{$g} = $syn0thbin."\t$syn_poly{$gene}";
			$syngenes{$g} = $gene;

			## Sum the sites and divergences:
			## Rounds the number of total replacement and polymorphism sites for that gene	
			$roundrepsites = sprintf "%.0f", $repsites{$gene};
			$roundsynsites = sprintf "%.0f", $synsites{$gene};
			$rep_sites_scalar += $roundrepsites;
			$syn_sites_scalar += $roundsynsites;
			$rep_div_scalar += $rep_div{$gene};						## Dn
			$syn_div_scalar += $syn_div{$gene};						## Ds
			$poly_rep_total_scalar += $poly_rep_total{$gene};		## Pn
			$poly_syn_total_scalar += $poly_syn_total{$gene};		## Ps
		}

		#############################################
		## Now we just need to sum the SFS for DFE-alpha, plus print our various out files
		foreach $gene (keys %rep_sfs) {
			my @sumrep = split(/\t/, $rep_sfs{$gene});
			my @sumsyn = split(/\t/, $syn_sfs{$gene});
		
# 			Add them to the resampled genes hash, so we can skip them later:
			$resampledgenes{$gene}++;
		
			foreach my $i (0..$#sumrep) {
				$rep_sfs[$i] += $sumrep[$i];
				$syn_sfs[$i] += $sumsyn[$i];
			}
		}

		## Series of print statements
		print HANDA "Dn\tDs\tDn/Ds\tPn\tPs\tPn/Ps\talpha\tN\n";
		if ($syn_div_scalar > 0) {
			$dnds = $rep_div_scalar/$syn_div_scalar;
		}
		else {
			$dnds = "NA";
		}

		if ($poly_rep_total_scalar > 0) {
			if ($poly_syn_total_scalar > 0) {
				$pnps = $poly_rep_total_scalar/$poly_syn_total_scalar;
			}
		}
		else {
			$pnps = "NA";
		}

		if ($dnds ne "NA" && $pnps ne "NA"){
			if ($dnds > 0) {
				$alpha = 1-(($pnps)/($dnds));
			}
			else {
				$alpha = "NA";
			}
		}
		else {
			$alpha = "NA";
		}
		print HANDA "$rep_div_scalar\t$syn_div_scalar\t$dnds\t$poly_rep_total_scalar\t$poly_syn_total_scalar\t$pnps\t$alpha\t$n1\n";
	
		print SFS "1\n";
		print SFS $sampsize."\n";
		print SFS join "\t", @rep_sfs;
		## Add an n+1-th bin of 0s -- this is xi/xi... but with Polymorphorama we do not have that information. So, we use 0 instead.
	# 	print SFS "\t0\n";
	# 	print SFS "\n";	## Just new line, no extra 0!	
	## Newest update: rather than 0 in the last bin, we're doing divergence instead (with the rationale being that is every site differs [xi/xi bin] then isn't this divergence? the only other option would be it's the total differences seen relative to the REFERENCE genome... and i still can't suss out what they mean exactly!)
		print SFS "\t$rep_div_scalar\n";
	
		print SFS join "\t", @syn_sfs;
		## Add an n+1-th bin of 0s -- this is xi/xi... but with Polymorphorama we do not have that information. So, we use 0 instead.
	# 	print SFS "\t0\n";	
	# 	print SFS "\n";	## Just new line, no extra 0!	
		print SFS "\t$syn_div_scalar\n";
	
		print DIV "1\t".int($rep_sites_scalar)."\t".$rep_div_scalar."\n";
		print DIV "0\t".int($syn_sites_scalar)."\t".$syn_div_scalar."\n";

		print CONFIG1 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG1 "sfs_input_file\t$sfs\n";
		print CONFIG1 "est_dfe_results_dir\t".$dir1."\n";
		print CONFIG1 "site_class\t0\n";
		print CONFIG1 "fold\t1\n";
		print CONFIG1 "epochs\t2\n";
		print CONFIG1 "search_n2\t1\n";
		print CONFIG1 "t2_variable\t1\n";
		print CONFIG1 "t2\t50\n";
	
		print CONFIG2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG2 "sfs_input_file\t$sfs\n";
		print CONFIG2 "est_dfe_results_dir\t".$dir2."\n";
		print CONFIG2 "est_dfe_demography_results_file\t".$dir1."/est_dfe.out\n";
		print CONFIG2 "site_class\t1\n";
		print CONFIG2 "fold\t1\n";
		print CONFIG2 "epochs\t2\n";
		print CONFIG2 "mean_s_variable\t1\n";
		print CONFIG2 "mean_s\t-0.1\n";
		print CONFIG2 "beta_variable\t1\n";
		print CONFIG2 "beta\t0.5\n";
	
		print CONFIG3 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG3 "divergence_file\t".$divout."\n";
		print CONFIG3 "est_alpha_omega_results_file\t".$alphafile."\n";
		print CONFIG3 "est_dfe_results_file\t".$dir2."/est_dfe.out\n";
		print CONFIG3 "neut_egf_file\t".$dir1."/neut_egf.out\n";
		print CONFIG3 "sel_egf_file\t".$dir2."/sel_egf.out\n";
		print CONFIG3 "do_jukes_cantor\t1\n";
		print CONFIG3 "remove_poly\t0\n";

		## Reset at the beginning of each iteration
		$rep_div_scalar = $syn_div_scalar = $rep_sites_scalar = $syn_sites_scalar = $roundrepsites = $roundsynsites = $rep0thbin = $syn0thbin = $poly_rep_total_scalar = $poly_rep_total_scalar = 0;
		%rep_sfs = %syn_sfs = %syngenes = %repgenes = ();
		@rep_sfs = @syn_sfs = ();

	##############################################################
	## Now we need to do all this again for the second gene set:	##
	##############################################################
	if(defined $labelshuffle) {
		foreach my $k (0..$#$resamp2_ref) {
			$gene = $$resamp2_ref[$k];
			print GENES2 "$gene\n";
		
			## To create the 0th bin of the array, we take the total number of sites and subtract those contributing to both divergence and polymorphism (the unmutated sites):
			$rep0thbin = sprintf "%.0f", ($repsites{$gene}-($rep_div{$gene}+$poly_rep_total{$gene}));
			$syn0thbin = sprintf "%.0f", ($synsites{$gene}-($syn_div{$gene}+$poly_syn_total{$gene}));	
#			print "LIST 2: $gene\t$repsites{$gene}\t$roundrepsites\t$rep0thbin\t$rep_div{$gene}\t$poly_rep_total{$gene}\n";
			## Sum all of the polymorphic sites for a given gene
	#		$rep0thbin = sprintf "%.0f", $repsites{$gene}-$poly_rep_total{$gene};
	#		$syn0thbin = sprintf "%.0f", $synsites{$gene}-$poly_syn_total{$gene};
			
			## Dual hashes because duplicates get overwritten:
			$g++;
			$rep_sfs{$g} = $rep0thbin."\t$rep_poly{$gene}";
			$repgenes{$g} = $gene;
			$syn_sfs{$g} = $syn0thbin."\t$syn_poly{$gene}";
			$syngenes{$g} = $gene;

			## Sum the sites and divergences:
			$roundrepsites = sprintf "%.0f", $repsites{$gene};
			$roundsynsites = sprintf "%.0f", $synsites{$gene};
			$rep_sites_scalar += $roundrepsites;		
			$syn_sites_scalar += $roundsynsites;
			$rep_div_scalar += $rep_div{$gene};						## Dn
			$syn_div_scalar += $syn_div{$gene};						## Ds
			$poly_rep_total_scalar += $poly_rep_total{$gene};		## Pn
			$poly_syn_total_scalar += $poly_syn_total{$gene};		## Ps
		}
		foreach $gene (keys %rep_sfs) {
			my @sumrep = split(/\t/, $rep_sfs{$gene});
			my @sumsyn = split(/\t/, $syn_sfs{$gene});

			foreach my $i (0..$#sumrep) {
				$rep_sfs[$i] += $sumrep[$i];
				$syn_sfs[$i] += $sumsyn[$i];
			}
		}
		
		print HANDA2 "Dn\tDs\tDn/Ds\tPn\tPs\tPn/Ps\talpha\tN\n";
		if ($syn_div_scalar > 0) {
			$dnds = $rep_div_scalar/$syn_div_scalar;
		}
		else {
			$dnds = "NA";
		}

		if ($poly_rep_total_scalar > 0) {
			if ($poly_syn_total_scalar > 0) {
				$pnps = $poly_rep_total_scalar/$poly_syn_total_scalar;
			}
		}
		else {
			$pnps = "NA";
		}

		if ($dnds ne "NA" && $pnps ne "NA"){
			if ($dnds > 0) {
				$alpha = 1-(($pnps)/($dnds));
			}
			else {
				$alpha = "NA";
			}
		}
		else {
			$alpha = "NA";
		}

		## Open OUTFILES		
		open(SFS2, "> $sfs2") or die "Cannot write outfile: $!";
		open(DIV2, "> $divout2") or die "Cannot write outfile: $!";
		open(CONFIG1_2, "> $config1_2") or die "Cannot write outfile: $!";
		open(CONFIG2_2, "> $config2_2") or die "Cannot write outfile: $!";
		open(CONFIG3_2, "> $config3_2") or die "Cannot write outfile: $!";
		open(HANDA2, "> $handalphafile2") or die "Cannot write outfile: $!";	
		open(GENES2, "> $genesfile2") or die "Cannot write outfile: $!";

		## Series of print statements
		print HANDA2 "$rep_div_scalar\t$syn_div_scalar\t$dnds\t$poly_rep_total_scalar\t$poly_syn_total_scalar\t$pnps\t$alpha\t$n2\n";
	
		print SFS2 "1\n";
		print SFS2 $sampsize."\n";
		print SFS2 join "\t", @rep_sfs;
		## Add an n+1-th bin of 0s -- this is xi/xi... but with Polymorphorama we do not have that information. So, we use 0 instead.
	# 	print SFS2 "\t0\n";
	## Newest update: rather than 0 in the last bin, we're doing divergence instead (with the rationale being that is every site differs [xi/xi bin] then isn't this divergence? the only other option would be it's the total differences seen relative to the REFERENCE genome... and i still can't suss out what they mean exactly!)
		print SFS2 "\t$rep_div_scalar\n";
	# 	print SFS2 "\n";	## Just new line, no extra 0!	

		print SFS2 join "\t", @syn_sfs;
		## Add an n+1-th bin of 0s -- this is xi/xi... but with Polymorphorama we do not have that information. So, we use 0 instead.
	# 	print SFS2 "\t0\n";	
	#	print SFS2 "\n";	## Just new line, no extra 0!	
		print SFS2 "\t$syn_div_scalar\n";

		print DIV2 "1\t".int($rep_sites_scalar)."\t".$rep_div_scalar."\n";
		print DIV2 "0\t".int($syn_sites_scalar)."\t".$syn_div_scalar."\n";

		print CONFIG1_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG1_2 "sfs_input_file\t$sfs2\n";
		print CONFIG1_2 "est_dfe_results_dir\t".$dir1_2."\n";
		print CONFIG1_2 "site_class\t0\n";
		print CONFIG1_2 "fold\t1\n";
		print CONFIG1_2 "epochs\t2\n";
		print CONFIG1_2 "search_n2\t1\n";
		print CONFIG1_2 "t2_variable\t1\n";
		print CONFIG1_2 "t2\t50\n";
	
		print CONFIG2_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG2_2 "sfs_input_file\t$sfs2\n";
		print CONFIG2_2 "est_dfe_results_dir\t".$dir2_2."\n";
		print CONFIG2_2 "est_dfe_demography_results_file\t".$dir1_2."/est_dfe.out\n";
		print CONFIG2_2 "site_class\t1\n";
		print CONFIG2_2 "fold\t1\n";
		print CONFIG2_2 "epochs\t2\n";
		print CONFIG2_2 "mean_s_variable\t1\n";
		print CONFIG2_2 "mean_s\t-0.1\n";
		print CONFIG2_2 "beta_variable\t1\n";
		print CONFIG2_2 "beta\t0.5\n";
	
		print CONFIG3_2 "data_path_1\t/home/shared/dfe/data/\n";
		print CONFIG3_2 "divergence_file\t".$divout2."\n";
		print CONFIG3_2 "est_alpha_omega_results_file\t".$alphafile2."\n";
		print CONFIG3_2 "est_dfe_results_file\t".$dir2_2."/est_dfe.out\n";
		print CONFIG3_2 "neut_egf_file\t".$dir1_2."/neut_egf.out\n";
		print CONFIG3_2 "sel_egf_file\t".$dir2_2."/sel_egf.out\n";
		print CONFIG3_2 "do_jukes_cantor\t1\n";
		print CONFIG3_2 "remove_poly\t0\n";		

		$i++;
}}

close(SFS); close(DIV); close(CONFIG1); close(CONFIG2); close(CONFIG3); close(SFS2); close(DIV2); close(CONFIG1_2); close(CONFIG2_2); close(CONFIG3_2); close(GENES1); close(GENES2);


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			RUN DFE-alpha		    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
## Currently not a function I am providing
# if(defined $run) {
# 	system("./run-dfe-alpha.pl $set1 $num_iter");
# 	if(defined $glist2) {
# 		system("./run-dfe-alpha.pl $set1 $num_iter");
# 	}
# }

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##    
##			File Organization	    	##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##
my $dir1 = "./dfe_alpha_files_".$set1;
my $dir2 = "./dfe_alpha_files_".$set2;
mkdir $dir1 or die "\nCould not open $dir1: $!\n\n";
system("mv $set1* $dir1");
	if(defined $glist2){
		mkdir $dir2 or die "\nCould not open $dir2: $!\n\n";
		system("mv $set1* $dir2");
	}

## Semi-Global variable assignment outside of subs
my @list;
my @shuffled_indexes;
my @resamp;
my $size;
my $size1;
my $size2;

############## SAMPLES WITHOUT REPLACEMENT FROM A SINGLE GENE SET
sub resample_without {
    ($size, @list) = @_;
	@shuffled_indexes = shuffle(0..$#list);  ## We shuffled the indexes rather than the whole array to speed this up.

	@resamp = @shuffled_indexes[ 0 .. $size-1];		

	my @resampledgenes = @list[ @resamp ];

	return (\@resampledgenes);
}

############## RESAMPLING SUBROUTINE; SAMPLES WITH REPLACEMENT
sub resample_with {
    my(@list) = @_;
	my @resampledgenes = ();
	my $totalgenes = $#list;		## The total number of genes in the array
	for my $i (0..$totalgenes) {
		my $random = rand($totalgenes+1);			## Choose a number between 0 and the total number of genes -- but have to add one because of the array size being shifted to the left
		push(@resampledgenes, $list[$random]);
#		$resampled{$i} = $list[$random];		## Store it into a hash with number as key and the gene as value
	}
	return (\@resampledgenes);
}

############## SAMPLES WITHOUT REPLACEMENT AND SHUFFLES LABELS BETWEEN TWO GENE SETS
sub label_shuffle {
    ($size1, $size2, @list) = @_;
	@shuffled_indexes = shuffle(0..$#list);  ## We shuffled the indexes rather than the whole array to speed this up.
	my %resampled1 = ();
	my %resampled2 = ();
	my $totalgenes = scalar(@list);		## The total number of genes in the array

	## Pick the first N shuffled genes
	my @resamp1 = @shuffled_indexes[ 0 .. $size1-1];		
	my @resamp2 = @shuffled_indexes[ $size1 .. $#list];		## Do NOT subtract 1! You will shortchange yourself by 1 gene if you do

	my @resampledgenes1 = @list[ @resamp1 ];
	my @resampledgenes2 = @list[ @resamp2 ];

	## Lastly, we need to make an array of arrays to return:
	@resamp = (@resampledgenes1, @resampledgenes2);

	return (\@resampledgenes1, \@resampledgenes2);
}

############## CHECKS THE TWO GENE SETS FOR OVERLAPPING ELEMENTS
sub overlap {
    my $l1ref = shift @_;
    my @l1 = @{$l1ref};
    my $l2ref = shift @_;
    my @l2 = @{$l2ref};
    
    my @comb = (@l1, @l2);
    my %seen;
    my %l1 = my %l2 = ();
    my $i = my $j = my $flag = 0;
    my $v = "";
    
    grep { !$seen{$_}++ } @comb;
    
    if (scalar(@comb) == scalar(keys %seen)) {
       print "I found no overlap between the two gene sets.\n\n";
    }

    else {
	## So, if they overlap, we need to return a second list that has the elements from the first list removed.
		print "\nYou gave me overlapping lists. I am deleting the overlap from the second list.\n\n";
	
			for ($i = 0; $i <= $#l1; $i++) {	 
				$l1{$l1[$i]}++;
#   				print "$l1[$i]\n";
			}
			for ($j = 0; $j <= $#l2; $j++) {	 
				$l2{$l2[$j]}++;
			}
			foreach my $key (keys %l1) {
				if (exists $l2{$key}) {
					delete $l2{$key};
				}	
			}
 	
 		##Turn the keys back into an array:
 		my @newl1 = my @newl2 = ();		## Empty them
 		foreach $i (keys %l1) { 
 			push(@newl1, $i);
 		}

 		foreach $i (keys %l2) { 
 			push(@newl2, $i);
 		}
 		@l1 = @newl1;
 		@l2 = @newl2;
 	}
    	return (\@l1, \@l2);			## If the lists had no overlap, they were returned exactly as given. If they did, list2 now has had the list1 elements deleted.
}
  
