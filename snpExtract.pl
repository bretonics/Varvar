#!/usr/bin/env perl

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use MyIO;

# ==============================================================================
#
#   CAPITAN: Andres Breton
#   FILE: snpExtract.pl
#   LICENSE:
#   USAGE: Extract desired SNPs using columns' values as bounds to call true SNPS
#
# ==============================================================================

#-------------------------------------------------------------------------------
# VARIABLES
my $outFile = "out.txt";
#-------------------------------------------------------------------------------
# COMMAND LINE
my $FILE;
my @SNP;
my @DEPTH;
my $usage = "

 [options]

Options:
    -snp        SNP percentage (float) lower and upper bound; [90.0 100.0]
    -depth      Depth of SNP coverage bound; [0 Inf]
    -help       Shows this message

";


# OPTIONS
GetOptions(
    'file=s'                => \$FILE,
    'snp=f{0,100}'          => \@SNP,
    'depth=i{0,}'           => \@DEPTH,
    help                    => sub{pod2usage($usage);}
)or pod2usage(2);

checkCLA(\@SNP, \@DEPTH); #check command line arguments passed
#-------------------------------------------------------------------------------
# CALLS
extractColumns($FILE, $outFile, \@SNP, \@DEPTH);

#-------------------------------------------------------------------------------
# SUBS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = argChecks();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function checks for arguments passed in command-line
# using global variables
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = Dies from errors
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub checkCLA {
    # unless ($FILE) {
    #     die "You did not provide an input file.", $!, $usage;
    # }
    unless (@SNP && @DEPTH) {
        die "Did not provide snp and/or depth bounds to use as limits.", $!, $usage;
    }
    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes arguments
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub extractColumns {
    my ($file, $outFile, $snp, $depth) = @_;
    my @snp = @$snp;
    my $snpMin = $snp[0];
    my $snpMax = 100.0;
    $snpMax = $snp[1] if ($snp[1]);

    my @depth = @$depth;
    my $depthMin = $depth[0];
    my $depthMax = "inf";
    $depthMax = $depth[1] if ($depth[1]);

    my $i = 1; my $snpCount = 0;
    my $FH = getFH("<", $file);
    while (<$FH>) {
        next if $. == 1; #skip header
        my @line = split /\t/; #split line by tab delim
        my ($snpPer) = $line[10] =~ /(^\d+\.\d+)\s\%/; #get percentage float
        if ($snpPer <= $snpMax && $snpPer >= $snpMin) { #check SNP % is within bounds
            my $depth = $line[13];
            if ($depth <= $depthMax && $depth >= $depthMin) { #check Depth is within bounds
                say "SNP $i is inside bounds- SNP %: $snpPer \t Depth: $depth";
                my $outFH = getFH(">", $outFile);
                say $outFH $_; #write passing SNP to file
                $snpCount++;
            } else { #Depth is outside of bounds
                # say "SNP $i is outside depth bounds- Depth: $depth";
            }
        } else { #SNP % is outside of bounds
            # say "SNP $i is outside SNP % bounds- SNP %: $snpPer";
        }
        $i++;
    }
    say "-----------------------------------------------------\nTotal number of SNPs: $i";
    say "Passing: $snpCount\n";
    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ();
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes arguments
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output =
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
