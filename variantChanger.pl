#!/usr/bin/env perl

use strict; use warnings; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";

use MyIO;

# ==============================================================================
#
#   CAPITAN: Andres Breton
#   FILE: variantChanger.pl
#   LICENSE:
#   USAGE: Change SNPs filtered through snpExtract.pl in sequence
#
# ==============================================================================


#-------------------------------------------------------------------------------
# COMMAND LINE
my $SEQFILE;
my $SNPFILE;
my $OUTFILE;
my $usage = "\n\n$0 [options]\n
Options:
    -seq        Sequence file to change filtered SNPs
    -snp        Filtered SNP file from snpExtract
    -out        Output file name
    -help       Shows this message

";

# OPTIONS
GetOptions(
    'seq=s'             => \$SEQFILE,
    'snp=s'             => \$SNPFILE,
    'out=s'             => \$OUTFILE,
    help                => sub{pod2usage($usage);}
)or pod2usage(2);

checkCLA(); #check command line arguments passed
#-------------------------------------------------------------------------------
# VARIABLES
my $variantsREF = getSNPinfo($SNPFILE);
#-------------------------------------------------------------------------------
# CALLS
changeSNPs($SEQFILE, $OUTFILE, $variantsREF);
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
    unless ($SEQFILE && $SNPFILE && $OUTFILE) {
        die "You did not provide all file names.", $!, $usage;
    }
    return;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($SNPFILE);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes one argument, the SNP file name given as
# a command-line argument. Returns reference of Hash of AoA
# with relevant variant information:
# $refPosition, $type, $calledBase, $variantLen
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = (\%variants);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub getSNPinfo {
    my ($snpFile) = @_;
    my $SNPFH = getFH("<", $snpFile);
    my %variants; #hash of array of arrays to hold all variant info
    while (<$SNPFH>) {
        next if $. == 1; #skip header
        my ($refPosition, $type, $calledBase, $variantLen) = _variantLogic($_); #handle variant logic, get type & length
        # Push anonymous array to %variants hash (Hash of Array of Arrays)
        # Deals with multiple variants at same reference postions, such as
        # a string of insertions
        push @{ $variants{$refPosition} }, [ $type, $calledBase, $variantLen ];
    } close $SNPFH;
    return(\%variants);
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($SEQFILE, $OUTFILE, $variantsREF);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes 3 arguments, the sequence file name to be
# modified, the outout file name, and the reference to Hash of
# AoA from getSNPinfo subroutine.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $output = $outfile with modified variants in sequence
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub changeSNPs {
    my ($seqFile, $outFile, $variantsREF) = @_;
    my %variants =  %$variantsREF; #dereference variants Hash of Array of Arrays
    my $OUTFH = getFH(">", $outFile);
    my $SEQFH = getFH("<", $seqFile);
    my @sortedRefPos = sort { $a<=>$b } keys %variants; #sort and trasnform keys (reference positions) to array
    my $numVars = @sortedRefPos; #variant count
    my $startPos = 1; #sequence start at position 1
    my $endPos = 0;

    while (my $sequence = <$SEQFH>) {
        if ($. == 1) { #deal with FASTA header
            print $OUTFH $sequence; next;
        }
        chomp($sequence);
        my $seqLen = length($sequence); #get FASTA sequence length per line
        $startPos = $startPos + $seqLen unless $. == 2; #seq line start position, ommit 1st sequence line redundancy
        $endPos = $endPos + $seqLen; #seq line end position
        my @seq = split('', $sequence); #split nucleotides from seq line into array
        my $IDXindel = 0; #manage Indel index fluxuations after inserting/deleting in sequence

        for (my $i = 0; $i < $numVars; $i++) { #itereate through variants in order of reference position
            my $refPosition = $sortedRefPos[$i];
            if ( $refPosition >= $startPos and $refPosition <= $endPos) {
                foreach my $varOccurance ( @{$variants{$refPosition}} ) { #get each occurance at reference position -> handle multiple for same position
                    # Array from Hash of AoA with $type, $calledBase, $variantLen
                    # Sets each accordingly
                    my @variantInfo = @$varOccurance;
                    my $type = $variantInfo[0];
                    my $calledBase = $variantInfo[1];
                    my $variantLen = $variantInfo[2];

                    # Handle SNP Type -> modify
                    my $offset = $refPosition - $startPos + $IDXindel; #handle position number within current substring sequence line
                    if ($type eq "insertion") {
                        $offset++; #insert 1 postion over $refPosition
                        splice(@seq, $offset, 0, $calledBase); #insert in sequence
                        $offset--; #reset offset to original
                        $IDXindel = $IDXindel + $variantLen;
                    } elsif ($type eq "deletion") {
                        splice(@seq, $offset, $variantLen); #delete in sequence
                        $IDXindel = $IDXindel - $variantLen;
                    } elsif ($type eq "SNP") {
                        $seq[$offset] =  $calledBase; #replace SNP in sequence
                    } else {
                        die "Something is really wrong. Could not modify variant.\n", $!;
                    }
                }
                $sequence = join('', @seq);
            }
        }
        say $OUTFH $sequence; #print line not needing changes to new seq file
    } close $SEQFH; close $OUTFH;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $input = ($_);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This function takes one argument, the current line from the SNP
# file. Processes file and returns relevant variant information:
# $refPosition, $type, $calledBase, $variantLen
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# $return = ($refPosition, $type, $calledBase, $variantLen);
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub _variantLogic {
    my ($line) = @_;
    my @items = split("\t", $line);
    my ($refPosition, $type, $refBase, $calledBase) = ($items[3], $items[4], $items[5], $items[6]); #get variant type Indel/SNP

    my $variantLen = 0;
    # Handle insertion vs. deletion logic, otherwise == SNP
    if ($type eq "Indel") {
        if ($refBase eq "-") { #probably an insertion in reference
            $type = "insertion";
            $variantLen = length($calledBase); #take called base string length being inserted
            say "INSERTION\t($calledBase) at position $refPosition in reference";
        } elsif ($calledBase eq "-") { #probably a deletion in reference
            $type = "deletion";
            $variantLen = length($refBase); #take reference string length being deleted
            say "DELETION\t($refBase) at position $refPosition in reference";
        } else {
            warn "Could not determine variant type ($refBase), continuing as SNP...";
        }
    } else {
        say "SNP\t\t($refBase -> $calledBase) at position $refPosition in reference";
    }
    return ($refPosition, $type, $calledBase, $variantLen);
}
