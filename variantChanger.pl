#!/usr/bin/env perl

use strict; use warnings; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use FindBin; use lib "$FindBin::RealBin/lib";
use Data::Dumper qw(Dumper);
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
    -snp        SNP file from snpExtract containing filtered SNPs to change
    -out        Out file name
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

sub getSNPinfo {
    my ($snpFile) = @_;
    my $SNPFH = getFH("<", $snpFile);
    my %variants; #hash of array of arrays with variant info
    while (<$SNPFH>) {
        next if $. == 1; #skip header
        my ($refPosition, $type, $calledBase, $variantLen) = _variantLogic($_); #handle variant logic, get type & length
        # Push anonymous array to %variants hash ==> Hash of Array of Arrays
        # Deals with multiple variants at same reference postions, such as
        # a string of insertions
        push @{ $variants{$refPosition} }, [ $type, $calledBase, $variantLen ];
    } close $SNPFH;
        # print Dumper \%variants; exit
    return(\%variants);
}

sub changeSNPs {
    my ($seqFile, $outFile, $variantsREF) = @_;
    my %variants =  %$variantsREF; #dereference variants Hash of Array of Arrays
    my $OUTFH = getFH(">", $outFile);
    my $SEQFH = getFH("<", $seqFile);
    my @sortedRefPos = sort { $a<=>$b } keys %variants; #sort and trasnform to array
    my $numVars = @sortedRefPos;
    my $startPos = 1;
    my $endPos = 0;

    while (my $sequence = <$SEQFH>) {
        if ($. == 1) {
            print $OUTFH $sequence; next;
        }
        chomp($sequence);
        # my ($sequence) = chomp($_);
        my $seqLen = length($sequence); #get FASTA sequence char/line count
        $startPos = $startPos + $seqLen unless $. == 2; #seq line start position, ommit 1st line redundancy
        $endPos = $endPos + $seqLen; #seq line end position
        my @seq = split('', $sequence); #split seq line into array
        my $IDXindel = 0;
# say "\nLine $. length: $seqLen start: $startPos \n$sequence";
        for (my $i = 0; $i < $numVars; $i++) { #itereate through variants in order of reference position
            my $refPosition = $sortedRefPos[$i];
            if ( $refPosition >= $startPos and $refPosition <= $endPos) {
                foreach my $varOccurance ( @{$variants{$refPosition}} ) { #get each occurance at reference position -> handle multiple for same position

                    my @variantInfo = @$varOccurance; #array with $type, $calledBase, $variantLen
                    my $type = $variantInfo[0];
                    my $calledBase = $variantInfo[1];
                    my $variantLen = $variantInfo[2];
                    # my @seq = split('', $seqLine); #split seq line into array
# say "INDX = $IDXindel";
                    my $offset = $refPosition - $startPos + $IDXindel; #handle position number within current substring line
# say "Offset Before: $offset";
                    # Handle SNP Type -> modify
                    if ($type eq "insertion") {
                        $offset++; #insert 1 postion over $refPosition
                        splice(@seq, $offset, 0, $calledBase); #insert in sequence
                        $offset--; #reset offset to original
                        $IDXindel = $IDXindel + $variantLen;
# say "index insert: $IDXindel";
                    } elsif ($type eq "deletion") {
                        splice(@seq, $offset, $variantLen); #delete in sequence
                        $IDXindel = $IDXindel - $variantLen;
# say "index deletion: $IDXindel";
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


sub _variantLogic {
    my ($line) = @_;
    my @items = split("\t", $line);
    my ($refPosition, $type, $refBase, $calledBase) = ($items[3], $items[4], $items[5], $items[6]); #get variant type Indel/SNP

    my $variantLen = 0;
    if ($type eq "Indel") { # Handle insertion vs. deletion logic, otherwise == SNP
        if ($refBase eq "-") { #probably an insertion in reference
            $type = "insertion";
            $variantLen = length($calledBase); #take called base string length
            say "INSERTION\t($calledBase) at position $refPosition in reference";
        } elsif ($calledBase eq "-") { #probably a deletion in reference
            $type = "deletion";
            $variantLen = length($refBase); #take reference string length
            say "DELETION\t($refBase) at position $refPosition in reference";
        } else {
            warn "Could not determine variant type ($refBase), continuing as SNP...";
        }
    } else {
        say "SNP\t\t($refBase -> $calledBase) at position $refPosition in reference";
    }
    return ($refPosition, $type, $calledBase, $variantLen);
}
