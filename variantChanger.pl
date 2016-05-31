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



#-------------------------------------------------------------------------------
# CALLS
my $variantsREF = getSNPinfo($SNPFILE);
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
    return(\%variants);
    # print Dumper \%variants; exit

}

sub changeSNPs {
    my ($seqFile, $outFile, $variantsREF) = @_;
    my $OUTFH = getFH(">", $outFile);
    my $SEQFH = getFH("<", $seqFile);
    my $startPos = 1;
    my $endPos = 0;
    while (<$SEQFH>) {
        if ($. == 1) {
            print $OUTFH $_; next;
        }
        my $lineLen = length($_); #get FASTA sequence char/line count
        $startPos = $startPos + $lineLen unless $. == 2; #seq line start position, ommit 1st line redundancy
        $endPos = $endPos + $lineLen; #seq line end position
        if ( $refPosition >= $startPos and $refPosition <= $endPos) {
            my @seq = split('', $_); #split seq line into array
            my $offset = $refPosition - $startPos; #handle position number within current substring line
            # Handle SNP Type -> modify
            if ($type eq "insertion") {
                splice(@seq, $offset, $variantLen, $calledBase); #insert in sequence
            } elsif ($type eq "deletion") {
                splice(@seq, $offset, $variantLen); #delete in sequence
            } elsif ($type eq "SNP") {
                $seq[$offset] =  $calledBase; #replace SNP in sequence
            } else {
                die "Something is really wrong. Could not modify variant.\n", $!;
            }
            print $OUTFH @seq; #print modified sequence line
            next; #skip to next sequence line
        }
        print $OUTFH $_; #print line to new seq file
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
            say "Variant type ($refBase) at $refPosition taken as an insertion in reference";
        } elsif ($calledBase eq "-") { #probably a deletion in reference
            $type = "deletion";
            $variantLen = length($refBase); #take reference string length
            say "Variant type ($refBase) at $refPosition taken as a deletion in reference";
        } else {
            warn "Could not determine variant type ($refBase), continuing as SNP...";
        }
    } else {
        say "SNP ($refBase -> $calledBase) at $refPosition";
    }
    return ($refPosition, $type, $calledBase, $variantLen);
}
