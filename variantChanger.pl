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
my $usage = "\n\n $0 [options]\n
Options:
    -seq    Sequence file to change filtered SNPs
    -snp    SNP file from snpExtract containing filtered SNPs to change
    -out    Out file name
    -help   Shows this message

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
changeSNPs($SEQFILE, $SNPFILE, $OUTFILE);

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
        die "You did not provide all file names.", $!;
    }
    return;
}

sub changeSNPs {
    my ($seqFile, $snpFile, $outFile) = @_;
    my $lineLen = _getLineLen($seqFile); #get FASTA sequence char/line count
    my $snpFH = getFH("<", $snpFile);
    my $outFH = getFH(">", $outFile);
    while (<$snpFH>) {
        next if $. == 1; #skip header
        my $position = $_ =~ /^?\t.+\t.+\t(\d+)\t.*/;
        my $type = _variantLogic($_); #handle variant logic

        # Open Sequence File
        my $seqFH = getFH("<", $seqFile);
        my $startPos = 1;
        my $endPos = $lineLen;
        while (<$seqFH>) {
            next if $. == 1;
            $startPos = $startPos + $lineLen unless $. == 2; #seq line start position count
            $endPos = $endPos + $lineLen; #seq line end position count
            if ( $startPos < $position > $endPos) {
                my @seq = split('', $_); #split seq line into array
                # Handle SNP Type
                if ($type eq "insert") {
                    splice(@seq, $position, $varLen, $snp); #insert in sequence
                } elsif ($type eq "delete") {
                    splice(@seq, $position, $varLen, $snp); #delete in sequence
                } elsif ($type eq "snp") {
                    splice(@seq, $position, $varLen, $snp); #replace snp in sequence
                } else {
                    #something is really wrong
                }
                return @seq; #return modified sequence line
            }

        }

    }


    _getLineNum();
}

sub _getLineLen {
    my ($file) = @_
    while (<getFH("<", $file)>) {
        next if $. == 1;
        my $lineLen = length($_);
        return $lineLen;
    }
}

sub _variantLogic {
    my ($line) = @_;
    my $type = $line =~ /^?\t.+\t.+\t\d+\t(.+)\t.*/; #get variant type INDEL/SNP
    my $varLen = 1;
    if ($type eq "Indel") {
        # Handle insertion vs. deletion logic
        $varLen = $line =~ /^?\t.+\t.+\t\d+\t.+\t.+\t(.+)\t.*/;
        $varLen = length($varLen); #length of insertion
    }
    return ($type, $varLen);
}

sub _getLineNum {
    my () = @_;
    return;
}
