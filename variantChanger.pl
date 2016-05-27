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
        die "You did not provide all file names.", $!, $usage;
    }
    return;
}

sub changeSNPs {
    my ($seqFile, $snpFile, $outFile) = @_;
    my $SNPFH = getFH("<", $snpFile);
    my $OUTFH = getFH(">", $outFile);
    while (<$SNPFH>) {
        next if $. == 1; #skip header
        my ($refPosition, $type, $variant, $variantLen) = _variantLogic($_); #handle variant logic, get type & length
say $refPosition, $type, $variant, $variantLen; exit;
        # Open Sequence File
        my $SEQFH = getFH("<", $seqFile);
        my $startPos = 1;
        my $endPos = 0;
        while (<$SEQFH>) {
            print $OUTFH $_; #print line to new seq file
            next if $. == 1;
            my $lineLen = length($_); #get FASTA sequence char/line count
            $startPos = $startPos + $lineLen unless $. == 2; #seq line start position, ommit 1st line redundancy
            $endPos = $endPos + $lineLen; #seq line end position
            if ( $refPosition >= $startPos and $refPosition <= $endPos) {
                my @seq = split('', $_); #split seq line into array
                my $offset = $refPosition - $startPos; #handle position number within current substring line
                # Handle SNP Type -> modify
                if ($type eq "insertion") {
                    splice(@seq, $offset, $variantLen,  $variant,); #insert in sequence
                } elsif ($type eq "deletion") {
                    splice(@seq, $offset, $variantLen); #delete in sequence
                } elsif ($type eq "SNP") {
                    $seq[$offset] =  $variant; #replace snp in sequence
                } else {
                    die "Something is really wrong. Could not modify variant ($variant).\n", $!;
                }
                print $OUTFH @seq; #print modified sequence line
            }
        } close $SEQFH; close $OUTFH;
    } close $SNPFH;
}


sub _variantLogic {
    my ($line) = @_;
    my ($refPosition, $type, $variant) = $line =~ /^\?\t.+\t\d+\t(\d+)\t(\w+)\t(\w+)\t\w+\t.*/; #get variant type Indel/SNP
    my $variantLen = 0;
    if ($type eq "Indel") { # Handle insertion vs. deletion logic, otherwise == SNP
        if ($variant eq "-") { #probably an insertion in reference
            $type = "insertion";
            # Length of deletion may be > 1
                        say "Variant type ($variant) taken as a deletion in reference";
        } elsif ($variant =~ /\w+/) { #probably a deletion in reference
            $type = "deletion";
            $variantLen = length($variant); #take reference string length
            say "Variant type ($variant) taken as an insertion in reference";
        } else {
            warn "Could not determine variant type ($variant), continuing as SNP...";
        }
    }
    return ($refPosition, $type, $variant, $variantLen);
}
