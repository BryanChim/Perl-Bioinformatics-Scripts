#!/usr/local/bin/perl
use warnings;
use strict;

my $inputfastadir = shift;
my $geneFile = shift;

if (!defined ($geneFile)) {
        die "USAGE perl $0 need inputs for fasta directory and gene file\n";
        }

my @files = <$inputfastadir/*.fsa>;
    foreach my $filename (@files) {
	if(-d $filename) {
	    listFiles($filename);
	} else {
	    processFile($filename);
	}
	}

sub processFile {
my ($inputfa) = @_;
`perl fastaStats.pl $inputfa`;
`perl fastaBLASTBP.pl $inputfa ecoli_16s $geneFile`;
}
