#!/usr/local/bin/perl
use warnings;
use strict;

my $prot_useMarkov_outputfile = "BryanChim_6803_prot-scores.txt";
my $filter_outputfile = "BryanChim_6803_prot-scores_filtered.txt";

open FH, $prot_useMarkov_outputfile or die "could not open wfeaesdf\n";
open FHO, ">$filter_outputfile" or die "could not open output\n";

my $line;
my @linearray;
my $score;
my $LF = "\n";

while ($line = <FH>) {
        chomp $line;
        @linearray = split("\t", $line);
        $score = $linearray[0];
        if ($score < -8000) {
        print FHO $line, $LF;
        }
        }
