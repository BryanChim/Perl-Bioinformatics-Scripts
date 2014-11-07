#!/usr/local/bin/perl
$|++;
use warnings;
use strict;

# accept inputs for filenames, from command line
my $inputFA = shift;
my $seqID = shift;
my $outputFA = shift;

if (!defined ($outputFA)) {
        die "USAGE: perl $0 need inputFA seqID and outputFA\n";
        }

open INFA, $inputFA or die "could not open input FA file";
open SEQ, $seqID or die "could not open sequence id file";
open OUTFA, ">$outputFA" or die "could not open output FA file";


my $line; # stores lines in SEQ and INFA for processing
my $subline;  # stores the truncated sequences
my $id; # stores accession IDs
my $start = 0; # stores desired sequence start values
my $end = 0; # stores desired sequence end values
my @linearray; # store lines of SEQ

$line = <SEQ>;

# read through and store lines of SEQ in @linearray
while ($line = <SEQ>) {
        chomp $line;
        push @linearray, $line;
        }

close SEQ;

        # iterate through each SEQ line, get the ID, start, end values
        foreach my $lin (@linearray) {
                ($id, $start, $end) = $lin =~ /(.*?)\s+(\d+)\s+(\d+)/;
                
                        # read and chomp lines from INFA, regex to search for
                        # accession ID - if found, go to next line for sequence
                        while ($line = <INFA>) {
                        chomp $line;
                        if ($line =~ /^>$id.*/) {
                                $line = <INFA>;

                                # use the start and end values to extract a substring from the sequence
                                $subline = substr($line, $start - 1, $end - $start + 1);
                                
                                # print to OUTFA: accession ID, start, end, and the extracted sequence
                                print OUTFA ">$id | start=$start | end=$end\n$subline\n";
                                }
                        }
                        close INFA;
                        open INFA, $inputFA or die "could not open input FA file";
                }

 close INFA;
 close OUTFA;
                


        




