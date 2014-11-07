#!/usr/local/bin/perl
$|++;
use warnings;
use strict;
use File::Basename;

my $assemblyFA = shift;
my $projectname = shift;
my $geneFA = shift;

if(!defined($geneFA)) {
        die "USAGE: perl $0 need assembly fasta, project and gene fasta file names\n";}

`formatdb -i $assemblyFA -p F`;
`blastall -d $assemblyFA -i $geneFA -o BlastOutput.blast -p blastn -e 1e15`;
`perl useParse.pl BlastOutput.blast`;

my $bp = "BlastOutput.bp";
my $outputFA = $assemblyFA . "_" . $projectname . "_extracted.fa";

open BP, $bp or die "could not open bp file\n";
open AF, $assemblyFA or die "could not open assembly fasta file\n";
open OFA, ">$outputFA" or die "could not open output file\n";


my ($line, $sID, $sStart, $sEnd, $seqlength);
my $FAbasename = basename($assemblyFA);
my $subline;

while ($line = <BP>) {
        chomp $line;

        if ($line =~ /SUB:(\w+.*?)\s*/) {
                $sID = $1;
                }


        if ($line =~ /.*sb:(\d+)\sse:(\d+)\s.*/) {
                ($sStart, $sEnd) = ($1, $2);
                while ($line = <AF>) {
                        chomp $line;
                        if ($line =~ /^>$sID\s.*/) {
                                $line = <AF>;
                                chomp $line;
                                if ($sEnd > $sStart) {
                                        $subline = substr($line, $sStart - 1, $sEnd - $sStart + 1);
                                        } else {$subline = substr($line, $sEnd - 1, $sStart - $sEnd + 1);
                                        $subline = reverse($subline);
                                        }
                                $seqlength = length($subline);
                                print OFA ">$FAbasename|$sID|$seqlength|$sStart|$sEnd\n$subline\n";
                                }
                        }
                }
                close AF;
                open AF, $assemblyFA or die "could not reopen assembly fasta file\n";

        }


close AF;
close BP;
close OFA;


