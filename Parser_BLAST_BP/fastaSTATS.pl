#!/usr/local/bin/perl
use warnings;
use strict;
$|++;

my $inputfasta = shift;

if (!defined ($inputfasta)) {
        die "USAGE: perl $0 need input fasta file";}

my $stats_outputfile = $inputfasta . "_stats.txt";

open AF, $inputfasta or die "could not open input fasta file\n";
open STO, ">$stats_outputfile" or die;

my ($line, $seqID, $gctotal, $gcpct, $avgcontig,
$seqlen, $n50lengthSum, $n50contigCount, $n50MIN, $n50AVG);
$n50lengthSum = 0;
my @seqlengths;
my %seqHash;
$seqHash{count} = 0;
$seqHash{length} = 0;
$seqHash{min} = 100000000;
$seqHash{max} = 0;


while ($line = <AF>) {

        chomp $line;

        if ($line =~ /^>(.*?)\s\[.*/) {
                $seqID = $1;
                $line = <AF>;
                chomp $line;
                $seqHash{$seqID} = $line;
                $seqHash{count}++;
                $seqHash{length} += length($line);
                if (length($line) > $seqHash{max}) {
                        $seqHash{max} = length($line);
                        }
                if (length($line) < $seqHash{min}) {
                        $seqHash{min} = length($line);
                        }

                for (my $i = 0; $i < length($line); $i++) {
                        my $seqchar = substr($line, $i, 1);
                        if ($seqchar eq 'g' || $seqchar eq 'G'
                                || $seqchar eq 'c' || $seqchar eq 'C') {
                                 $gctotal++;
                                 }
                         }
                
                }
                }
$gcpct = ($gctotal/$seqHash{length}) * 100;

$avgcontig = $seqHash{length}/$seqHash{count};

foreach my $seq (keys %seqHash) {
        $seqlen = length($seqHash{$seq});
        push @seqlengths, $seqlen;
        }
foreach my $len (sort {$b <=> $a} @seqlengths) {
        if ($n50lengthSum < ($seqHash{length}/2)) {
                $n50lengthSum += $len;
                $n50MIN = $len;
                $n50contigCount++;
                }
        }
$n50AVG = $n50lengthSum/$n50contigCount;


print STO "No of Sequences = $seqHash{count}\n";
print STO "Total Sequence Length = $seqHash{length}\n";
print STO "GC Total = $gctotal\n";
print STO "GC Perc = $gcpct;\n";
print STO "Avg Contig Size =  $avgcontig\n";
print STO "Min Contig Size = $seqHash{min}\n";
print STO "Max Contig Size = $seqHash{max}\n";
print STO "Total N50 contigs = $n50contigCount\n";
print STO "N50ContigSize = $n50MIN\n";
print STO "N50Avg = $n50AVG";
                
close AF;
close STO;
