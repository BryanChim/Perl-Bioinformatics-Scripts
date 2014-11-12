#!/usr/local/bin/perl

#### fastaSTATS.pl

#    PURPOSE: Calculate basic contig statistics for a given FASTA file, including n50
#
#    INPUT: From $inputfasta, via command-line input 
#		(NOTE: to be passed in via fastaMAIN.pl)
#
#    OUTPUT: From $stats_outputfile, generated as "$inputfasta_stats.txt"
#			Will contain the following metrics, calculated from input FASTA file: 
#				- Number of Contigs Found 
#				- Summed Length of All Contigs
#				- GC Total Count
#				- GC Total Percent
#				- Average Contig Length
#				- Smallest Contig Length
#				- Largest Contig Length
#				- Total Number of n50 Contigs*
#				- Length of Smallest n50 Contigs, aka the "n50"
#				- Average Length of n50 Contigs*
#
#	  *NOTE: n50 Contigs refer to the minimal set of largest contigs that, when summed together,
#		exceed more than 50% length of ALL contigs found in the FASTA file
#		-- serves as a measure of length-weighted distribution

# libraries and pragmas
use warnings;
use strict;
$|++;

# command-line variable -- $inputfasta -- to be passed in via fastaMAIN.pl 
my $inputfasta = shift;

# check for existence of $inputfasta, else die
if (!defined ($inputfasta)) 
	{die "USAGE: perl $0 need input fasta file";}

# scalar variables, used to hold statistical information while parsing
my ($line, $seqID, $gctotal, $gcpct, $avgcontig,
$seqlen, $n50lengthSum, $n50contigCount, $n50MIN, $n50AVG);
$n50lengthSum = 0;

# array variables - @seqlengths holds lengths of all sequences found
my @seqlengths;

# declaration and initialization of %seqHash
# holds count of sequences, sum of sequence lengths and tracks min/max lengths
my %seqHash;
$seqHash{count} = 0;
$seqHash{length} = 0;
$seqHash{min} = 100000000;
$seqHash{max} = 0;

# output file
my $stats_outputfile = $inputfasta . "_stats.txt";

# open file-handles for input FASTA and output stats file
open AF, $inputfasta or die "could not open input FASTA file\n";
open STO, ">$stats_outputfile" or die;

# start parsing input FASTA file
while ($line = <AF>) 
{

	chomp $line;

	# get sequence ID, go to next line to get sequence
	if ($line =~ /^>(.*?)\s\[.*/) 
	{
		$seqID = $1;
		$line = <AF>;
		chomp $line;
		$seqHash{$seqID} = $line;
		
		# increment sequence count, add sequence length to summed total length
		$seqHash{count}++;
		$seqHash{length} += length($line);
		
		# check if length of sequence exceeds current max - if so update the max
		if (length($line) > $seqHash{max}) 
		{
			$seqHash{max} = length($line);
		}
		
		# check if length of sequence is lower than current min - if so update the min
		if (length($line) < $seqHash{min}) 
		{
			$seqHash{min} = length($line);
		}

		# iterate through each nucleotide, increment GC total as necessary
		for (my $i = 0; $i < length($line); $i++) 
		{
			my $seqchar = substr($line, $i, 1);
					
			if ($seqchar eq 'g' || $seqchar eq 'G'
			|| $seqchar eq 'c' || $seqchar eq 'C') 
			{
				$gctotal++;
			}
		}			
	}
} # finish parsing input FASTA file


# calculate GC percent
$gcpct = ($gctotal/$seqHash{length}) * 100;

# calculate average contig length
$avgcontig = $seqHash{length}/$seqHash{count};

# extract all the sequence lengths, push into @seqlengths array
foreach my $seq (keys %seqHash) 
{
	$seqlen = length($seqHash{$seq});
	push @seqlengths, $seqlen;
}

# For loop to calculate n50 statistics:
# - Sorts and iterates through all contig lengths from longest to shortest
# - Keep running sum of contig lengths UNTIL the sum exceeds %50 of length of ALL contigs
# - $n50MIN ultimately tracks the last and smallest contig to be summed AKA n50
# - Increment $n50contigCount for each contig that is summed up in $n50lengthSum
foreach my $len (sort {$b <=> $a} @seqlengths) 
{
	if ($n50lengthSum < ($seqHash{length}/2)) 
	{
		$n50lengthSum += $len;
		$n50MIN = $len;
		$n50contigCount++;
	}
}

# calculate average length of those contigs that contribute to $n50lengthSum
$n50AVG = $n50lengthSum/$n50contigCount;

# output to the stats output file ($inputfasta_stats.txt)
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
        
# end of program -- close files        
close AF;
close STO;
