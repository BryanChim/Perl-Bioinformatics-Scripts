#!/usr/local/bin/perl

#### SeqExtraction
#
#	 VERSION 1: Memory-inefficient -- closes reopens filehandle of FASTA
#			for EVERY sequence to be extracted
#
#    PURPOSE: Extracts specific fragments of sequences from a FASTA file
#			using a read-in list of Accession IDs, start coordinates and 
#			end coordinates.
#		
#    INPUT: From $inputFA, a nucleotide sequence file in FASTA format
#			From $seqID, contains IDs of desired sequences to be extracted:
#				Format: <Accession ID>	<Start> <End>
#
#    OUTPUT: From file $outputFA, contains extracted sequences in the format:
#				<Accession ID> | <start> | <end>
#				<sequence>
#

############### LIBRARIES AND PRAGMAS ###############
$|++;
use warnings;
use strict;

####################### FILES #######################
my $inputFA = shift;
my $seqID = shift;
my $outputFA = shift;

if (!defined ($outputFA)) 
	{die "USAGE: perl $0 need inputFA seqID and outputFA\n";}

open INFA, $inputFA or die "could not open input FA file";
open SEQ, $seqID or die "could not open sequence id file";
open OUTFA, ">$outputFA" or die "could not open output FA file";

#################### VARIABLES ######################
my $line; 	# stores lines in SEQ and INFA for processing
my $subline;  	# stores the truncated sequences
my $id; 	# stores accession IDs
my $start = 0; 	# stores desired sequence start values
my $end = 0; 	# stores desired sequence end values
my @linearray; 	# store lines of SEQ


################### MAIN PROGRAM ####################

$line = <SEQ>;

# read through and store lines of SEQ in @linearray
while ($line = <SEQ>) 
{
	chomp $line;
	push @linearray, $line;
}

close SEQ;

# iterate through each SEQ line, get the ID, start, end values
foreach my $lin (@linearray) 
{
	($id, $start, $end) = $lin =~ /(.*?)\s+(\d+)\s+(\d+)/;
		
	# read and chomp lines from INFA, regex to search for
	# accession ID - if found, go to next line for sequence
	while ($line = <INFA>) 
	{
		chomp $line;
		
		if ($line =~ /^>$id/) 
		{
			$line = <INFA>;

			# use the start and end values to extract a substring from the sequence
			$subline = substr($line, $start - 1, $end - $start + 1);
			
			# print to OUTFA: accession ID, start, end, and the extracted sequence
			print OUTFA ">$id | start=$start | end=$end\n$subline\n";
		}
		
	}
	
	# close and reopen input FASTA file for next sequence
	close INFA;
	open INFA, $inputFA or die "could not open input FA file";
}

#### END OF PROGRAM; CLOSE FILES
 close INFA;
 close OUTFA;
                


        




