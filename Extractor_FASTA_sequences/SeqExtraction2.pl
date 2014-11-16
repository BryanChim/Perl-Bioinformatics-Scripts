#!/usr/local/bin/perl

#### SeqExtraction2
#
#	 VERSION 2 MODIFICATIONS: 
#			** Opted to store each set of sequenceID-start-end
#			   in a 2-dimensional array using array references for easier iteration.
#			** Iteration through this array (@array_of_seq_info) is performed
#			   whenever a new header line and sequence is found, then writes to output as normal	  
#			** Much more memory efficient than VERSION 1 -- no longer opens/closes INFA repeatedly
#
#    PURPOSE: Extracts specific fragments of sequences from a FASTA file
#			using a read-in list of Accession IDs, start coordinates and 
#			end coordinates.
#		
#    INPUT: From $seqID, contains IDs of desired sequences to be extracted:
#				Format: <Accession ID>	<Start> <End>
#			From $inputFA, a nucleotide sequence file in FASTA format
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
my $seqID = shift;
my $inputFA = shift;
my $outputFA = shift;

if (!defined ($outputFA)) 
	{die "USAGE: perl $0 need inputFA seqID and outputFA\n";}

open SEQ, $seqID or die "could not open sequence id file";
open INFA, $inputFA or die "could not open input FA file";
open OUTFA, ">$outputFA" or die "could not open output FA file";

#################### VARIABLES ######################
my $line; 	# stores lines in SEQ and INFA for processing
my $header_line;	# stores previous header line
my $sequence_line;	# stores previously read sequence line
my $subline;  	# stores the truncated sequences
my $id = 0; 	# preassigned array index for accession ID
my $start = 1; 	# preassigned array index for start coordinate
my $end = 2; 	# preassigned array index for end coordinate
my $current_id;	# current sequence ID to try to match
my $current_start;	# current start coordinate used when extracting truncated sequence
my $current_end;	# current end coordinate used when extracting truncated sequence
my @array_of_seq_info;	# stores a array to a list of references to seq info
		# ie. the Accesssion IDs, Starts and Ends from $seqID file

################### MAIN PROGRAM ####################

# ignore column headers
$line = <SEQ>;

# read through $seqID file and store each set of (ID, Start, End) into @array_of_seq_info
while ($line = <SEQ>) 
{
	chomp $line;
	my @linearray = split('\t', $line);
	push @array_of_seq_info, \@linearray;
}

close SEQ;

# read through input FASTA file
while ($line = <INFA>) 
{
	chomp $line;
	
	# check for new header-- if found, store the header line and the sequence line
	if ($line =~ /^>/) 
	{
		$header_line = $line;
		$sequence_line = <INFA>;
		
		# iterate through each ID in our seq info array, write output FASTA if ID match found
		foreach my $seq_info_array_ref (@array_of_seq_info)
		{
			$current_id = ${$seq_info_array_ref}[$id];
			
			if ($header_line =~ /^>$current_id/)
			{
				# get start and end coordinates for the currently matched accession ID
				$current_start = ${$seq_info_array_ref}[$start];
				$current_end = ${$seq_info_array_ref}[$end];
				
				# use the start and end values to extract a substring from the sequence
				$subline = substr($sequence_line, $current_start - 1, $current_end - $current_start + 1);
				
				# print to OUTFA: accession ID, start, end, and the truncated sequence
				print OUTFA ">$current_id | start=$current_start | end=$current_end\n$subline\n";			
			
			}
		}
	}
}
	
#### END OF PROGRAM; CLOSE FILES
 close INFA;
 close OUTFA;