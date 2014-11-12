#!/usr/local/bin/perl

#### fastaMAIN.pl

#    PURPOSE: Main program for a BLASTn pipeline which calls programs, fastaSTATS.pl and fastaBLASTBP.pl
#		- Calculates and outputs contig statistics, including n50
#		- Runs BLASTn, using input gene ($geneFile) against all FASTA files in input directory
#		- Converts BLASTn output to condensed .bp format
#		- Parses .bp for sequence IDs and start/end coordinates 
#		- Extract raw sequence of hits from original FASTA, outputs data in a FASTA format 
#
#    INPUT: From $inputfastadir (FIRST argument - directory containing FASTA files to be analyzed)
#				- each FASTA file should have the .fsa extension
#			From $geneFile (SECOND argument - FASTA file containing gene of interest)
#
#    OUTPUT: 
#		- .fsa_stats for each FASTA (produced by fastaSTATS.pl)
#		- BLAST database files for each input FASTA (.nhr, .nin, .nsq)
#		- BLAST outputs and .bp files (BlastOutput.blast and BlastOutput.bp)
#		- Final FASTA output for each run ("<FASTA_basename>_<project_basename>_extracted.fa")
#		** See fastaSTATS.pl and fastaBLASTBP.pl for more details
#
#	 NOTE: Command-line BLAST needs to be installed on your computer for fastaBLASTBP.pl to run
#

# libraries and pragmas
use warnings;
use strict;

# command-line inputs
my $inputfastadir = shift;
my $geneFile = shift;
my $appropriate_FASTA_input = 0;

# check for definition of $geneFile or die
if (!defined ($geneFile)) {
        die "USAGE perl $0 need inputs for fasta directory and gene file\n";
        }


# if 1st input is not a directory, but is instead a single FASTA file, process it
if (! -d $inputfastadir && $inputfastadir =~ /.*\.fsa/) 
{
	processFile($inputfastadir);
	$appropriate_FASTA_input = 1;
	print "FASTA file BLASTED and processed!";
}


# otherwise, compile a list of all the FASTA files within the input directory 
else
{
	my @files = <$inputfastadir/*.fsa>;

	# for each FASTA file found, run processFile subroutine
	foreach my $filename (@files) 
	{

		#print "File: $filename\n";
		 processFile($filename);
		 $appropriate_FASTA_input = 1;
		
	}
	
	# notification of success
	if ($appropriate_FASTA_input) 
	{
		print "All FASTA files BLASTED and processed!";
	}
	
}

# notification of inappropriate input
if (!$appropriate_FASTA_input)
{
	die "USAGE perl $0 input FASTA or FASTA directory is invalid\n";
}	

#### processFile #### runs fastaSTATS.pl and fastaBLASTBP.pl	 
# 	fastaSTATS.pl: calculates contig statistics using 
#	fastaBLASTBP runs BLASTn, converts output to .bp, parses .bp to extract raw nucleotide hits 
sub processFile 
{
	my ($inputfa) = @_;
	`perl fastaStats.pl $inputfa`;
	`perl fastaBLASTBP.pl $inputfa ecoli_16s $geneFile`; #### specify your project basename in place of "ecoli_16s"
}
