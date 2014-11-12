#!/usr/local/bin/perl

#### fastaBLASTBP.pl

#    PURPOSE: BLASTn pipeline-- runs blast on a given query against a given target
#		- Extracts condensed hit info into a .bp file, with useParse.pl
#		- Extracts each sequence ID, start and end from .bp file
#		- Using the ID, start and end, extracts the hit sequence from original query
#		- Outputs a custom-generated FASTA output
#
#    INPUT: NOTE -- to be passed in directly via fastaMAIN.pl
#			$assemblyFA - input used to create the BLAST database
#			$projectname - takes the query gene's FASTA file
#			$geneFA - query sequence to BLASTed against $assemblyFA
#
#    OUTPUT: From $outputFA, generated as "$assemblyFA_$projectname_extracted.fa"
#			Rows will contain: 
#			>FASTA base-name|seq ID|seq length|seq start|seq end|sequence


# libraries and pragmas
$|++;
use warnings;
use strict;
use File::Basename;

# command-line variables (passed directly via fastaMAIN.pl)
my $assemblyFA = shift;
my $projectname = shift;
my $geneFA = shift;

# check for existence of $geneFA -- die if not found
if(!defined($geneFA)) 
{die "USAGE: perl $0 need assembly fasta, project and gene fasta file names\n";}


# variable declaration/initialization
my ($line, $sID, $sStart, $sEnd, $seqlength);
my $FAbasename = basename($assemblyFA);
my $subline;
my $bp = "BlastOutput.bp";
my $outputFA = $assemblyFA . "_" . $projectname . "_extracted.fa";

#### BLAST commands to generate database from $assemblyFA 
# BLAST a particular gene, $geneFA, against the database with e-value threshold of 1e15
# run useParse.pl on the output to generate .bp file 
`formatdb -i $assemblyFA -p F`;
`blastall -d $assemblyFA -i $geneFA -o BlastOutput.blast -p blastn -e 1e15`;
`perl useParse.pl BlastOutput.blast`;

# open .bp, assembly FA, and output file-handles
open BP, $bp or die "could not open bp file\n";
open AF, $assemblyFA or die "could not open assembly fasta file\n";
open OFA, ">$outputFA" or die "could not open output file\n";

# parse .bp file
while ($line = <BP>) 
{
	chomp $line;

	# get sequence ID
	if ($line =~ /SUB:(\w+.*?)\s*/) 
	{
			$sID = $1;
	}

	# get sequence start and end
	if ($line =~ /.*sb:(\d+)\sse:(\d+)\s.*/) 
	{
			($sStart, $sEnd) = ($1, $2);
			
			# find the sequence ID in the original assembly FA
			while ($line = <AF>) 
			{
				chomp $line;
				if ($line =~ /^>$sID\s.*/) 
				{
					$line = <AF>;
					chomp $line;
					
					# extract the substring of the nucleotide sequence from the hit's start to end		
					if ($sEnd > $sStart) 
					{
							$subline = substr($line, $sStart - 1, $sEnd - $sStart + 1);
					} 
					
					else 
					
					# if start coordinate is greater than end, hit is on reverse strand
					{	
						$subline = substr($line, $sEnd - 1, $sStart - $sEnd + 1);
						$subline = reverse($subline);
					}
					
					# finally, get length of the hit and print data to output file
					$seqlength = length($subline);
					print OFA ">$FAbasename|$sID|$seqlength|$sStart|$sEnd\n$subline\n";
				}
			}
	}
	
	# close and reopen assembly file for next hit
	close AF;
	open AF, $assemblyFA or die "could not reopen assembly fasta file\n";

}

# end of program, close all file-handles
close AF;
close BP;
close OFA;


