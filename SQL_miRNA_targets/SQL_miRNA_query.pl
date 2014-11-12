#!/usr/local/bin/perl

#### SQL_miRNA_query

#	 VERSION 1 (does not filter by minimum TOOL number)
#    PURPOSE: An exercise in querying an SQL database (MAPI_ts)
#		directly from within a Perl script. MAPI_ts stores
#		gene and miRNA target information in 4 tables:
#		- tbl_genes 
#		- tbl_predicted_targets
#		- tbl_expression
#		- tbl_unigene_refseq
#		** See accompanying README for more details
#		
#    INPUT: Query parameters, $miRNA and $tissue, via command-line shift
#		- SQL User name, from $username, via command-line shift
#		- SQL password, from $password, via command-line shift 
#
#    OUTPUT: From file $outputfile, via command-line shift
#		- Output is a SQL table:
#			- miRNA gene unigene_id type tissue chromosome ch_start ch_end description

############### LIBRARIES AND PRAGMAS ###############
use warnings;
use strict;
use DBI;

# command-line inputs for micro-RNA, tissue, output, user name and password
my $miRNA = shift;
my $tissue = shift;
my $outputfile = shift;
my $username = shift;
my $password = shift;

# variables
my @data;		#
my %geneHash;

# check for appropriate password input, otherwise die
if (!defined ($password)) 
{
	die "USAGE perl $0 need input for miRNA, tissue, output file, username and password\n";
}

# open output file
open OUT, ">$outputfile" or die "could not open output file\n";

# attempt to connect to MySQL database, "MAPI_ts" using the inputted username/password
my $db = DBI->connect("DBI:mysql:MAPI_ts", $username,$password) or
    die "Error in connecting database bnfo600\n";

# prepare and execute an SQL query	
my $sql = "select p.miRNA, u.refseq, u.unigene, e.type, e.tissue, g.chromosome, g.ch_start, g.ch_end, g.description, p.tool  from tbl_genes g, tbl_expression e, tbl_predicted_targets p, tbl_unigene_refseq u where u.unigene = e.gene and u.refseq = p.gene and u.refseq = g.refseq and p.miRNA = '$miRNA' and e.tissue = '$tissue'";
my $sql_statement = $db->prepare($sql);
$sql_statement->execute();

# print header line
print OUT "miRNA\tgene\tunigene_id\ttype\ttissue\tchromosome\tch_start\tch_end\tdescription\n";#tool
#            0      1        2         3       4    5       6         7          8           9

# fetch results of the above SQL query, print to output file
while(@data = $sql_statement->fetchrow_array()) 
{
    if ($data[0] eq "miRNA") 
	{
		next;
	}
	
    print OUT "@data\n";

}

close OUT;

