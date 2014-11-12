#!/usr/local/bin/perl

#### SQL_miRNA_plusTOOLS_query

#	 VERSION 2 (filters by minimum TOOL number)
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
#			- miRNA gene unigene_id type tissue chromosome ch_start ch_end description, tool

############### LIBRARIES AND PRAGMAS ###############
use warnings;
use strict;
use DBI;

################ COMMAND-LINE INPUTS ################
my $miRNA = shift;
my $tissue = shift;
my $outputfile = shift;
my $username = shift;
my $password = shift;

# check for appropriate password input, otherwise die
if (!defined ($password)) {
        die "USAGE perl $0 need input for miRNA, tissue, output file, username and password\n";
        }

#################### VARIABLES ######################

my $i = 0;			# counter for number of targets to output
my $omiRNA;			# miRNA ID
my $unigene;		# unigene ID
my $type;			# gene type
my $otissue;		# tissue type
my $chromosome;		# chromosome number
my $ch_start;		# chromosome coordinate start
my $ch_end;			# chromosome coordinate end
my $description;	# gene description
my $toolcount;		# tool counter
my @tools;			# array for tool names for a given gene
my @data;			# stores results of running fetchrow_array() on the SQL query
my %geneHash;		# stores gene information acquired from SQL query
my %toolHash;		# stores tools for a given gene

###################### FILES ########################
open OUT, ">$outputfile" or die "could not open output file\n";

################## MAIN PROGRAM #####################

# attempt to connect to MySQL database, "MAPI_ts" using the inputted username/password
my $db = DBI->connect("DBI:mysql:MAPI_ts", $username,$password) or
    die "Error in connecting database MAPI_ts\n";

# prepare and execute an SQL query
my $sql = "select p.miRNA, u.refseq, u.unigene, e.type, e.tissue, g.chromosome, g.ch_start, g.ch_end, g.description, p.tool  from tbl_genes g, tbl_expression e, tbl_predicted_targets p, tbl_unigene_refseq u where u.unigene = e.gene and u.refseq = p.gene and u.refseq = g.refseq and p.miRNA = '$miRNA' and e.tissue = '$tissue';";
my $sql_statement = $db->prepare($sql);
$sql_statement->execute();

# print header line
print OUT "miRNA\tgene\tunigene_id\ttype\ttissue\tchromosome\tch_start\tch_end\tdescription\n";#tool
#  indices  0      1        2         3       4       5         6         7          8            9


# fetch results of the above SQL query, store in %geneHash, with gene as key
# collect tools that claim to predict miRNA for given gene
while(@data = $sql_statement->fetchrow_array()) {

    if ($data[0] eq "miRNA") {next;}

    $geneHash{$data[1]}{miRNA} = $data[0];
    $geneHash{$data[1]}{unigene} = $data[2];
    $geneHash{$data[1]}{type} = $data[3];
    $geneHash{$data[1]}{tissue} = $data[4];
    $geneHash{$data[1]}{chromosome} = $data[5];
    $geneHash{$data[1]}{ch_start} = $data[6];
    $geneHash{$data[1]}{ch_end} = $data[7];
    
    if (!defined ($data[8])) {
        $geneHash{$data[1]}{description} = "NULL";
        } else {
        $geneHash{$data[1]}{description} = $data[8];
        };
        
    $toolHash{$data[1]}{$data[9]} = $data[9];
    }


# iterate through genes that are predicted miRNA targets
# if number of tools that corroborate on prediction meets
# or exceeds 3, output gene information
foreach my $gene (keys %toolHash) 
{
	if ((scalar (keys %{$toolHash{$gene}})) >= 3 && ($i < 10)) 
	{
			$omiRNA = $geneHash{$gene}{miRNA};
			$unigene = $geneHash{$gene}{unigene};
			$type = $geneHash{$gene}{type};
			$otissue = $geneHash{$gene}{tissue};
			$chromosome = $geneHash{$gene}{chromosome};
			$ch_start = $geneHash{$gene}{ch_start};
			$ch_end = $geneHash{$gene}{ch_end};
			$description = $geneHash{$gene}{description};
			push @tools, (keys %{$toolHash{$gene}});
			$toolcount = scalar(keys %{$toolHash{$gene}});
			
			print OUT "$omiRNA\t$gene\t$unigene\t$type\t$otissue\t$chromosome\t$ch_start\t$ch_end\t$description\n";
			print OUT "Tools: ";
			
			foreach my $tool (@tools) 
			{
				print OUT "$tool, ";
			}
			
			print OUT "\tTool Count: $toolcount\n\n";

			 $i++;
			 @tools = ();
	}

}

close OUT;

