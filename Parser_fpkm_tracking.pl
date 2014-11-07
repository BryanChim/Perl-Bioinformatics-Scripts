#!/usr/local/bin/perl
$|++;
use warnings;
use strict;

####################### FILES #######################

my $file_gene_input = shift;
my $file_gene_output = shift;

if (!defined($file_gene_output))
        {die "USAGE: perl $0 need file-gene-input and file-gene-output";}

open FGI, $file_gene_input or die "could not open file_gene_input";
open FGO, ">$file_gene_output" or die "could not open file_gene_output";

#################### VARIABLES ######################
 my $line;				# line read variable
 my $chr;				# chromosome number
 my $id;				# gene id
 my $fpkm;				# gene FPKM value
 my $total_genes;		# gene counter
 my $totalnz_genes;		# counter for non-zero read-count genes
 my $minFPKM;			# tracker variable for smallest FPKM value
 my $maxFPKM;			# tracker variable for largest FPKM value
 my $avgFPKM;			# average FPKM
 my $minGene;			# tracker variable for gene with smallest FPKM value
 my $maxGene;			# tracker variable for gene with largest FPKM value
 my $printstring;		# holds joined data for each chromosome
 my @cols;				# holds file data by column
 my %chrom;				# hash for storing and updating chromosome data
 
 <FGI>;

# read in fpkm_tracking file
while ($line = <FGI>) 
{
	# chomp and split line,
	chomp $line;
	@cols = split("\t", $line);
	
	# get chromosome number
	if ($line =~ /(\d+):/) 
	{
		$chr = $1;
	}
	
	# get gene ID and FPKM value
	$id = $cols[0];
	$fpkm = $cols[10];
	
	# initialize DUMMY values for max and min FPKM
	if (!defined $chrom{$chr}{"maxFPKM"}) 
	{
		$chrom{$chr}{"maxFPKM"} = 0;
	}
	
	if (!defined $chrom{$chr}{"minFPKM"}) 
	{
		$chrom{$chr}{"minFPKM"} = 1000000;
	}
	
	# add to FPKM sum if non-zero, and increment count for non-zero FPKM genes  
	if (!($fpkm eq '0')) 
	{
		$chrom{$chr}{"fpkmsum"} += $fpkm;
		$chrom{$chr}{"nonzero_fpkmcount"}++;
		
		#update mix/max FPKM
		if ($fpkm < $chrom{$chr}{"minFPKM"}) 
		{
			$chrom{$chr}{"minFPKM"} = $fpkm;
			$chrom{$chr}{"minGene"} = $id;
		}
		
		if ($fpkm > $chrom{$chr}{"maxFPKM"}) 
		{
			$chrom{$chr}{"maxFPKM"} = $fpkm;
			$chrom{$chr}{"maxGene"} = $id;
		}
	}
	
	# increment gene counter for current chromosome
	$chrom{$chr}{"idcount"}++;
}
        
		
# print header line
print FGO "Chr\t#Genes\t#Non-zero Genes\tMin FPKM\tMax FPKM\tAvg FPKM\tMin Gene\tMax Gene\n";
   
# query %chrom hash and print out data from lowest to highest chromosome number   
foreach my $chrnum (sort {$a <=> $b} keys %chrom) 
{        
	$total_genes = $chrom{$chrnum}{"idcount"};
	$totalnz_genes = $chrom{$chrnum}{"nonzero_fpkmcount"};
	$minFPKM = $chrom{$chrnum}{"minFPKM"};
	$maxFPKM = $chrom{$chrnum}{"maxFPKM"};
	$avgFPKM = $chrom{$chrnum}{"fpkmsum"}/$chrom{$chrnum}{"nonzero_fpkmcount"};
	$minGene = $chrom{$chrnum}{"minGene"};
	$maxGene = $chrom{$chrnum}{"maxGene"};

	$printstring = join("\t", $chrnum,$total_genes,$totalnz_genes,$minFPKM,$maxFPKM,$avgFPKM,$minGene,$maxGene);
		 
	print FGO "$printstring\n";
             
}

# close filehandles	
close FGI;
close FGO;
exit;

