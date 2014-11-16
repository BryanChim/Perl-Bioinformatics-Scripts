#!/usr/local/bin/perl
# BryanChim_Linkage.pl
#
#	 ORIGINAL AUTHORS: Paul Fawcett and Jeff Elhai (Fall 2004)
#
###*** MODIFICATIONS BY: Bryan Chim -- Version 3 (Fall 2013) ***###############################
#		** Implemented Pearson Coefficient calculation for determining node distances (sub CalculatePearson)
#		** Implemented WEIGHTED Average Linkage Clustering functionality (sub AverageVectors)
#		** Implemented  both Single Linkage and Complete Linkage Clustering functionality 
#			(the BULK of the changes ~ sub GetArrayData and sub MakeNode)
#
###*** NOTE: Specify your desired linkage clustering type in the variable $clusteringtype (line 53)
#
#    PURPOSE: Clusters Microarray data using either:
#				- Weighted average linkage clustering 
#				- Single linkage clustering or 
#				- Complete linkage clustering
#             Output files are suitable for use with Michael Eisen's "Java TreeView"
#
#    INPUT FILES:
#             1. Tab Delimited log-transformed Spotted Microarray Ratio data
#                first column is gene ID, all subsequent represent microarray
#                ratio data.  The first row is assumed to be column headers--- ~ (base)$filename, line 55
#
#    OUTPUT FILES:
#             1. *.cdt file
#             2. *.gtr file

############## LIBRARIES AND PRAGMAS ################
use strict;
use warnings;

############## Subroutine Declarations ##############

sub CalculatePearson (\@\@);  # Calculates the similarity between two vectors

# note that the (\@\@) acts as a PROTOTYPE for the subroutine declaration
# specifying that we intend to pass the subroutine two separate
# references to lists.  This avoids the 'list collapse' that occurs
# when you try to directly pass two lists.  See page 262 of Cozens.

sub AverageVectors(\@\@$$); # Averages two vectors on an element-by-element basis

sub MakeNode ();       # joins all the genes and nodes

sub OrderTheGenes ($); # Depth-first search of the node tree to output genes    
                                                                                
sub GetArrayData;      # Reads the array data from an input file               
                                                                            
#################### CONSTANTS ######################

##*** Specify clustering type -- either "average" "single" or "complete" ***###                                                                    
my $clusteringtype = "complete"; ############################################## 

my $filename  = "BacillusData2";
my $inputfile = $filename . ".txt";
my $gtr_file  = ($filename . "_" . $clusteringtype . ".gtr");
my $cdt_file  = ($filename . "_" . $clusteringtype . ".cdt");
my $DUMMY     = -999;


my $weight;

my $genelist;
my $name      = 0;  ## the position of the "name" field in the master hash %nodehash
my $data      = 2;  ## the position of the reference to the gene expression vector in %nodehash
my $taken     = 1;  ## the position of the field marking if a gene or node is paired
my $right     = 4;  ## the position of the field containing the name of the right node
my $left      = 3;  ## the position of the field containing the name of the left node
my $correlate = 5;  ## the position of the field containing the correlation between right and left
if ($clusteringtype eq "average") {
         $weight = 6;
        }
if ($clusteringtype eq "single" || $clusteringtype eq "complete") {
         $genelist = 6;
        }

#################### VARIABLES ######################

my %similarity;       # master hash of node/gene and node/node correlations
my $columns;          # The number of arrays in the input file
my $rows;             # The number of genes in the input file
my $i;                # General purpose loop iterator
my %nodehash;         # The master hash used in version 2
my $mastercount = 0;  # for debugging - keep track of number of recursions

################### MAIN PROGRAM ####################

open CDT_FILE, ">$cdt_file" or die "Can't write $cdt_file: $!\n";
open GTR_FILE, ">$gtr_file" or die "Can't write $gtr_file: $!\n";

($rows, $columns) = GetArrayData($inputfile);

        ## input the microarray data
        # Although a shoddy practice, this subroutine also
        # silently populates the global variable @genelist.


print "$columns Microarrays were present in the input file.\n";
print "$rows Genes were present in the input file.\n";


MakeNode();  ## Not so much a true subroutine as much as a way to clarify the program logic

OrderTheGenes("NODE".($rows-1));
        # The first call of a recursive depth first search
        # that results in the nodes being printed out in the order required


close CDT_FILE;
close GTR_FILE;

print "All done!\n";

################################################################################
## GetArrayData:                                                               #
##                                                                             #
##         This function reads the tab-delimited gene expression data into     #
##         a multidimensional array called @genelist                           #
################################################################################

sub GetArrayData
        {
        my ($path) = @_;
        my $line;
        my @tempgene;
        my $i = 0;

        print "Parsing the input file...\n";

        open DATA_FILE, "<$path" or die "Can't open $path: $!\n";

        $line = <DATA_FILE>;        # Read first line (array names)

        while (defined $line)
                {
                chomp $line;             # Remove line break
                @tempgene = split(/\t/, $line); #break the line into tab-delimited chunks and put in a temporary list

                if ($i == 0)  ## if 0, we are reading the array names line, and this is special
                        {

                        print CDT_FILE "GID\t","CLID\t","NAME\t","GWEIGHT"; # for compatibility with Eisen formated .cdt file
                        for (my $j = 0; $j <scalar(@tempgene)-1; $j++)
                                {
                                print CDT_FILE "\t$tempgene[$j+1]";
                                }
                        print CDT_FILE "\nEWEIGHT\t\t\t";
                        for (my $j = 0; $j <scalar(@tempgene)-1; $j++)
                                {
                                print CDT_FILE "\t1";
                                }
                        print CDT_FILE "\n"; ## the whole mess above just writes the first two lines of the .cdt file i
                        }
                else
                        {
                        my $size     = @tempgene - 1;
                        my @templist = @tempgene[1..$size];
                        my @tempgenelist;                                                                              # for single or complete
                        push @tempgenelist, "GENE".$i;                                                          ######## clustering, each gene
                        $nodehash{"GENE".$i} = [$tempgene[$name], "unpaired", \@templist,"","","", 1];                 # pushes its own identity
                                                                                                                       # into a list, which will
                        if ($clusteringtype eq "single" || $clusteringtype eq "complete") {                            # later be collapsed with
                        $nodehash{"GENE".$i} = [$tempgene[$name], "unpaired", \@templist,"","","", \@tempgenelist];    # its future clusters
                        }
                        }

                $line = <DATA_FILE>;
                $i++;
                }

        close DATA_FILE;

        my $rows     = $i-1;
        my $columns = scalar( @{$nodehash{"GENE1"}->[$data]});

        return $rows, $columns;
        }


################################################################################
## CalculatePearson:                                                           #
##                                                                             #
##         This function calculates the Pearson Correlation Coefficient        #
##         between two gene expression vectors.                                #
##         the function assumes the lists are of equal size                    #
##         the zero th index is skipped, since it contains the gene name       #
################################################################################

sub CalculatePearson (\@\@)
{
	 my $deltas;
	 my $r;
	 my $xiyi = 0;
	 my $xi = 0;
	 my $yi = 0;
	 my $xi2 = 0;
	 my $yi2 = 0;
	(my $arrayX, my $arrayY) = @_;
	 my $size = scalar(@{$arrayX});

	for (my $i = 0; $i < $size; $i++)
	{
		$xiyi += (${$arrayX}[$i] * ${$arrayY}[$i]);
		$xi += ${$arrayX}[$i];
		$yi += ${$arrayY}[$i];
		$xi2 += (${$arrayX}[$i]**2);
		$yi2 += (${$arrayY}[$i]**2);
	}

	$r = (($size * $xiyi) - ($xi * $yi))/
	(sqrt(($size * $xi2) - $xi**2) *
	(sqrt(($size * $yi2) - $yi**2))
	);



 
	return $r;

}

################################################################################
## AverageVectors:                                                             #
##                                                                             #
##         This subroutine calculates the pairwise average                     #
##         of two gene expression vectors.                                     #
##         the function assumes the lists in the arguments are of equal size   #
################################################################################

sub AverageVectors (\@\@$$)    ## This is suitable for Average Linkage Clustering           ######### I left this (and all other
        {                    ## would need to be modified for complete linkage.                     # functions for average linkage
        my @returnlist;                                                                             # clustering) intact should you
                                                                                                    # wish to use it for comparison
       (my $arrayX, my $arrayY, my $weightX, my $weightY) = @_;

        for (my $i = 0; $i < scalar(@{$arrayX}); $i++)
                {
                $returnlist[$i] = ((${$arrayX}[$i]*$weightX) + (${$arrayY}[$i]*$weightY))
                                  /($weightX + $weightY);
                }
        return @returnlist;
        }


################################################################################
## MakeNode creates all of the nodes by joining genes to genes, genes to nodes,#
##          or nodes to nodes.  All nodes are recorded in the master hash      #
##          %nodehash.                                                         #
################################################################################

sub MakeNode()
{
	print "Building the initial similarity table\n";
	my $bestPearson;
	my $bestKeyLeft;
	my $bestKeyRight;

	 ### First we make an initial table of Gene similarities
	 ### Keep track of the best Pearson observed (and between which genes)

	 $bestPearson = -999;  ## Start off with an obvious bogus value

	foreach my $key (keys %nodehash)  ## iterate over all the genes
	{
		foreach (keys %nodehash)  ## compare each gene to all others
		{
			if ( ($_ eq $key) || (defined $similarity{$key.$_}))
			{
				next; ## Don't calculate self similarity or reciprocal keys
			}

			my $tempPearson = CalculatePearson
							(
							@{$nodehash{$key}->[$data]},
							@{$nodehash{$_}  ->[$data]}
							);

			$similarity{$key.$_} = $tempPearson;
			$similarity{$_.$key} = $tempPearson;

			if ($tempPearson >= $bestPearson)
			{
				$bestPearson  = $tempPearson;
				$bestKeyLeft  = $key;
				$bestKeyRight = $_;
			}
		}
	}

	## Now we have built an initial hash of similarities, that we need only
	## update when a new node is created. Now to the job of making nodes...

	my $newweight;

	for ($i = 1; $i < $rows; $i++) ## we must create ($rows-1) nodes.
	{
		print GTR_FILE "NODE$i","X\t$bestKeyLeft","X\t$bestKeyRight","X\t$bestPearson\n";
		print "Creating NODE$i","\ from $bestKeyLeft"," and $bestKeyRight "," correlation $bestPearson\n";

		## First, mark our bestKeyLeft and bestKeyRight nodes as "paired"
		$nodehash{$bestKeyRight}->[$taken] = "paired";
		$nodehash{$bestKeyLeft }->[$taken] = "paired";

		## Now create a new node based on the bestKeyLeft and bestKeyRight

####### if $clusteringtype chosen is "average," take the weighted average of the two vectors		
		if ($clusteringtype eq "average" ) 
		{  
			my @tempAverage = AverageVectors                                                         
						(                                                                         
						@{$nodehash{$bestKeyRight}->[$data]},
						@{$nodehash{$bestKeyLeft }->[$data]},
						$nodehash{$bestKeyRight}->[$weight],
						$nodehash{$bestKeyLeft }->[$weight]
						);
						
			my $newweight = $nodehash{$bestKeyRight}->[$weight] + $nodehash{$bestKeyLeft }->[$weight];
			$nodehash{"NODE".$i} =
			[
			"NODE".$i,      # Names really only needed for genes, but anyway..
			"unpaired",     # All new nodes are initially unpaired
			\@tempAverage,  # Remember, we must store a reference to a list, not a list
			$bestKeyLeft,   # It is was made from left...
			$bestKeyRight,  #    ...and right nodes
			$bestPearson,    # Keep track of the pearson for reporting later
			$newweight
			];
		}
		
####### if single or complete linkage is specified, collapse their gene lists together into one								
		if ($clusteringtype eq "single" || $clusteringtype eq "complete" ) 
		{   
		
			push @newgenelist, (@{$nodehash{$bestKeyRight}->[$genelist]},
								@{$nodehash{$bestKeyLeft}->[$genelist]}
							   );

			$nodehash{"NODE".$i} =
			[
			"NODE".$i,      # Names really only needed for genes, but anyway..
			"unpaired",     # All new nodes are initially unpaired
			"",  # Remember, we must store a reference to a list, not a list
			$bestKeyLeft,   # It is was made from left...
			$bestKeyRight,  #    ...and right nodes
			$bestPearson,    # Keep track of the pearson for reporting later
			\@newgenelist
			];
		}

		# Now we must update the similarity hash
		
####### If $clusteringtype chosen is "average" then all pairwise Pearson values must be recalculated
		if ($clusteringtype eq "average") 
		{                                  
			foreach (keys %nodehash)  ## iterate over all the extant nodes for the update                     
			{                                                                                         
				#print "working on key ", $_, "\n";

				if (  ($nodehash{$_}->[$taken] eq "paired")  #don't bother with paired nodes
				   || ($_ eq "NODE".$i)                      #don't bother with itself
				   || (defined($similarity{$_."NODE".$i}))   #don't bother with reciprocal keys
				   )
				{
					#print "skipped key ", $_, "\n";
					next;
				}
				
				else
				{
					my $tempPearson = CalculatePearson         # find Pearson between
							(                                  # just created nodes
							@{$nodehash{$_}->[$data]},         # and the remaining nodes/genes
							@{$nodehash{"NODE".$i}->[$data]}
							);

					$similarity{"NODE".$i.$_} = $tempPearson;  # record for both the normal
					$similarity{$_."NODE".$i} = $tempPearson;  # and reciprocal key
				}
			}
		}
			
############ If $clusteringtype chosen is "single" or "complete" then we must iterate through:
		#### 	- all keys of nodehash (outer for loop)
		#### 	- all genes in NODE$i's new gene list (middle for loop)
		#### 	- unpaired/nonidentical genes in all the gene lists of all nodes in nodehash (inner for loop)
		#### Then find the best Pearson out of all possible pairs between NODE$i and each of the other nodes
		#### * (single linkage - closest distance - highest Pearson)
		#### * (complete linkage - furthest distance - lowest Pearson)
		if (($clusteringtype eq "single") || ($clusteringtype eq "complete")) 
		{      
			### outer for loop																		  
			foreach (keys %nodehash) 
			{                                                   
				### set dummy values 																	  
				if ($clusteringtype eq "single") 
					{$bestPearson = -999;}                             
				if ($clusteringtype eq "complete") 
					{$bestPearson = 999;}                                         
				
				### middle for loop
				foreach my $geneInNODE(@{$nodehash{"NODE".$i}->[$genelist]}) 
				{   
					### inner for loop
					foreach  my $nodeGene (@{$nodehash{$_}->[$genelist]}) 
					{      

						if ( ($nodehash{$_}->[$taken] eq "paired")
								|| ($_ eq "NODE".$i)
								|| ($nodeGene eq $geneInNODE)
						   ) {next;}
										 
						my $tempPearson = $similarity{$nodeGene.$geneInNODE};

							
						if ($clusteringtype eq "single") 
						{
							if ($tempPearson > $bestPearson) 
								{$bestPearson = $tempPearson;}
						}
						
						if ($clusteringtype eq "complete") 
						{	if ($tempPearson < $bestPearson) 
								{$bestPearson = $tempPearson;}
						}

						$similarity{"NODE".$i.$_} = $bestPearson;
						$similarity{$_."NODE".$i} = $bestPearson;

					}
				}
			}
		}

		### OK, now we must iterate through the similarity hash
		### to find our next 'bestPearson', 'bestRightKey', & 'bestLeftKey' for
		### creating the next node.
		### Note that this is a little different than the way we did it in the
		### initial table, which contained only Genes, and no Nodes.

		$bestPearson = $DUMMY;  ## reinitialize with an obviously bogus value

		foreach (keys %similarity)  ## iterate over the similarity keys
		{

			/([A-Z]{4}[0-9]+)/;            ## recover the individual keys from $_
			my $keyLeft = $1;              ## Parenthesised regexp groups are stored in $1
			/((?<=[0-9])[A-Z]{4}[0-9]+)/;  ## a lookbehind regex to obtain keyRight
			my $keyRight = $1;

				#print "Now matching for $_\n";

			if (  $nodehash{$keyLeft }->[$taken] eq "paired"
			   || $nodehash{$keyRight}->[$taken] eq "paired" )
			{
				next;  #as usual, consider only unpaired combos
			}
		
			if ($similarity{$_} >= $bestPearson)
			{
			  #  print "Best Pearson reassigned to $_!!!\n";
				$bestPearson  = $similarity{$_};
				$bestKeyLeft =  $keyLeft;
				$bestKeyRight = $keyRight;
			  #  print "Best keys: L: $bestKeyLeft  R: $bestKeyRight\n";
			}
		}
    }
}


################################################################################
## OrderTheGenes:                                                              #
##                                                                             #
##   This function recursively works its way down a subtree of nodes using a   #
##   a depth-first search following the left-hand nodes first. When it runs    #
##   into a gene, it prints the gene and pops up a level.                      #
################################################################################

sub OrderTheGenes ($)      # get in the habit of always using a function prototype!
{

        (my $currentnode) = @_;

        print "Visiting ", $currentnode, "\n";

        if ($nodehash{$currentnode}->[$left] =~ /NODE\d+/)
		{
			OrderTheGenes($nodehash{$currentnode}->[$left])
        } else 
			{
				PrintNode($currentnode, $left);
			}

        if ($nodehash{$currentnode}->[$right] =~ /NODE\d+/) 
		{
			OrderTheGenes($nodehash{$currentnode}->[$right])
        } else 
			{		
				PrintNode($currentnode, $right);
			}

}

################################################################################
## PrintNode:                                                                  #
##    Simple subroutine to print out the information on a node or gene         #                                                                         #
##    $currentnode tells the starting node                                     #
##    $direction tells whether to print data associated with L or R gene       #
##    Format compatible with Eisen's 'TreeView'                                #
################################################################################

sub PrintNode
{
	(my $currentnode, my $direction) = @_;

	print CDT_FILE $nodehash{$currentnode}->[$direction], "X\t";   ## print GENEnumberX
	print CDT_FILE $nodehash{$nodehash{$currentnode}->[$direction]}->[$name], "\t"; #gene name duplicated
	print CDT_FILE $nodehash{$nodehash{$currentnode}->[$direction]}->[$name], "\t";  #
	print CDT_FILE "1";  # for Eisen's EWEIGHT column

	foreach (@{$nodehash{$nodehash{$currentnode}->[$direction]}->[$data]}) # my goodness!!
			{                                                              # Isn't Perl code lovely?
			print CDT_FILE "\t$_";  # Get at, print & tab-delimit each data item
			}
	print CDT_FILE "\n";
	return;
}
        

