#!/usr/local/bin/perl

#### BryanChim_Contig_Hit_Parse

#	 VERSION 1 - does not filter by e-value or number contig hits
#
#    PURPOSE: Extracts and outputs BLASTp hits
#		
#    INPUT: From file $filein, FHT 
#
#    OUTPUT: To standard output:
#		- contig ID
#		- gene name
#		- protein description
#		- evalue


$|++;
use warnings;
use strict;

my $eval_cutoff = "1e-50"; # <-- can be changed; minimum is currently 1e-200

# open filehandle for read-in
my $filein = "BryanChim_Contig_Hit_Parse";
open FHT, $filein or die "could not open file";

# declare variables to temporarily store regex captures
my $line;
my $gene;
my $evalue;
my $contig_id;
my $protein = "";
my $proteinline2;

# declare hashes to store relevant data for retrieval and output
my %geneHash;
my %contigHash;

# read and chomp
while ($line = <FHT>) {
        chomp $line;
        
        # regex to store contig names
        if ($line =~ /^Query= (VCU\d+contig\d+)/) {
                $contig_id = $1;
                }

        # regex to store e-value in 2nd level contigHash, corresponding to the
        # correct gene and contig
        if ($line =~ (/.*\[gene=(.*)\].*\.{3}\s+\d{3,4}\s+(.*)/)) {
                $gene = $1;
                $evalue = $2;
                $contigHash{$contig_id}{$gene} = $evalue;
                }

         # Regex to store protein description
         # checks if protein description spills over to next line
                #if so, employ another regex to capture and concatenate the rest
                # if not, remove ']' from description
         # store description in geneHash, corresponding to the current gene key
         #______________________________________________________________________
         if ($line =~ /^>\slcl.*\[gene=(.*)\]\s+\[protein=(.*)/){
                $gene = $1;
                $protein = $2;
                if (substr ($protein, -2, 1) eq ']'){
                chop $protein;
                chop $protein;
                } else {
                        $line = <FHT>;
                        if ($line =~ /^(.*?)\]/) {
                                $proteinline2 = $1;
                                $protein .= $proteinline2;
                                }
                        }
                $geneHash{$gene}{"protein_description"} = $protein;
                }
                
        }
        
# sort and print out info by nested foreach loop
# extract each contig from contigHash, and obtain matched genes
# get e-values from contigHash, using all stored combinations of contig and gene
# get protein-info from geneHash for each gene via key - protein_description
# check e-value -- print all info if it passes the cutoff
#_______________________________________________________________________________
my $hprotein;
foreach my $contig (sort keys %contigHash) {
        foreach my $hgene (sort keys %{$contigHash{$contig}}) {
          $evalue = $contigHash{$contig}{$hgene};
                if (defined($geneHash{$hgene}{"protein_description"})) {
                        $hprotein = $geneHash{$hgene}{"protein_description"};
                        }
                if (sprintf("%.200f", $evalue) < sprintf("%.200f", $eval_cutoff)){
                        print "$contig\t$hgene\t$hprotein\t$evalue\n";
                        }
                }
          }

close FHT;


                 


