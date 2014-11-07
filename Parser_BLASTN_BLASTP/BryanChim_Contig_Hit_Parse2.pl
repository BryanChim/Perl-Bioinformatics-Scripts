#!/usr/local/bin/perl
$|++;
use warnings;
use strict;

# e-value and contig-match-count cutoffs -- can be changed.
my $eval_cutoff = "1e-50"; # <-- current minimum is 1e-200
my $contigmatch_cutoff = 1; # <-- filter out genes that do not match MORE than this number of contigs

#open filehandle for read-in
my $filein = "BryanChim_Contig_Hit_Parse";
open FHT, $filein or die "could not open file";

# declare variables to temporarily store the relevant regex captures
my $line;
my $gene;
my $evalue;
my $contig_id;
my $protein = "";
my $proteinline2;

# declare hashes to store contig-matches, e-values and match-count
my %geneHash;
my %contig_matchcountHash;


# read and chomp file line by line
while ($line = <FHT>) {
        chomp $line;

        # regex to capture and store contig name
        if ($line =~ /^Query= (VCU\d+contig\d+)/) {
                $contig_id = $1;
                }

        # regex to capture and store gene names (as 1st level key in both hashes)
        # and e-values (stored in 2nd level of geneHash, keyed by contig_id)
        # match count is incremented for the respective gene in contig_matchcountHash
        # ONLY if it passes the user-specified cutoff value
        #_______________________________________________________________________
        if ($line =~ (/.*\[gene=(.*)\].*\.{3}\s+\d{3,4}\s+(.*)/)) {
                $gene = $1;
                $evalue = $2;
                $geneHash{$gene}{$contig_id} = $evalue;
                if (sprintf("%.200f", $evalue) < sprintf("%.200f", $eval_cutoff)){
                        if (defined($contig_matchcountHash{$gene})) {
                                $contig_matchcountHash{$gene}++;
                                } else {
                                        $contig_matchcountHash{$gene} = 1;
                                        }
                                }
                }
        }

# Check and print out values via nested foreach loops
# First foreach loop extracts gene names from geneHash, runs them through
# contig_matchcountHash to acquire match count values
# if match count value exceeds the cutoff, print out the gene name, then open up
# 2nd level of geneHash to acquire and print contigs and the respective e-values
#_______________________________________________________________________________
foreach my $hgene (sort keys %geneHash) {
        if (defined($contig_matchcountHash{$hgene})) {
                if (($contig_matchcountHash{$hgene}) > $contigmatch_cutoff) {
                print "\n$hgene\t ";
                foreach my $hcontig (sort keys %{$geneHash{$hgene}}) {
                        $evalue = $geneHash{$hgene}{$hcontig};

                                print "$hcontig\t$evalue\t";
                                }
                        }
                }
        }

close FHT;

