#!/usr/local/bin/perl
use warnings;
use strict;
use DBI;

my $miRNA = shift;
my $tissue = shift;
my $outputfile = shift;
my $username = shift;
my $password = shift;
if (!defined ($password)) {
        die "USAGE perl $0 need input for miRNA, tissue, output file, username and password\n";
        }

open OUT, ">$outputfile" or die "could not open output file\n";



my $db = DBI->connect("DBI:mysql:MAPI_ts", $username,$password) or
    die "Error in connecting database MAPI_ts\n";

my $sql = "select p.miRNA, u.refseq, u.unigene, e.type, e.tissue, g.chromosome, g.ch_start, g.ch_end, g.description, p.tool  from tbl_genes g, tbl_expression e, tbl_predicted_targets p, tbl_unigene_refseq u where u.unigene = e.gene and u.refseq = p.gene and u.refseq = g.refseq and p.miRNA = '$miRNA' and e.tissue = '$tissue';";

my $sql_statement = $db->prepare($sql);

$sql_statement->execute();

print OUT "miRNA\tgene\tunigene_id\ttype\ttissue\tchromosome\tch_start\tch_end\tdescription\n";#tool
#  indices  0      1        2         3       4       5         6         7          8            9



my %geneHash;
my %toolHash;
my @data;

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


my $i = 0;
my $omiRNA;
my $unigene;
my $type;
my $otissue;
my $chromosome;
my $ch_start;
my $ch_end;
my $description;
my @tools;
my $toolcount;

foreach my $gene (keys %toolHash) {


        if ((scalar (keys %{$toolHash{$gene}})) >= 3 && ($i < 10)) {

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
                        foreach my $tool (@tools) {
                        print OUT "$tool, ";
                        }
                        print OUT "\tTool Count: $toolcount\n\n";

                 $i++;
                 @tools = ();
                }

        }

close OUT;

