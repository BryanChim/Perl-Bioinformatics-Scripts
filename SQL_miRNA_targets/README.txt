Program finds 10 predicted miRNA targets for requested miRNA by aleast three prediction tools. 
Command-line input and output parameters:
    1) miRNA = miRNA of interest (ex. hsa-miR-17-3P)
    2) tissue =  tissue of interest  (ex. liver tumor) 
    3) Output file
    4) username = mysql username
    5) password = mysql password

This program will generate  the output below in tab-delimited format.
- miRNA
- gene
- unigene id
- type
- tissue
- gene chromosome
- gene chromosome start
- gene chromosome end
- gene description

Print only those genes which are predicted to be targetted by atleast 3  prediction tools. 
use: 
tbl_expression (gene, type, tissue, expression), 
tbl_predicted_targets (miRNA, gene, score, tool),
tbl_unigene_refseq (unigene, refseq) 
tbl_genes (refseq, contig, chromosome, ch_start, ch_end, direction, symbol, description)
* unigene = gene in expression
* refseq = gene in predicted_targets, refseq in gene, 

** join expression and refseq via unigene
** join predicted_targets and refseq via refseq

mysql> select p.miRNA, u.refseq, u.unigene, e.type, e.tissue, g.chromosome, g.ch_start, g.ch_end, g.description, p.tool  from tbl_genes g, tbl
_expression e, tbl_predicted_targets p, tbl_unigene_refseq u where u.unigene = e.gene and u.refseq = p.gene and u.refseq = g.refseq and p.miRN
A = "hsa-miR-17-3P" and e.tissue = "fetus" limit 10;
