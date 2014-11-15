A perl program designed to log in to, authenticate and directly query a MySQL Database, MAPI_ts,
in order to find predicted miRNA gene targets.

* Command-line input and output parameters:
    1) miRNA = miRNA of interest (ex. hsa-miR-17-3P)
    2) tissue =  tissue of interest  (ex. liver tumor) 
    3) output file
    4) username = mysql username
    5) password = mysql password

* The database in question is comprised of the following tables (irrelevant or unused tables omitted):
tbl_expression (gene, type, tissue, expression), 
tbl_predicted_targets (miRNA, gene, score, tool),
tbl_genes (refseq, contig, chromosome, ch_start, ch_end, direction, symbol, description),
tbl_unigene_refseq (unigene, refseq) 
    * where unigene = gene in tbl_expression
    * and where refseq = gene in tbl_predicted_targets 

* The program is currently set up to find top 10 predicted miRNA targets for the requested miRNA 
    BY AT LEAST 3 PREDICTION TOOLS ~ accessed and counted via joins with tbl_predicted_targets 

* Tab-Delimited Output:
- miRNA
- gene
- unigene id
- type
- tissue
- gene chromosome
- gene chromosome start
- gene chromosome end
- gene description

* Example of a MySQL Query that could be called within the script:
mysql> select p.miRNA, u.refseq, u.unigene, e.type, e.tissue, g.chromosome, g.ch_start, g.ch_end, g.description, p.tool  
> from tbl_genes g, tbl_expression e, tbl_predicted_targets p, tbl_unigene_refseq u 
> where u.unigene = e.gene and u.refseq = p.gene and u.refseq = g.refseq and p.miRNA = "hsa-miR-17-3P" and e.tissue = "fetus" limit 10;
