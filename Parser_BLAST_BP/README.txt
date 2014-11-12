BLASTn pipeline that accepts:
- a single of directory of FASTA files, with .fsa file extension (1st arg)
- a single FASTA file containing the sequence of some gene of interest (2nd arg)
- calculates statistics and n50 metric for each FASTA file (outputs to .fsa_stats)
- BLASTs gene of interest against each FASTA, extracts and outputs raws hit sequences (outputs to _extracted.fa)

* Pipeline begins by runnning fastaMAIN.pl with above arguments
* fastaMAIN.pl collects input FASTA files and further calls fastaSTATS.pl and fastaBLASTBP.pl to complete the pipeline

** See comments and documentation in fastaMAIN.pl header for more information
