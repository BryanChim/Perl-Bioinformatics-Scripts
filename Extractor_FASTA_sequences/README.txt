A basic Perl implementation of FASTA sequence extraction.

Program Flow:

1) Read in a tab delimited text file that contains, on each line:
<AccessionID>	<Start>	<End>

~ $seqID, EXAMPLE: 16S_sequences_to_be_extracted.txt


2) Attempts to find AccessionID matches in the input FASTA file

~ $inputFA, EXAMPLE: 16S_MICROBIOME_REF_FILE.FA


3) Once (if) a match is found for the ID, pulls out the specific
nucleotide sequence using substr() and writes it to an output FASTA file:
<Accession ID> | <Start> | <End>
<Sequence> 

~ $outputFA


Note: VERSION 2 (SeqExtraction2.pl) was redesigned for greater memory
efficiency. All of the ID/start/end's are stored in a 2-dimensional array.
FASTA file is read through only ONCE, instead of being repeatedly 
opened/closed (as in VERSION 1). All possible ID checks are performed 
whenever a new header is found.  