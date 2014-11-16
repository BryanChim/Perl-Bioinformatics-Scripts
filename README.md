Perl-Bioinformatics-Scripts
===========================
A collection of Perl programs written over the past few years for bioinformatics classes and/or work.

- **BLASTP_homegrown** - a basic interpretation and implementation of NCBI BlastP that outputs alignments between 2 input protein sequences

- **Cluster_Analysis** - performs weighted average, single, or complete linkage clustering on an input microarray dataset and outputs .cdt and .gtr files that are compatible with Java TreeView

- **Extractor_FASTA_sequences** - given an input file containing lines of <Accession ID> <Start> <End>, attempts to search an input FASTA file to find those sequences by <Accession ID>, and output truncated sequences from their <Start> to <End>

- **Markov_Modeling** - using an input training set of DNA sequences, generates a Markov Model of user-specified order for either nucleotide or amino acid transition states and then scores an input DNA or protein sequence

- **Parser_BLASTP_hits** - a simple parser that extracts hits from NCBI BlastP output.. version 2 can filter hits based on maximum e-value or minimum match count

- **Parser_BLAST_STATS_BP** - an extensive BLAST pipeline that performs a series of tasks including: calculation of FASTA statistics and n50, BLASTing of an input query sequence against an input target FASTA, and extraction of resultant hits from original input target FASTA

- **Position_Specific_Scoring_Matrix** - using an input training set of consensus DNA sequences, constructs a PSSM that is used to score each possible window within an input sequence

- **SQL_miRNA_targets** - logs into a MySQL miRNA target database, runs a query, fetches the resultant table, and parses it to determine the most likely miRNA targets for a user-specified miRNA within a user-specified tissue

- **Parser_BLASTn_allhits_nohits.pl** - parses through BlastN output and extracts two sets of data: one comprised of information on every single query that yielded a hit against the target/subject, and another comprised of information on every single query that yielded NO hits against the target/subject

- **Parser_fpkm_tracking.pl** - parses through the *.fpkm_tracking output file that is generated from Cufflinks, and extracts/calculates some basic statistics (ie. total # genes, smallest/largest/average fpkm, most/least expressed genes etc.)
