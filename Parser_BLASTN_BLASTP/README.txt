A simple BLAST output parser that filters BLASTP hits

Version 2 filters by a minimum e-value and match count

- Would need to modify the $filein variable to match the desired BLASTp input
- Would need to modify the regular expression that recognizes the contig header line 
  (will depend on BLAST version and your custom naming convention0
