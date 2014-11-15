An exercise in building and using a Position-Specific-Scoring Matrix.

* FindMotif_mod.pl accepts a training set of motifs 
~ $motif_input (EXAMPLE: 71NpNtSm.txt) and then, for each position:
-- counts occurrence of each nucleotide BY POSITION
-- calculates frequency of each nucleotide at each position
-- calculates information content at each position
-- applies a pseudocount to adjust for any "non-occurrences"
-- log-transforms the frequencies to generate a summable scoring matrix

* The program then accepts an inputted sequence ~ $sequence_file (EXAMPLE: lef.txt), 
runs a sliding window through it of same length as our PSSM, and scores each window.

** For each of these windows, ONLY the positions deemed "informational" will be
scored and considered in the running total. 
	Hits are outputted ~ $motif_hits (EXAMPLE: motif.hits)

** You may adjust certain constants (near top of program, ~ line 50) to customize your PSSM
	- $B:  total pseudocounts to be distributed when adjusting scores
	- $info_threshold:  information threshold for a window position to be scored
	- $threshold: threshold for minimum hit score
	- $number_to_save: number of top scoring hits to save to output
	- %bg_frequency: background frequencies used in calculating log-odds ratios
		* current frequencies are based on Anabaena intergenic sequences