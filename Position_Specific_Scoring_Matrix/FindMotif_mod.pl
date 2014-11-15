#!/usr/bin/perl -w

# FINDMOTIF_mod
#
#    ORIGINAL AUTHORS: Jeff Elhai and Paul Fawcett (Fall 2003)
#
#	 MODIFICATIONS BY: Bryan Chim (Fall 2013 - VERSION 5)
#		** SEE "DIFFERS FROM VERSION 4" section BELOW
#
#    PURPOSE: Calculates position-specific scoring matrix
#		       from aligned sequence
#             Finds sequences in genome sequence that PSSM
#                fits best
#
#    INPUT FILES:
#             motif_input: input file of aligned motifs
#                          alignments may have gaps, represented by "."
#             sequence_file: input file containing sequences in FA format
#                          May contain multiple contigs
#
#    OUTPUT FILE: Description of format of output files use strict;
#             motif_hits: output file of best motifs within sequence file
#                Each line is name of sequence, score, start_coord, sequence
#
#    DIFFERS FROM VERSION 3 in that the subroutine Calc_motif_frequency has been split
#             into two subroutines to facilitate later modification of the program.
#             The two are:
#                 Calc_nucleotide_counts
#                 Calc_nucleotide_frequencies
#
#    DIFFERS FROM VERSION 4 in that information content is calculated for each position in the 
#			motif training set ~ $motif_input
#		* Functionality is combined with the subroutine CALC_NUCLEOTIDE_COUNTS () into:
#			CALC_NUCLEOTIDE_COUNTS_ANDINFO (motif sequences) - LINES 230-283
#		* If information content exceeds a user-defined threshold ~ $info_threshold
#			then it gets scored in the score_of_region() subroutine. 
#			(other positions are ignored)
#
# REMOVE FOLLOWING LINE IN READ_CONTIGS AFTER TESTING !!!!
#
#      if (length($sequence) < 500000) {$line=""}  ### TESTING ONLY! ###
#
# This line causes only the first 50,000 nucleotides of the genomic sequence
# to be read, greatly speeding up execution time. But of course the results
# of searches are no longer valid.

############## LIBRARIES AND PRAGMAS ################
  use strict;

#################### CONSTANTS ######################

  my $true = 1;
  my $false = 0;
  my $LF = "\n";            # Line feed
  my $B = 0.1;              # Weighting factor for pseudocounts
  my $info_threshold = 1.2;    # Threshold used when deciding whether or not push into @informational  
  my $threshold = -2**30;   # Save motifs with scores greater than $threshold
                            #   Initialize to very small number so that all
                            #   scores initially are saved
  my $number_to_save = 10; # Number of top scores and coordinates to save
  my %bg_frequency =        # Nucleotide frequencies of overall sequences
     ("A" => 0.33264,       # Taken from Anabaena intergenic sequences
      "T" => 0.33264,
      "G" => 0.16736,
      "C" => 0.16736);

#################### VARIABLES ######################

  my $motif_length;         # Length of motifs
  my @motif_sequences;      # All motifs read
  my @contigs;              # All contigs read
  my $contig;               # Current contig
  my $contig_index;         # Index of current contig
  our @headers;             # All contig headers read
  my %counts;               # {nucleotide}[position in motif]
                            #    Number of occurrences of nucleotides at each position

  my @informational;        # Those positions in the motif that contribute
                            #    to the motif -- i.e. those that are not
                            #    don't-care positions.
  my %log_q;                # $log_q{$nucleotide}[$position] is the log of the
                            #    frequency of the nucleotide at that position
  my %log_p;                # $log_p{$nucleotide} is the background frequency
                            #    of the nucleotide
  my @hit_info;             # Information about each high-scoring sequence

###################### FILES ########################
#    Aligned sequences in the file defined by $motif_input are analyzed.
#    The resulting PSSM is used to search through the file defined by $sequence_file.
#    The sequence fragments within the large sequence that best match the PSSM
#       are output to the file defined by $motif_hits

  my $motif_input = '71NpntSm.txt';
  my $sequence_file = 'lef.txt';
  my $motif_hits = 'motif.hits';

################### MAIN PROGRAM ####################
#    Read motifs
#    Calculate PSSM
#       Calculate counts of nucleotides in aligned motifs
#       Determine informational columns
#       Calculate frequencies (modified by pseudocounts)
#    Read contigs (or single sequence)
#    Use PSSM to score window running through contigs
#    Print top hits

  print "Getting motifs...", $LF;
  @motif_sequences = Get_motif($motif_input);

  print "Making PSSM...", $LF;

  @informational = Calc_nucleotide_counts_andINFO (@motif_sequences);
  Calc_nucleotide_frequencies(@motif_sequences);
  Calc_background_frequencies();
  
     my $N = Size(@motif_sequences);
     foreach my $nucleotide ("A", "C", "G", "T") {
        foreach my $position (@informational) {
           my $count = ($counts{$nucleotide}[$position] or 0);
           $log_q{$nucleotide}[$position] =
              log(($count + $B * $bg_frequency{$nucleotide}) / ($N + $B));
		   my $score = $log_q{$nucleotide}[$position] - $log_p{$nucleotide};
		   print "$position\t$nucleotide\t$score\n";
		   
        }
		
	print "\n";
     }
  
  print "Reading contigs...", $LF;
  @contigs = Read_contigs($sequence_file);

  print "Analyzing genome...", $LF;
  foreach $contig_index (0 .. Size(@contigs)-1) {
     Analyze_sequence($headers[$contig_index], $contigs[$contig_index]);
  }

  Print_hits();

#################### SUBROUTINES ####################

         #### INPUT/OUTPUT SUBROUTINES ####

#### READ_CONTIGS
#    Reads file of sequences in FastA format
#    May be more than one sequence per file
#    Saves header in global variable @header
#    Returns array of sequences (all upper case)
#    Saves total length of sequences read in global variable $total_length
  sub Read_contigs {
     my ($file_path) = @_;
     open GENOME, "<$file_path" or die "Can't open $file_path: $!\n";

     our @headers;              # Stores FastA header in global variable
                                # Must be declared in main program to be accessed
     our $total_length = 0;     # Stores number of nucleotides read in global variable
                                # Must be declared in main program to be accessed
     my $sequence;              # Individual sequence
     my @sequences;             # Stores all sequences

     my $line = <GENOME>;
     while (defined $line) {
        chomp $line;
        if (substr($line,0,1) eq ">") {      # Found header
           push @headers, $line;             # Save header
           if (defined $sequence) {          # If sequence exists, save it
              push @sequences, $sequence;
              $total_length = $total_length + length($sequence);
           }
           $sequence = "";                   # Initialize new sequence
        } else {                             # More sequence
           if (length($sequence) > 50000) {$line=""}  ### TESTING ONLY! ###
           $sequence .= uc($line);           # Append line to sequence
        }
        $line = <GENOME>;
     }
     push @sequences, $sequence;             # Save last sequence
     return @sequences;
  }

#### GET_MOTIF
#    Return array of aligned sequences
#    Set the global variables $motif_length
#    Quality check: All sequences must be the same length
  sub Get_motif {
     my ($motif_file_path) = @_;
     my @motif_sequences;
     my $line;
     my $first_time = $true;
     my $motif;
     my $first_motif;

     open MOTIF_INPUT, "<$motif_file_path" or die "Can't open $motif_file_path: $!\n";

     $line = <MOTIF_INPUT>;
     while (defined $line) {
        if ($line =~ /\b([ACGTacgt]+)\b/) {
           $motif = uc($1);
           if ($first_time) {
              $first_motif = $motif;
              $motif_length = length($first_motif);
              $first_time = $false;
           }
           if (length($motif) == $motif_length) {
              push @motif_sequences, $motif;
           } else {
              die "Length of motif $motif in $motif_input is not $motif_length $LF";
           }
        }
        $line = <MOTIF_INPUT>;
     }
     close MOTIF_INPUT;
     return (@motif_sequences);
  }

#### PRINT_HITS
#    Sort hits for final time
#    Go through each hit from best to worst
#    Formatted print of each hit to screen and to file
  sub Print_hits {
     sort_and_discard_excess_values();
     open HITS, ">$motif_hits" or die "Can't open $motif_hits: $!\n";

     foreach my $info (@hit_info) {
        my ($score, $name, $start, $subregion) = @$info;
        printf "%9s %6.1f %6d %s\n", $name, $score, $start, $subregion;
        printf HITS "%9s %6.1f %6d %s\n", $name, $score, $start, $subregion;
     }
     close HITS
  }
       #### PSSM CONSTRUCTION SUBROUTINES ####

#### CALC_NUCLEOTIDE_COUNTS_ANDINFO (motif sequences)
#    Accumulate nucleotide totals at each position
#    counts{nucleotide}[position] is the number of times the nucleotide appears
#          at the position
# 	 Calculates information content at each position in the motifs
#    Pushes sufficiently "informational" positions to @informational
 sub Calc_nucleotide_counts_andINFO 
 {
     my (@motif_sequences) = @_;
     my $motif;
     my $position;
     my $nucleotide;
     my $count;
     my $N = Size(@motif_sequences);
     my @informational;
     my $pseudoct_score;
	 my $prob;
     my $uncertainty = 0;
     my $information_content = 0;

     foreach $motif (@motif_sequences) {
        foreach $position (0 .. $motif_length - 1) {
           $nucleotide = substr($motif, $position, 1);
           if (defined $counts{$nucleotide}[$position]) {
              $counts{$nucleotide}[$position] = $counts{$nucleotide}[$position] + 1;
           } else {
              $counts{$nucleotide}[$position] = 1;
           }
        }
     }
 
     foreach $position (0 .. $motif_length - 1) 
	 {
        foreach $nucleotide ('A', 'G', 'C', 'T') 
		{

           if (defined (($counts{$nucleotide}[$position])))
		   {
			   $count = ($counts{$nucleotide}[$position]);
			   $prob = $count/$N;
			   #$pseudoct_score = (($count + $B * $bg_frequency{$nucleotide}) / ($N + $B));
			   $uncertainty += -($prob * (log($prob)/log(2)));          
           }
        }
		
			$information_content = 2 - $uncertainty;
			print "$position\t$information_content\n";
			
			if ($information_content > $info_threshold) 
			{
				push (@informational, $position);
			}
			
      $uncertainty = 0;
     }

     print "@informational\n";
     return @informational;

}

#### CALC_NUCLEOTIDE_FREQUENCIES (motif sequences)
#    Calculates the adjusted frequency at each position.
#
#    The actual frequency is
#
#               counts{nucleotide}[position] / N   ,
#    where:
#       counts{nucleotide}[position] is the number of times the nucleotide appears
#          at the position
#       N is the number of motifs found in the training session,
#
#    The adjusted formula is a little more complicated:
#
#                              counts{nucleotide}[position] + B * bg_pct[nucleotide]
#    q{nucleotide}[position] = ------------------------------------------------
#                                                   N + B
#    where counts and N are as before, and:
#       q{nucleotide}[position] is the adjusted frequency of the nucleotide
#          at the position
#       B is the weighting constant for pseudocounts
#       bg_pct[nucleotide] is the fraction that the nucleotide appears in total
#          DNA
#
#    We actually store $log_q{$nucleotide}[position], the log of the adjusted
#    frequency, rather than the frequency itself.

  sub Calc_nucleotide_frequencies {

     my $motif;
     my $position;
     my $nucleotide;
     my $count;
     my $N = Size(@motif_sequences);

     foreach $nucleotide ("A", "C", "G", "T") {
        foreach $position (@informational) {
           $count = ($counts{$nucleotide}[$position] or 0);
           $log_q{$nucleotide}[$position] =
              log(($count + $B * $bg_frequency{$nucleotide}) / ($N + $B));
        }
     }
  }

#### CALC_BACKGROUND_FREQUENCIES
#    Calculate background percentage (for pseudocounts)
  sub Calc_background_frequencies {
     foreach my $nucleotide ("A", "C", "G", "T") {
        $log_p{$nucleotide} = log($bg_frequency{$nucleotide});
     }
  }

           #### SCORING SUBROUTINES ####

#### ANALYZE_SEQUENCE
#    Runs motif-sized window from beginning to end of
#       sequence
#    Scores each window
#    If score exceeds threshold score, then store info
  sub Analyze_sequence {
     my ($name, $sequence) = @_;
     my $region_start;
     my $subregion;

     foreach $region_start (0 .. length($sequence) - $motif_length) {
        if ($region_start % 20000 == 0) {print $region_start, " "}
        $subregion = substr($sequence, $region_start, $motif_length);
        my $score = score_of_region($subregion);
        if ($score > $threshold) {
           store_hit_information($name, $score, $region_start, $subregion)
        }
     }
  }

#### SCORE_OF_REGION
#    The score is Sum [log_q - log_p] at each position
#    This means that the closer the sequence is to the set that
#       produced q, the higher the score
#    If any region has a gap (".") give it a lethal score
  sub score_of_region {
    my ($region) = @_;
    my $score = 0;
    foreach my $position (@informational) {
       my $nucleotide = substr($region, $position, 1);
       if ($nucleotide eq ".") {return $threshold - 1}
       $score += $log_q{$nucleotide}[$position] - $log_p{$nucleotide};
    }
    return $score
  }

#### STORE_HIT_INFORMATION (name, score, start, subregion)
#    Store the hit information:
#        Name of contig
#        Score of window according to PSSM
#        Start coordinate of window within contig
#        Sequence of window
#    If we've accumulated too many hits, then discard all but
#        the highest-scoring ones.
#    Also reset the threshold to the worst score saved

  sub store_hit_information {
     my ($name, $score, $start, $subregion) = @_;

     push @hit_info, [$score, $name, $start, $subregion];
     if (Size(@hit_info) > 9 * $number_to_save) {
        sort_and_discard_excess_values();
        $threshold = $hit_info[-1][0];
     }
  }

#### SORT_AND_DISCARD_EXCESS_VALUES
#    Sort saved scores highest to lowest
#    Reset last index of hit_info so that only number_to_save are retained
  sub sort_and_discard_excess_values {
     @hit_info = reverse sort { $$a[0] <=> $$b[0] } @hit_info;
     $#hit_info = $number_to_save - 1 if @hit_info > $number_to_save;
  }

              #### MISCELLANEOUS SUBROUTINES ####

##### SIZE (array)
#     Returns size of array
#     Synonym for scalar (array)
  sub Size {return scalar(@_)}
