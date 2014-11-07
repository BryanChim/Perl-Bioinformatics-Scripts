#!/usr/bin/perl -w
# UseMarkov
#    VERSION: Version 2 (8 December 2004)
#    PURPOSE:  Applies Markov model to sets of gene sequences
#         Outputs scores of fitness to model
#    INPUT FILES: 
#         Model file:
#              Hash made by MakeMarkov, stored by store method of Storable
#              Sums are stored in hash with key ending with chr(255)
#              Program will accept model of any order
#         Sequence file:
#              FastA formated DNA sequence
#    OUTPUT FILE: Tab-delimited file, with following quantitites:
#         score: Log of probability of gene sequence occuring according to model, 
#                   vs according to background nucleotide frequencies
#                Score is increasingly positive for sequences increasingly close 
#                   to sequences of training set
#         orf_name: First field of FastA header
#         orf_length: Number of nucleotides in FastA sequence
#         orth_eval: E-value of closest cyanobacterial ortholog 
#                   or "none" if no ortholog
#         orf_descr: Remainder of FastA header
#
#    ### INCOMPLETE! ###
#       See Calc_Markov_score

############## LIBRARIES AND PRAGMAS ################

  use strict;
  #use Warnings;
  use Storable;
  use FastA_module;

#################### CONSTANTS ######################

  my $true = 1;	  	        # Perl doesn't have logical variables so
  my $false = ! $true;          #    we simulate them
  my $LF = "\n";                # Linefeed
  my $tab = "\t";               # Tab
  my $sum_symbol = chr(255);    # Key representing sum of elements of hash
  my $fill_char = chr(1);       # Used to pad keys to full length
  my $window_size = 200;        # Number of nucleotides used to calculate score
  my $window_increment = 25;    # Number of nucleotides by which to shift window

#################### VARIABLES ######################

  my %model;                    # Markov Model for all text
  my %aa_frequencies;          # Frequency in sequence of each
  my %model_bg;
  my $order;		        # Order of Markov Model used
  my $N;                        # Number of sequences in training set
  my $number_of_keys;           # Number of letter combinations in model
  my $number_of_windows;  	# Number of windows considered in making the model
  my @protein_info;                 # (Extracted from header) Name, length, Eval of ortholog, annotation

  my $header;                   # Header line in current FastA-formatted orf
  my $sequence;                 # Given sequence 

  my %bg_frequency;             # nucleotide frequencies of overall sequences
  my %log_q;                    # $log_q{$nucleotide}[$position] is the log of the
                                #    frequency of the nucleotide at that position
  my %log_p;                    # $log_p{$nucleotide} is the background frequency
                                #    of the nucleotide
  my $score;                    # Score representing fit of gene with model 
                                #    based on training set
  my @scores;                   # Set of scores derived from window moving 
                                #    through sequence

###################### FILES ########################

  my $model_file_name = 'model.dat';
  %model = Get_Model($model_file_name);        # Restore Markov model

  my $input_file_name = 'BryanChim_6803_prots.fa';
  my $number_of_sequences = Read_FastA_sequences($input_file_name);

  my $output_file_name = 'BryanChim_6803_prot-scores.txt';
  open OUTPUT, ">$output_file_name" or die "Can't open $output_file_name: $!\n";

################### MAIN PROGRAM ####################
####  Calculate amino acid frequencies (may or may not be used)
####  Calculate background frequencies (log_p) based on ...
####  Calculate transition frequencies (log_q) based on model 
####  Go through given sequence, bit by bit, calculating score for each bit,
####     and save the scores

   %aa_frequencies = Calc_aa_frequencies();
   %log_p = Calc_background_frequency();
   %log_q = Calc_motif_frequency(%model);
 #  foreach my $keyss (keys %log_q) {
#   print $log_q{$keyss}, $LF;
 #  }
   print $LF, "Calculating scores...";
   foreach my $seq (0 .. ($number_of_sequences -1))
      {($header, $sequence) = Get_sequence_info($seq);

       @protein_info = Parse_protein_annotation($header);
       $score = Calc_Markov_score($sequence, \%log_p, \%log_q);  # Calculate fit with Markov model
       Output_protein_info();
      }

   print "   DONE!", $LF;

#################### SUBROUTINES ####################

##### GET_MODEL (file name)
#     Retrieves model stored in given file as hash
#     Uses retrieve method of module Storable
#     Calls for calculation of statistics 
#        (only one, the order of the model, is used later)

  sub Get_Model {
     my ($file_name) = @_;
     my $model_reference = retrieve $file_name;
     my %model = %$model_reference;
     print $LF, "MODEL from $file_name:", $LF;
     $order = Stats_of_model(%model);
     return %model;
  }

##### STATS_OF_MODEL (hash)
#     Determines order of retrieved model and number of sequences in training set
#     Strategy: 

  sub Stats_of_model {
     my (%model) = @_;
     my $biggest_length = 0;
     my $current_order;
     my $number_of_instances = 0;
     my $number_of_keys = scalar (keys %model);
     my $number_of_sequences = 0;

     foreach my $key (keys %model) {
        if (substr ($key, -1, 1) eq $sum_symbol) 
           {$number_of_instances = $number_of_instances + $model{$key}}
        elsif (length ($key) == 1) 
           {$number_of_sequences = $number_of_sequences + $model{$key}}
        elsif (length ($key) > $biggest_length) 
           {$biggest_length = length ($key)}
     }
     $current_order = $biggest_length - 1;

     print "   order $current_order", $LF;
     print "   $number_of_sequences sequences in training set", $LF;
     print "   $number_of_keys keys and $number_of_instances instances";
     print $LF, $LF;

     return $current_order;
  }

                ########## SCORING SUBROUTINES ##########

##### CALC_AA_FREQUENCIES
#     Counts nucleotides in given sequence
#     Calculates and returns nucleotide frequencies

  sub Calc_aa_frequencies {
     my ($seq) = @_;         # Given sequence
     my %aa_sum;            # Temporary hash to count amino acids
     my %frequency;          # Frequency of each amino acids
     my $count;              # Temporary count of specific nucleotide 
                             #     in specific sequence
     my $amino_acid;         # Current amino acid
     my $orf;                # Current number of orf


     foreach $orf (1 .. $number_of_sequences) {
        ($header, $sequence) = Get_sequence_info($orf);
        foreach $amino_acid ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") {
           $count = ($sequence =~ s/$amino_acid/$amino_acid/g);
           $aa_sum{$amino_acid} += $count;
           $aa_sum{$sum_symbol} += $count;
        }
     }

     print "AMINO ACID FREQUENCIES:", $LF;
     foreach $amino_acid ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") {
        $frequency{$amino_acid} = $aa_sum{$amino_acid}/$aa_sum{$sum_symbol};
        print "   $amino_acid: ";
        printf "%5.3f$LF", $frequency{$amino_acid};
     }
     return %frequency;
  } 

##### CALC_BACKGROUND_FREQUENCY
#     Counts amino acids in all sequences
#     Calculates amino acid frequencies
#     Calculates log of frequencies, used in later calculations

  sub Calc_background_frequency {
     my %log_p;

     foreach my $amino_acid ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") {
        $log_p{$amino_acid} = log($aa_frequencies{$amino_acid});
        }

     return %log_p;
  } 

##### CALC_MOTIF_FREQUENCY
#     Calculates adjusted frequency (q) of nucleotide at last position 
#        of Markov chain, expressed as log(q)
    
  sub Calc_motif_frequency {
     my (%model) = @_;
     my %log_adj_freq;
     my $given_letters;    # First letters in letter combination
     my $predicted_letter; # Last letter in letter combination
     my $letters_key;      # Full letter combination
     my $instances;        # Number of times letters occur
     my $sum;              # Number of times given letters appear, followed by
                           #   any letter
     my $sum_key;          # Key accessing sum 

     foreach my $key (keys %model) {
        if (not (substr($key, -1) eq $sum_symbol)) {
           ($instances, $sum) = Instances_of_key($key);
           $log_adj_freq{$key} =
                log( $instances/$sum )
        }
     }
     return %log_adj_freq;
  }

##### INSTANCES_OF_KEY
#     Retrieves from model instances of a given key 
#        and the sum of all instances of the key minus the last letter
#     Example:  Instances_of_key("GAGC") would return
#        (number_of_instances_of_GAGC, number_of_instances_of_GAG_anything)  
  sub Instances_of_key {
     my ($key) = @_;                 # all letters in key
     my $given_letters = substr ($key, 0, length($key)-1);
                                     # First letters in key
     my $predicted_letter = substr($key, -1, 1);
                                     # Last letter in key
     my $instances = $model{$key};   # Number of times key occur
     my $sum_key = $given_letters . $sum_symbol;
     my $sum = $model{$sum_key};     # Number of times given letters appear, followed by
                                     #   any letter

    return ($instances, $sum)
  }
  
##### CALC_MARKOV_SCORE (sequence)
#     Calculates score by comparing probability of arriving at given sequence
#        by model frequencies vs background amino acid frequencies
#     Runs a window through the sequence, where the size of the window is the
#        order of the model

  sub Calc_Markov_score {
     my ($seq, $logPref, $logQref) =  (@_);       # Sequence of which to calculate score
     my %logP = %$logPref;
     my %logQ = %$logQref;
  # foreach my $keyss (keys %logP) {
  # print $keyss, $LF;
  # }

     $seq = (" " x $order) . $seq;
                            # Blank letters inserted to set first letter as special
     my $letters;           # Letters (length order + 1) under consideration
     my $given_letters;     # First letters
     my $predicted_letter;  # Last letter
     my $score;             # log_q - log_p for given letters
     my $score_sum = 0;     # Will be sum of scores for each letter group

     foreach my $position (0 .. length($seq) - ($order+1)) {    
                   # Moves order+1th nucleotide through sequence
        $given_letters = substr ($seq, $position, $order);
        $predicted_letter = substr ($seq, $position + $order, 1);
        $letters = substr ($seq, $position, $order+1);
        $letters =~ s/ //g;        # Remove beginning blanks

  if (not defined ($logQ{$letters})) {
  my $QxScore = -9001;
  $score = $QxScore - $logP{$predicted_letter};
  }
     if (defined ($logQ{$letters})) {
        $score = $logQ{$letters} - $logP{$predicted_letter};                                                                                    ### SOMETHING DIFFERENT NEEDS TO GOES HERE

        }
        $score_sum += $score;
     }
     return $score_sum;
  }
##### PARSE_PROTEIN_ANNOTATION (header line)
#     Extract various information from header line of FastA-formatted file
#     Returns protein sequence name, description and gene
  sub Parse_protein_annotation {
     my ($annotation) = @_;

     $annotation =~ /\w{2}\|(.*?)\|.*?\s(.*?)\sOS=.*?GN=(.*?)\s.*/;

     my $protein_seq_name = $1;
     my $protein_seq_description = $2;
     my $gene = $3;

     return ($protein_seq_name, $protein_seq_description, $gene);
  }

##### OUTPUT_ORF
#     Output score and orf info to screen and to file
  sub Output_protein_info {
#    printf "%6.1f", $score;
#    print $tab, join($tab, @orf_info), $LF;
     printf OUTPUT "%6.1f", $score;
     print OUTPUT $tab, join($tab, @protein_info), $LF;
   }
   
##### SIZE (array)
#     Returns size of array (synonymous with scalar)
  sub Size {return scalar(@_)}  
