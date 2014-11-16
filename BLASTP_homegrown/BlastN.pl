#!/usr/bin/perl -w
# BlastN 
#    VERSION: Version 3.2 (10 October 2003)
#    PURPOSE: Basic local alignment sequence tool
#             Local alignment of two DNA sequences
#    INPUT FILES:  
#             1. FastA file used as query
#             2. FastA file used as target (=subject)
#    OUTPUT FILES: Output to screen
#             Gives each extended match found, including coordinates

############## LIBRARIES AND PRAGMAS ################
  use strict;

#################### CONSTANTS ######################

  my $true = 1;                # Perl has no specific value for true
                               #     so we supply one
  my $false = not $true;       # Equivalent to $false = 0
  my $gap = "-";               # Used to indicate gap in alignment
  my $bar = "|";               # Used to indicate match between sequences
  my $LF = "\n";               # Linefeed
  my $backwards = -1;          # Indicates desire to extend match left to right
  my $forwards = +1;           # Indicates desire to extend match right to left
  my $max_segment_length = 50; # Nucleotides per line in output

#### BLAST PARAMETERS

  my $match_reward = 1;
  my $mismatch_penalty = 2;
  my $gap_begin_penalty = 5;
  my $gap_extend_penalty = 2;
  my $word_length = 5;

#################### VARIABLES ######################

  my $query_name;             # Header of FastA file used as query for Blast
  my $query;                  # Sequence used as query
  my $query_word_start;       # Position within query where word match found
  my $query_word;             # Sequence of word match
  my $query_word_pattern;     # 
  my $effective_query_end;    # Length of query, less word length
                              #    (actually last coordinate - word length + 1)
  my $target_name;            # Header of FastA file used as target (subject) for Blast
  my $target;                 # Sequence used as target
  my $target_word_start;      # Position within target where word match found

#### USED ONLY IN SUBROUTINES
#    These variables are used in multiple subroutines and so are declared
#    here as global variables, available to all subroutines

  my $direction;              # Forward if right of word match, otherwise backwards         
  my $target_boundary;
  my $query_boundary;
  my $max_score;              # High score in table
  my $max_i;                  # Position of the high score in table (i is y-axis)
  my $max_j;                  # Position of the high score in table (j is x-axis)
  my @score;                  # Scoring table, used for both directions of extensions
  my @traceback_i;            # Table keeps track of source of scores in scoring table (i)
  my @traceback_j;            # Table keeps track of source of scores in scoring table (j)
  my %seen;                   # Keeps track of which full matches have already been
                              #   encountered, to avoid duplications

###################### FILES ########################

  my $query_path = "q.txt";
  my $target_path = "t.txt";

################### MAIN PROGRAM ####################

#### INITIALIZE QUERY AND TARGET
#    Read query and target (subject) from files
#    Calculate last coordinate in query that could support a word match

  ($query_name, $query) = Get_FastA_sequence($query_path);
  ($target_name, $target) = Get_FastA_sequence($target_path);
  $effective_query_end = (length($query) - 1) - $word_length + 1; 

  print "Query: $query_name", $LF;
  print "Target: $target_name", $LF;
  print $LF;

#### RUN WINDOW THROUGH QUERY, NUCLEOTIDE BY NUCLEOTIDE
#    Move window from beginning to (effective) end of query
#    Extract word from query. Compiling the pattern ($query_word_pattern) 
#       saves execution time.
#    For each word match, find position of target
#       and try to extend match in both directions (Process_match)
#    Begin next search for a match just beyond current match, to allow
#       for the possibility of overlapping matches

  for $query_word_start (0 .. $effective_query_end) {
     $query_word = substr($query, $query_word_start, $word_length);
     $query_word_pattern = qr/$query_word/;
     while ($target =~ /$query_word_pattern/g) {
        $target_word_start = pos($target) - $word_length;
        Process_match($target_word_start, $query_word_start);
        pos($target) = $target_word_start+1;
     }
  }

#################### SUBROUTINES ####################

########## Initialization subroutines ##########

#### GET_FASTA_SEQUENCE (path to file)
#    Opens sequence file in FastA format
#    Returns header and upper cased sequence
  sub Get_FastA_sequence {

   my ($path) = @_;
   my $line;
   my $header;
   my $sequence = "";
   
   open SEQUENCE_FILE, "<$path" or die "Can't open $path: $!\n";
   $line = <SEQUENCE_FILE>;        # Read first line (must be header)
   chomp $line;                    # Remove line break
   $header = substr($line, 1);     # Saves header, discards ">"

   $line = <SEQUENCE_FILE>;        # Read next line
   while (defined $line) {         # While file still has lines...
      chomp $line;                 # Remove line break
      $line =~ s/\r//;             # And carriage return (if any)
      $sequence .= $line;          # Append the line to the sequence
      $line = <SEQUENCE_FILE>;     # Read next line
   }
 
   close SEQUENCE_FILE;
   return ($header, uc($sequence));
}

########## SCORING AND MATCHING SUBROUTINES ##########

#### PROCESS_MATCH (starting coordinate of target word, starting coordinate of query word)
# Calculates query and target words
# Extends matches in forwards then backwards directions
# Assembles full target hit and query hit
# Prints hit
  sub Process_match {
     ###  Calculate query and target words
     my ($target_word_start, $query_word_start) = @_;
     my $target_word_end = $target_word_start + $word_length - 1;
     my $query_word_end = $query_word_start + $word_length - 1;
     my $target_word = substr($target, $target_word_start, $word_length);
     my $query_word = substr($query, $query_word_start, $word_length);
   
     ###  Extend match in each direction
     my ($target_max, $query_max, $forward_target, $forward_query)
        = Extend_match_in_one_direction($forwards, $target_word_end, $query_word_end);

     my ($target_min, $query_min, $backward_target, $backward_query)
        = Extend_match_in_one_direction($backwards, $target_word_start, $query_word_start);

     ###  Assemble full target hit and query hit
     my $target_hit = $backward_target . $target_word . $forward_target;
     my $query_hit = $backward_query . $query_word . $forward_query;

     Print_hit($target_min, $target_hit, $query_min, $query_hit);
  }

#### EXTEND_MATCH_IN_ONE_DIRECTION 
#          (direction to extend, boundary of target word, boundary of query word)
#    Initialize scoring table and traceback table
#    Calculate scoring table (Go_through_diagonal), starting with first element (0,0)
#    Construct and return extension sequences
#    Return reverse of sequence of extension was going backwards
#    Returns:
#      last index in target of target extension sequence
#      last index in query of query extension sequence
#      target extension sequence
#      query extension sequence
#
  sub Extend_match_in_one_direction {
     ($direction, $target_boundary, $query_boundary) = @_;
     @traceback_i = ();          # Initialize traceback table...
     @traceback_j = ();          #
     $traceback_i[0][0] = -1;    # ... setting first element (corner) to -1
     $traceback_j[0][0] = -1;    #     to indicate that the score has no source
     $max_i = 0;                 # Initialize position of best score to 
     $max_j = 0;                 #     first element
     @score = ();                # Initialize scoring table
     $score[0][0] = $match_reward * $word_length;
     $max_score = $score[0][0];

     Go_through_diagonal(0,0);

     my ($target_extension, $query_extension) = Make_extension_sequences();
     if ($direction eq $forwards) {
        return ($target_boundary + $max_i, $query_boundary + $max_j, 
                $target_extension, $query_extension);
     } else {
        return ($target_boundary - $max_i, $query_boundary - $max_j, 
                reverse($target_extension)."", reverse($query_extension)."");
     }
  }

#### MAKE_EXTENSION_SEQUENCES
#    Construct $target_extension and $query_extension by tracing back from
#      the high score (at position $max_i,$max_j) found in Go_through_diagonal
#    Strategy is to start at the high score and go back to the beginning of
#      the scoring table.
#    If a scoring element was reached by a match or mismatch, then put nucleotide 
#      in the extensions.
#    If a scoring element was reached by a insertion/deletion, then put a gap
#      symbol in the appropriate extension and a nucleotide in the other.
#    Extensions grow backwards (from high score to beginning), so the string
#      has to be reversed at the end of the process
  sub Make_extension_sequences {
     my $target_extension = "";
     my $query_extension = "";
     my $i = $max_i;
     my $j = $max_j;
     my $i0;
     my $j0;

     while ($i > 0 or $j > 0) {        # Loop until both dimensions back to beginning
        $i0 = $traceback_i[$i][$j];    # ($i0,$ij) is the parent of ($i,$j)
        $j0 = $traceback_j[$i][$j];    # 

        while ($i > $i0 or $j > $j0) {       # Loop once unless a dimension has a gap,
                                             # indicated if index = parent index
           if ($i == $i0) { $target_extension .= $gap}
           else {
              $target_extension .= substr($target, $target_boundary + $i*$direction, 1);
              $i = $i-1;
           }
           if ($j == $j0) { $query_extension .= $gap }
           else {
              $query_extension .= substr($query, $query_boundary + $j*$direction, 1);
              $j = $j-1;
           }
        }  
     }
     return (reverse($target_extension)."", reverse($query_extension)."");
  }

#### GO_THROUGH_DIAGONAL
# Fills in scoring table for one side flanking the word match
# Finds extension by considering each query-target nucleotide pair
#    along the diagonal
# Branches out vertically and horizontally out from diagonal
  sub Go_through_diagonal {
     my ($starting_t, $starting_q) = @_;
     my $current_t = $starting_t;
     my $current_q = $starting_q;
     my ($max_in_row, $max_in_column);
     my $new_score = 0;
     do {
       $max_in_row = Fill_row_of_scoring_table($current_t, $current_q);
       $max_in_column = Fill_column_of_scoring_table($current_t, $current_q);     
       $current_t = $current_t + 1;
       $current_q = $current_q + 1;
       my $vertical_score = Vertical_score_of($current_t,$current_q);
       my $horizontal_score = Horizontal_score_of($current_t,$current_q);
       my $diagonal_score = Diagonal_score_of($current_t,$current_q);
       $new_score = Max_of($vertical_score, $horizontal_score, $diagonal_score, 0);
       if ($new_score == $vertical_score) {
          Record_score_change($current_t-1,$current_q, $current_t,$current_q, $new_score);
       } elsif ($new_score == $horizontal_score) {
          Record_score_change($current_t,$current_q-1, $current_t,$current_q, $new_score);
       } elsif ($new_score == $diagonal_score) {
          Record_score_change($current_t-1,$current_q-1, $current_t,$current_q, $new_score);
       } else {
          Record_score_change(0,0, $current_t,$current_q, 0);
       }
     } until ($new_score <=0 and $max_in_row <=0 and $max_in_column <=0);
  }

#### FILL_COLUMN_OF_SCORING_TABLE (starting i, starting j)
#    Fills column of scoring table from given starting point
#    Does this in two passes
#    In first pass, considers scores from all possible directions
#    In second pass, considers scores only from vertical direction
#    Second pass 
  sub Fill_column_of_scoring_table {
     my ($i, $j) = @_;
     my $vertical_score = $score[$i][$j] || 0;
     my @first_position = ($i);
     my $penalty = $gap_begin_penalty;
     my $last_i = 0;
     my $previous_last_i = 0;
     my $new_score;
     my $old_score;
     my $max_in_column = 0;
  
     do {
        $i = $i + 1;
        $old_score = $score[$i][$j] || 0;
        $vertical_score = $vertical_score - $penalty;
        my $horizontal_score = Horizontal_score_of($i,$j);
        my $diagonal_score = Diagonal_score_of($i,$j);
        $new_score = Max_of($old_score, $vertical_score, 
                               $horizontal_score, $diagonal_score, 0);
        if ($new_score >0) {
           if ($new_score == $vertical_score) {
              Record_score_change($i-1,$j, $i,$j, $new_score);
           } elsif ($new_score == $horizontal_score) {
              Record_score_change($i,$j-1, $i,$j, $new_score);
              push @first_position, $i;
           } elsif ($new_score == $diagonal_score) {
              Record_score_change($i-1,$j-1, $i,$j, $new_score);
              push @first_position, $i;
           }
           $max_in_column = Max($max_in_column, $new_score);
        }
        $penalty = $gap_extend_penalty;
     } until ($new_score == 0 & $i > $previous_last_i);
     while (@first_position) {
        $i = shift(@first_position);
        $penalty = $gap_begin_penalty;
        $vertical_score = $score[$i][$j] || 00;
        do {
           $i = $i + 1;
           $old_score = $score[$i][$j] || 0;
           $vertical_score = $vertical_score - $penalty;
           if ($vertical_score > $old_score) {
              Record_score_change($i,$j-1, $i,$j, $vertical_score);
           }
           $last_i = Max($i, $last_i);
           $penalty = $gap_extend_penalty;
        } until ( ($vertical_score <= $penalty) or ($vertical_score <= $old_score) )
     }
     $previous_last_i = $last_i;
     return $max_in_column;
  }
  
#### FILL_ROW_OF_SCORING_TABLE (starting i,j positions)
#
  sub Fill_row_of_scoring_table {
     my ($i, $j) = @_;
     my $horizontal_score = $score[$i][$j] || 00;
     my @first_position = ($j);
     my $penalty = $gap_begin_penalty;
     my $last_j = 0;
     my $previous_last_j = 0;
     my $new_score;
     my $old_score;
     my $max_in_row = 0;
  
     do {
        $j = $j + 1;
        $old_score = $score[$i][$j] || 0;
        $horizontal_score = $horizontal_score - $penalty;
        my $vertical_score = Vertical_score_of($i,$j);
        my $diagonal_score = Diagonal_score_of($i,$j);
        $new_score = Max_of($old_score, $vertical_score, 
                               $horizontal_score, $diagonal_score, 0);
        if ($new_score > 0) {
           if ($new_score == $horizontal_score) {
              Record_score_change($i,$j-1, $i,$j, $new_score);
           } elsif ($new_score == $vertical_score) {
              Record_score_change($i-1,$j, $i,$j, $new_score);
              push @first_position, $j;
           } elsif ($new_score == $diagonal_score) {
              Record_score_change($i-1,$j-1, $i,$j, $new_score);
              push @first_position, $j;
           }
           $max_in_row = Max($max_in_row, $new_score);
        }
        $penalty = $gap_extend_penalty;
     } until (($new_score == 0) and ($j > $previous_last_j));
  
     while (@first_position) {
        $j = shift(@first_position);
        $penalty = $gap_begin_penalty;
        $horizontal_score = $score[$i][$j] || 00;
        do { 
           $j = $j + 1;
           $old_score = $score[$i][$j] || 0;
           $horizontal_score = $horizontal_score - $penalty;
          if ($horizontal_score > $old_score) {
              Record_score_change($i,$j-1, $i,$j, $horizontal_score);
           }
           $last_j = Max($j, $last_j);
           $penalty = $gap_extend_penalty;
        } until ( ($horizontal_score <= $penalty or $horizontal_score <= $old_score) )
     }
     $previous_last_j = $last_j;
     return $max_in_row;
  }
  
#### VERTICAL_SCORE_OF (current target, query position in scoring table)
#    Returns score of upper element, less gap begin penalty
#    Returns 0 if upper element has no score
#    N.B. This is not quite accurate: Function should check traceback to
#      see which kind of gap penalty is appropriate
  sub Vertical_score_of {
     my ($current_t, $current_q) = @_;
     if (defined $score[$current_t-1][$current_q]) {
        return $score[$current_t-1][$current_q] - $gap_begin_penalty;
     } else {
        return 0;
     }
  }

#### HORIZONTAL_SCORE_OF (current target, query position in scoring table)
#    Returns score of adjacent element, less gap begin penalty
#    Returns 0 if adjacent element has no score
#    N.B. This is not quite accurate: Function should check traceback to
#      see which kind of gap penalty is appropriate
  sub Horizontal_score_of {
     my ($current_t, $current_q) = @_;
     if (defined $score[$current_t][$current_q-1]) {
        return $score[$current_t][$current_q-1] - $gap_begin_penalty;
     } else {
        return 0;
     }
  }

#### DIAGONAL_SCORE_OF (current i, j positions)
#    Looks diagonally backwards from current element
#    Returns 0 if current element is on a boundary
#    Otherwise returns diagonal element plus reward (if match)
#      or penalty (if mismatch)
  sub Diagonal_score_of {
     my ($i, $j) = @_;
     if ($i==0 | $j==0) {return 0}
     my $target_pos = $target_boundary + $i*$direction;
     if ($target_pos < 0 or $target_pos >= length($target)) {
       return 0;
     }

     my $query_pos = $query_boundary + $j*$direction;
     if ($query_pos < 0 or $query_pos >= length($query)) {
       return 0;
     }   

     my $change_amount;
     if (substr($target, $target_pos, 1) eq substr($query, $query_pos, 1)) {
        $change_amount = $match_reward 
     } else { $change_amount = -$mismatch_penalty }

     if (defined $score[$i-1][$j-1]) {
        return $score[$i-1][$j-1] + $change_amount;
     } else {
        return 0;
     } 
  }

#### RECORD_SCORE_CHANGE (parent's i coord, parent's j coord, i coord, j coord, score)
#    Saves current score in scoring table
#    Also records the direction from which the score came
#      (i.e. the parent's coordinates)
  sub Record_score_change {
     my ($parent_i, $parent_j, $i, $j, $new_score) = @_;

     $score[$i][$j] = $new_score;
     $traceback_i[$i][$j] = $parent_i;
     $traceback_j[$i][$j] = $parent_j;

     if ($new_score > $max_score) {
        $max_score = $new_score;
        $max_i = $i;
        $max_j = $j;
     }
  }


########## PRINTING SUBROUTINES ##########

#### PRINT_HIT (target coordinate, target sequence
#                query coordinate,  query sequence)
#    Converts Perl coordinates to human coordinates (+1)
#    Keeps track of what full matches has been seen (in %seen)
#    Splits full match to be printed into segments determined by 
#       $max_segment_length
  sub Print_hit {
     my($target_min, $target_hit, $query_min, $query_hit) = @_;
     $target_min = $target_min+1; # Perl counts from 0
     $query_min = $query_min+1;   # Perl counts from 0

     my $key = "$target_min:$target_hit:$query_min:$query_hit";
     return if $seen{$key};
     $seen{$key} = $true;

     my $start = 0;
     while ($start < length($target_hit)) {
        my $segment_length
           = Min($max_segment_length, length($target_hit) - $start);
        my $query_segment = substr($query_hit, $start, $segment_length);
        my $target_segment = substr($target_hit, $start, $segment_length);
        my $query_length = Count_bases($query_segment);
        my $target_length = Count_bases($target_segment);

        Print_segment("Query", $query_segment, $query_min, $query_length);
        Print_alignment($target_segment, $query_segment);
        Print_segment("Target", $target_segment, $target_min, $target_length);
        print $LF;

        $target_min = $target_min + $target_length;
        $query_min = $query_min + $query_length;
        $start = $start+$segment_length;
     }
  }

#### COUNT_BASES (sequence with gaps)
#    Counts the number of nucleotides, ignoring gaps
#    Substitutes nothing for gaps then counts what remains
  sub Count_bases {
     my ($segment) = @_;
     $segment =~ s/$gap//go;
     return length($segment);
  }

#### PRINT_SEGMENT (Query: or Target:, sequence to print, 
#                   starting coordinate, length of sequence)
#    Prints one line of output, formatted
  sub Print_segment {
     my ($name, $segment, $start, $length) = @_;
     printf "%-7s %7d %s %7d\n",
        "$name:", $start, $segment, $start+$length-1;
  }

#### PRINT_ALIGNMENT (query string, target string)
#    Compares query and targets and constructs alignment line
#      Does comparison by splitting strings into arrays and comparing
#      element by element
#    When a match, puts bar in alignment line
#    Otherwise, puts a space in alignment line
#    Prints line after its constructed
  sub Print_alignment {
     my ($q_string, $t_string) = @_;
     my @q = split("", $q_string);
     my @t = split("", $t_string);
     my $bars = "";
     my $last_index_of_q = @q - 1;
     my $k;

     for $k (0 .. $last_index_of_q) {
        if ( ($q[$k] eq $t[$k]) and ($q[$k] ne $gap) ) {
           $bars .= $bar;               # Insert | for aligned nucleotide
        } else { 
           $bars .=  " ";               # Insert space for unaligned nucleotide
        }
     }
     print " " x 16, $bars, $LF;        # indent then print alignment line
  }

########## UTILITIES ##########

#### MAX (m and n)
#    Returns the larger of the two numbers
  sub Max {
     my ($a, $b) = @_;
     if ($a > $b) { return $a }
     else { return $b } 
  }

#### MIN (m and n)
#    Returns the lesser of the two numbers
  sub Min {
     my ($a, $b) = @_;
     if ($a < $b) { return $a }
     else { return $b }
  }

#### MAX_OF (list of numbers)
#    Returns the largest of the numbers in the list
  sub Max_of {
     my @numbers = @_;
     my $max_of;
     for my $n (0 .. scalar(@numbers)-1) {
        if (defined $max_of) {
           if ($numbers[$n] > $max_of) {$max_of = $numbers[$n]}
        } else {
           $max_of = $numbers[$n];
        }
     }
     return $max_of; 
  }
   

