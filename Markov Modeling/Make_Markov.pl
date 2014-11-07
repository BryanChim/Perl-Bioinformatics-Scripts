#!/usr/bin/perl -w
# Hamlet
#    VERSION: Version 9 (15 June 2004)
#    PURPOSE:  Makes Markov model of Hamlet's speeches
#         Outputs new speeches based on model
#    INPUT FILES: Text file containing training text, in format:
#              {text}
#              {text}
#         Each bracketful of text represents one utterance
#         In creating new speeches, the program will begin with {
#            and end with } (but won't print either symbol)
#    OUTPUT FILES: None
#    DIFFERS FROM VERSION 9:
#         Clarifies program somewhat by adding new variables

############## LIBRARIES AND PRAGMAS ################

  use strict;
  use Storable;
  use FastA_module;

#################### CONSTANTS ######################

  my $true = 1;			# Perl doesn't have logical variables so
  my $false = ! $true;          #    we simulate them
  my $LF = "\n";                # Linefeed
  my $tab = "\t";               # Tab
  my $linewidth = 60;		# Print no more than this many characters/line
  my $begin_speech = "ATG";	# Symbol representing beginning of speech
  my $end_speech = "}";		# Symbol representing end of speech
  my $order = 3;		# Order of Markov Model used
  my $sum_symbol = chr(255);

#################### VARIABLES ######################

  $| = 1;
  my %model;                    # Markov Model for all text
  my @model_keys;		# List of characters encountered in text
	  			#   (i.e., letters not people!)

  my $speech;			# Contains speech currently considered
  my $line;                     # One line of text read from input file
  my @letters;                  # Individual lettes within line of text
  my $starting_letter;          # Beginning letter of moving window 
                                #         through line of text
  my $size_of_line;             # Number of letters in current line of text
  my $last_starting_letter;     # Last position in line of starting letter
  my @window;                   # Current characters, length = order + 1
  my %letters_seen;             # Hash of characters encountered in text 
  my $do_more = $true;		# Variable to tell program when you're through

###################### FILES ########################


my $input_file_name = '6803PHX.nt';
  #my $input_file_name = 'HamletSpeech.txt';
  open INPUT_FILE, "<$input_file_name" or die "Can't open $input_file_name: $!\n";

  my $model_file_name = 'model.dat';




my $seqcount = Read_FastA_sequences ("6803PHX.nt");

#print "SeqCount = $seqcount\n";

my ($header, $sequence);
my $i = 1;
my @seqarray;
while ($i < $seqcount) {
($header, $sequence) = Get_sequence_info ($i);
 print $sequence, $LF;
 push @seqarray, $sequence;
 $i++;
}
################### MAIN PROGRAM ####################
### Three parts:
###   Analyze text: Extracts letters from lines of text 
###                 and adds information to growing Markov model
###   Produce text: Uses Markov model to create new text in the
###                 style of the training set
###   Store model:  Puts Markov model as a hash into a file
   
  ### ANALYZE TEXT
  ###   Read a line of text
  ###   Extract letters from the text
  ###   Pass a window from one end of the line to the other
  ###   Update Markov model according to letters within window
  
  print "Please wait a few seconds for Hamlet to be created...", "\n";
  foreach my $seq (@seqarray) {
  $line = $seq;
  print "sequence\n";
  print "$seq\n";
  # Analyze text
     @letters = Break_up($line);   # Make sure text conforms to proper format
	  			   # Also split line into individual letters
     $size_of_line = Size(@letters);
     $last_starting_letter = $size_of_line - 1 - $order;

                                   # Move window of size order+1 over line
                                   #   incrementing every window encountered
     foreach $starting_letter (0 .. $last_starting_letter) {
        @window = @letters[$starting_letter .. $starting_letter+$order];
        Update_entry(@window);
     }
     }


  print $LF, "UPDATE COMPLETE", $LF;

  Make_keys();                   # Makes list of all characters encountered in text

  ### PRODUCE NEW TEXT
  ###    Use Markov model to make new text
  ###    Print the text
  ###    Ask user if we should do it again

 # while ($do_more == $true) {        # Produce text
 #    $speech = Make_speech ();       # Make rest of text, given first word
 #                                    #     (if specified)
 #    Print_speech ($speech);
#     $do_more = Ask_if_more ();	     # Ask user if another speech wanted
#  }

  ### STORE MODEL

  Output_model ($model_file_name);
  
#################### SUBROUTINES ####################

     ########### SUBROUTINES FOR ANALYZING TEXT ##########

##### BREAK_UP (line of text)
#     Processes and checks line of text
#     Splits text up into individual letters 
#     Returns letters, preceded by null strings to mark beginning
#        The number of null strings added = order of Markov chain
#        Doing this avoids the need to write special code to 
#        handle the beginning of lines
  sub Break_up {
     my ($line) = @_;
     chomp $line;
     $line =~ s/^\s+//; 		# Trim leading blanks
     $line =~ s/\s+$//;		        # Trim trailing blanks
     @letters = split("", $line);       # Separates letters within line
  
     # Die unless the speech is in the correct format of { [text] }
    
#     $letters[0..3] eq $begin_speech            # First letter must be "{"
#        or die "Line $line doesn't start with $begin_speech: $_\n";
#     if ($letters[@letters-1..3] eq "TGA || TAA || TAG)     # Last letter must be "}"
#        or die "Line $line doesn't end with $end_speech: $_\n";
  
#     shift @letters;      # Remove the first character, which is "{";
     return (("") x $order, @letters);
  }
  
##### UPDATE_ENTRY (current window)
#     --> Currently set up for 3rd order Markov <--
#     Splits window up into separate letters
#     Uses letters as key to increment count of window
#     Also increments total count of window - 1 letter
#     Flags window as having been seen

#     Note: The statement
#           $model{$key}++;
#       adds 1 to $model{$key}, if it is defined, and sets it to 1 if
#       it is undefined.

  sub Update_entry {
     my $key = join ("", @_);
     my $current_letter = substr ($key, -1, 1);
     my $old_letters = substr ($key, 0, length($key)-1);
     my $key_sum = $old_letters . $sum_symbol;
     if (not defined ($old_letters)) {die};
     $model{$key}++;
     $model{$key_sum}++;
     $letters_seen{$current_letter} = 1;
  }

##### MAKE_KEYS
#     Makes a list of all characters ever seen in analyzing text  
#     Puts list in @model_keys in arbitrary order
  sub Make_keys {
     @model_keys = keys(%letters_seen);
  }
  
     ########### SUBROUTINES FOR PRODUCING TEXT ###########

##### MAKE_SPEECH  
#     Takes as optional input characters to begin speech
#     Uses Markov hash to add on to initial characters
#     Strategy to initialize window: 
#          Define window initially as empty characters, 
#          Push onto the window the input characters
#          Trim the window until it is as big as the Markov order
  sub Make_speech {
     my ($first_word) = @_;
     my @window = (("") x $order);
     if (defined $first_word) {push @window, split("",$first_word)}
     while (Size(@window) > $order) {shift @window}
     my $speech = $first_word;
     my $next_letter = Find_next_letter (@window);
     while (not $next_letter eq $end_speech) {
        $speech = $speech . $next_letter; # The dot connects first variable to next
  
           # If @window is ("A", "B", "C") and $next_letter is "Z", change
           # @window to ("B", "C", "Z")
           # shift moves every letter over one position to the left and 
           #       then trashes the leftmost letter, $letter[0]
  
        push @window, $next_letter;
        shift @window;			
        $next_letter = Find_next_letter (@window);
     }
     return $speech;
  }

##### FIND_NEXT_LETTER (window)
#     Returns next letter to add on to speech
#     Strategy to find next letter:
#          Rememeber number of times the window has been seen
#          Choose a random integer n between 1 and that number
#          Go through instances window has been seen, in arbitrary order,
#             until nth instance is encountered. 
#          Return letter of that nth instance
  sub Find_next_letter {
     my ($key_base) = join ("", @_);
     my $key_sum = $key_base . $sum_symbol;
     my $sum;
     my $target;
     my $running_sum = 0;
     my $occurrences;
     $sum = $model{$key_sum};
     if (defined($sum)) {    # $sum should ALWAYS be defined, but you never know
        $target = Random_integer (1,$sum);
        $running_sum = 0;
        foreach my $letter (@model_keys) {
            my $key = $key_base . $letter;
            my $occurences = $model{$key};
            if (defined($occurences)) {
               $running_sum = $running_sum + $occurences;
               if ($running_sum >= $target) {return $letter}
            }
        }
     } else {                   # If $sum not defined something horrible happened
        Error_in_find("$key_sum");
        return $end_speech;
     }
     return $end_speech;
  }
  
##### PRINT_SPEECH (speech to print)
#     Strategy: 
#         * Print one line at a time, delete printed text once printed
#         * Make sure line breaks at a space, not in middle of word
#         * To do this, consider n characters of text and work your way
#           inward, one character at a time, until you encounter a space.
#         * n = $linewidth (see CONSTANTS section)
  sub Print_speech {
     my ($speech) = @_;
     print $LF;
     my $right_character = $linewidth - 1;
     while (length($speech) > $linewidth - 1) {
        if (substr($speech, $right_character, 1) eq " ") {
           print substr($speech, 0, $right_character), "\n";
           substr($speech, 0, $right_character+1) = "";
           $right_character = $linewidth - 1;
        } else {$right_character = $right_character - 1}
     }
     print $speech, "\n";
  }
   
     ########### MISCELLANEOUS SUBROUTINES ##########

##### ASK_IF_MORE
  sub Ask_if_more {
     print $LF;
     print "Want another? (Yes/No)  ";
     print "";
     my $key = substr(<STDIN>,0,1);	# Takes just first character of response
     if (uc($key) eq "Y") {return $true}
     else {return $false}
  }
  
##### ERROR_IN_FIND
#     This subroutine should never be called
  sub Error_in_find {
     my ($word) = $_;
     print "\nMISTAKE IN PROGRAM!  Please note following:\n";
     print "     No values for ($word)\n";
  }
  
  
##### RANDOM_INTEGER (low range, high range)
#     Returns a random integer somewhere between low range and high range, inclusive
  sub Random_integer {
    my ($low, $high) = @_;
    my $range = $high - $low + 1;
    my $random_fraction_of_range = rand($range);
    return $low + int($random_fraction_of_range);
  }

##### SIZE (array)
#     Returns size of array (synonymous with scalar)
  sub Size {return scalar(@_)}

##### OUTPUT_MODEL (name of file)
#     Stores %model in given file
#     Uses store method of module Storable
  sub Output_model { 
     my ($model_file_name) = @_; 
     print "Storing model"; 
     store \%model, $model_file_name; 
  } 
