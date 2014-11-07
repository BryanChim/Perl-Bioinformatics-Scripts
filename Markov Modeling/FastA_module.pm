#!/usr/bin/perl -w
# FASTA_MODULE
#    VERSION: Version 2 (1 Dec 2004)
#    PURPOSE:  Reads single or multiple FastA files
#
#        Read_FastA_sequence (path)
#           Opens file
#           Reads one or more contiguous sequences in FastA format, along with header
#           Saves header and (upper case) sequence
#           Used in conjunction with Get_sequence_info
#           Closes file
#           Returns number of sequences read
#          Example:
#             my $number_of_sequences 
#                 = Read_FastA_sequence("/BNFO601/6803orfs.nt");
#            [All header/sequence pairs will be read in.   
#             You can access them using Get_sequence_info]
#
#        Get_sequence_info (number of sequence)
#           Returns the header and sequence of the n'th sequence stored
#               by Read_FastA_sequence
#           Returns undefined if sequence doesn't exist, perhaps because
#               Read_FastA_sequence hasn't been successfully invoked
#          Example:
#             my ($header, $sequence) = Get_sequence_info(3);
#            [assigns the third header to $header and third sequence
#             to $sequence, from the file read in using Read_FastA_sequence]
#
#        Read_FastA_sequence (path)
#           Opens file
#           Reads first (or only) header and sequence in FastA format
#           Closes file
#           Returns header and (upper case) sequence
#          Example:
#             my ($header, $sequence) = Read_FastA_sequence("sequence.nt");
#            [Assigns the header and sequence in the file to the
#             two variables]
#
#    DIFFERENCES FROM PREVIOUS VERSION: Just a change in documentation

  package FastA_module;
  require Exporter;
  our @ISA = qw(Exporter);
  our @EXPORT = qw(Read_FastA_sequences Get_sequence_info Read_FastA_sequence);

############## LIBRARIES AND PRAGMAS ################

  use strict;

#################### CONSTANTS ######################

  my $true = 1;

#################### VARIABLES ######################

  my @seq_info;            # Stores all headers and sequences:
                           #  2D array: $seq_info[n,0] = nth header
                           #            $seq_info[n,1] = nth sequence

################# INITIALIZATION ####################

  return $true;

#################### SUBROUTINES ####################

##### READ_FASTA_SEQUENCES (path)
#     Reads file consisting of one or more FastA-formatted sequences
#     Saves an array, @seq_info, containing header and sequence
#     Returns number of sequences read
  sub Read_FastA_sequences {
     my ($file_path) = @_;
     open SEQUENCES, "<$file_path" or die "Can't open $file_path: $!\n";

     my $header;              # Current header
     my $sequence;            # Current sequence

     my $line = <SEQUENCES>;
     while (defined $line) {
        chomp $line;
        if (substr($line,0,1) eq ">") {              # Found header
           if (defined $header) {                    #   If not first sequence
              push @seq_info, [$header, $sequence];  #     Save current sequence
           }
           $header = substr($line,1);                #   Save header of next sequence
           $sequence = "";                           #   Initialize next sequence
        } else {
           $sequence .= uc($line);                   # Append line to sequence
        }
        $line = <SEQUENCES>;
     }
     if (defined $header) {                          # If current sequence exists
        push @seq_info, [$header, $sequence];        #   Save it
     }

     close SEQUENCES;
     return scalar(@seq_info);                       # Return number of sequences read
  }

##### GET_SEQUENCE_INFO (n)
#     Returns header and sequence of nth sequence in @seq_info
  sub Get_sequence_info {
     my ($n) = @_;
     my $header = ${$seq_info[$n-1]}[0];
     my $sequence = ${$seq_info[$n-1]}[1];
     return $header, $sequence;
  }

#### READ_FASTA_SEQUENCE (path to file)
#    Opens sequence file in FastA format
#    Returns header and upper cased sequence
  sub Read_FastA_sequence {

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