#!/usr/bin/perl -w
# Display_hash
#    VERSION: Version 6 (12 June 2004)
#    PURPOSE: Prints out the keys and values of a stored hash
#    INPUT FILES: Must be stored using the store method of the module Storable
#          The key "|" has a special meaning of summing previous values
#    OUTPUT FILES: Output is to screen or file
#    DIFFERS FROM previous version in:
#          - Allows output to file
#          - Models must be a one-dimensional hash
#          - Sorts keys properly, taking into account those that are
#               not full length

############## LIBRARIES AND PRAGMAS ################

  use strict;
  use Storable;

#################### CONSTANTS ######################

  my $LF = "\n";
  my $true = 1;
  my $false = 0;
  my $sum_symbol = chr(255);    # Identifies hash key as holding sum
  my $null_char = chr(1);       # used to fill hash keys less than full length

#################### VARIABLES ######################

  my %model;           # model whose keys are to be displayed
  my @key_info;        # hash keys in two forms:
                       #   alphabetizable form  (1st element of list)
                       #   actual key           (2nd element of list)

###################### FILES ########################
# Requests name of file using STDERR so as not to contaminate output to
#    STDOUT, in case you want to redirect it to a file

  print STDERR 'Enter file containing model (default = "Model.dat"):  ';
  my $input_file_name = <STDIN>;
  chomp $input_file_name;
  if ($input_file_name eq "") {
      $input_file_name = "Model.dat"}
  Get_Model($input_file_name);

  print STDERR 'Enter file to receive output (default = SCREEN):  ';
  my $OUTPUT;
  my $output_file_name = <STDIN>;
  chomp $output_file_name;
  if ($output_file_name eq "") {
     $OUTPUT = *STDOUT;}
  else {
     open OUTPUT, ">$output_file_name"
           or die "Can't open $output_file_name: $!\n";
     $OUTPUT = *OUTPUT;}
           
################### MAIN PROGRAM ####################

  @key_info = Sort_keys (\%model);
  Print_keys (\%model);

#################### SUBROUTINES ####################

##### GET_MODEL (path to file containing model)
#     Reads model using retrieve method from module Storable
#     Retrieve stores model in a hash and gives the reference to $model_reference
#     Reference allows us to define %model
  sub Get_Model {
     my ($model_file_name) = @_;
     my $model_reference = retrieve $model_file_name;
     %model = %$model_reference; 
  }

##### SORT_KEYS (hash reference)
#     Creates two-dimensional list of keys
#         1st element is the alphabetizable representation of the key
#         2nd element is the real key
#     Sorted by 1st element

  sub Sort_keys {
    my ($hash_ref) = @_;
    my %hash = %$hash_ref;
    my @raw_keys = keys %hash;
    my @key_info;
    my $alpha_key;
    my $max_length = 0;
    
    foreach my $key (@raw_keys) {        # find order
      if (length($key) > $max_length) {$max_length = length($key)}}

    foreach my $real_key (@raw_keys) {
      if (length($real_key) < $max_length) {
         $alpha_key = substr($null_char x $max_length . $real_key,
                                        -$max_length, $max_length)}
      else {$alpha_key = $real_key;}
      push @key_info, [$alpha_key, $real_key];}

    return sort {$a->[0] cmp $b->[0]} @key_info;
  }
  
##### PRINT_KEYS (hash reference)
#     Prints each key and number of instances
#     if key indicates sum, then print underscore and offset sum

  sub Print_keys {
    my ($hash_ref) = @_;
    my %hash = %$hash_ref;
    my $last_key;

    foreach my $key (@key_info) {
       my $alpha_key = $$key[0];
       my $real_key = $$key[1];
       $last_key = substr ($alpha_key, -1, 1);
       if ($last_key eq $sum_symbol) {
          print $OUTPUT "-" x (length($alpha_key) + 13), " "}
       else  {
          print $OUTPUT $alpha_key, "            "}
       print $OUTPUT $hash{$real_key}, $LF;}
  }
