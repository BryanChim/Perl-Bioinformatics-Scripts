
#!/usr/bin/perl -w

#### Parser_BLASTn_allhits_nohits

#    ORIGINAL AUTHORS: Jeff Elhai and Paul Fawcett (Fall 2004)
#
#	 MODIFICATIONS BY: Bryan Chim (Fall 2013) 
#		** Update regexes to account for change in BLASTn output format
#		** Additional output file for queries that made no hits against the target
#		** Account for query and contig description "spillover" across lines
#
#    PURPOSE: Extract information from BLAST output
#	 		  Compatible with output from BLASTn 2.2.28+
#
#    INPUT: From file $blast_path, with format:
#       (header)
#       ...
#       Query= [query_name] [query_description]
#                 ([query_length] letters)
#       ...
#       >ContigMMM.revised.geneNNN.protein [subject1_description]
#                Length = [subject1_length]
#
#        Score = [garbage], Expect = [expectation1]
#       ...
#       >ContigMMM.revised.geneNNN.protein [subject2_description]
#                Length = [subject2_length]
#
#        Score = [discard], Expect = [expectation2]
#       ...
#
#       NOTE: The Query= line and >Contig line may spill over if the description
#             is too long
#       NOTE: The >Contig line gives the Subject name: MMM.NNN
#
#    OUTPUT: From file $result_path, with format:
#       query_name query_description query_length subject1_name
#              subject1_description subject1_length
#
#****OUTPUT 2: From file $result_path_nohits:
#		- Same format as $result_path output -- restricted to information for
#			queries that yielded no hits against target
#
#       NOTE: If there is more than one matching subject, then
#             subject2_name subject2_description subject2_length is appended
#             to the end of the line.
#	      IF there is no matching subject, then the line ends with query_length
#
#       NOTE: There is a line for every query submitted to Blast
#
#       NOTE: Each item is separated from the next by a TAB

############## LIBRARIES AND PRAGMAS ################
  use strict;
  use warnings;

#################### VARIABLES ######################

  my $query_name;             # Name of query (e.g. all0001)
  my $query_description;      # Description of query (e.g. unknown protein)
  my $query_length;           # Length of query (e.g. 173 [bp])
  my $subject_name;           # Name of subject (e.g. 656.037)
  my $subject_description;    # Description of subject (database doesn't have any)
  my $subject_length;         # Length of query (in bp)
  my $expectation;            # E-value (e.g. 2e-017)
  my @query_info;             # Stores all above values
  my $line;                   # One line from input file
  my $more_description;

###################### FILES ########################

#### OPEN THE BLAST FILE

  my $blast_path = "C:/Users/Public/O157K12COMP";
  open BLAST, "<$blast_path" or die "Can't open $blast_path: $!\n";

#### OPEN OUR OUTPUT FILES

  my $result_path = "C:/Users/Public/O157K12COMPparse";
  open RESULT, ">$result_path" or die "Can't open $result_path: $!\n";
  
  my $result_path_nohits = "C:/Users/Public/O157K12COMPnohits";
  open RESULTS_NOHITS, ">$result_path_nohits" or die "Can't open $result_path_nohits: $!\n";

################### MAIN PROGRAM ####################
#
# We'll accumulate the items of information for each query in the
# array @query_info
#
# Then, for each line of the BLAST file, we check for two things: the start
# of a query, and the start of a match.
#
# When we find the start of a query, we print the accumulated information for
# the previous query. (But when we see the first query there won't be a previous
# query.  What then?  See below at PRINTING THE PREVIOUS QUERY.)
#
# The other possibility is that we might see the start of a match.  If so,
# we record it in @query_info.

$line = <BLAST>;                      
while (defined $line) 
{

	if ($line =~ /^[*****]/)	# In BLASTn 2.2.28, "No-Hit" queries begin with "***** NO HITS ..."***
	{
		print_previous_query_nohits();  
	}

                    
	if ($line =~ /^Query=/) 	# If 'Query' matches beginning of line, print previous query data
	{          
        print_previous_query();
        start_new_query();
    }

    if ($line =~ /^>Contig/) 	# If '>Contig' matches beginning of line, get subject data
	{ 
	record_subject() 	 
	}
                                       
    $line = <BLAST>;                   # Read next line
}

# Now for the last query.  All the others were printed when we were
# about to process the query after them; but there is no query after the
# last query, so we make an extra call to print_previous_query().

  print_previous_query();

#################### SUBROUTINES ####################

#### PRINT_PREVIOUS_QUERY and PRINT_PREVIOUS_QUERY_NOHITS
# Prints info stored in @query_info
# If there's no info there (i.e., the first time through), print nothing

sub print_previous_query 
{
   if (@query_info) 
   {
      my $result_line= join(",", @query_info); # Put comma between data
      print $result_line, "\n";                 # Print to monitor
      print RESULT $result_line, "\n";          # Print to file
#      @query_info = ();                   # Reinitialize storage to nothing
   }
}

sub print_previous_query_nohits 
{
   if (@query_info) 
   {
		my $result_line= join(",", @query_info); # Put comma between data
		print RESULTS_NOHITS $result_line, "\n";          # Print to file
                             # Reinitialize storage to nothing
   }
}

#### START_NEW_QUERY
# Extracts query information, saves in @query_info

sub start_new_query 
{

   ($query_name, $query_description) = ($line =~ /^Query=\s(\w+)\W+(.+)/);

   $line = <BLAST>;
    chomp $line;
    
	if (defined ($line))
    {
    $query_description = ($query_description . " " . $line);
    }
    
	unless (defined $query_name) 
    {
      warn "Didn't find query name on line $.\n";
      return;
    }

   $line = <BLAST>;
   
	while (defined $line) 
	{

		($query_length) = ($line =~ /Length=(\d+)/);
		
		if (defined $query_length) 
		{
			@query_info = ($query_name, $query_description, $query_length);
			return;
		}
		
      $line = <BLAST>;
	  
	}
	
}

#### RECORD_SUBJECT
# Extracts subject info, stores in @query_info

sub record_subject 
{

   find_subject_name();
   unless (defined $subject_name) 
   {
      warn "Didn't find subject_name on line $.\n";
      return;
   }

   find_subject_length();
   unless (defined $subject_length) 
   {
      warn "Didn't find subject_length on line $.\n";
      return;
   }

   find_expectation();
   unless (defined $expectation) 
   {
      warn "Didn't find expectation of chance match on line $.\n";
      return;
   }

   push @query_info, $subject_name, $subject_description, $subject_length,
        $expectation;

}

#### FIND_SUBJECT_NAME
# Extracts name of subject

sub find_subject_name 
{
   $subject_name = undef;
   if ($line =~ /^>Contig(\d+)\.revised\.gene(\d+)\.protein\W+(.*)/) 
   {
     $subject_name = sprintf "%03d.%03d", $1, $2;
     $subject_description = $3;
   }
}

#### FIND_SUBJECT_LENGTH
# Extracts length of subject

sub find_subject_length 
{
   $subject_length = undef;
   $line = <BLAST>;
   
   while (defined $line) 
   {
      ($subject_length) = ($line =~ /Length= (\d+)/);
	  
      if (defined $subject_length) 
	  {return;}
	  
      $line = <BLAST>;
   }
}

#### FIND_EXPECTATION
# Extracts E-value for comparison

sub find_expectation 
{

   $expectation = undef;
   $line = <BLAST>;
   
   if (defined $line) { $line = <BLAST> }
   
   unless (defined $line) 
   {
      warn "Incomplete query at $.\n";
      return;
   }
   
   ($expectation) = ($line =~ /Expect = (\S+)/);

}


