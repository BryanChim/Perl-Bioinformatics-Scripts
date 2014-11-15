Generates a Markov Model of user-specified order from a user-inputted nucleotide training set.
The resultant model can then be used to score a set of contigs/ORFs/whatever-you-like to check 
for orthology or xenology.

Basic workflow is as follows:

1) Training set is fed into NUCLEOTIDE_Make_Markov.pl ~ FROM $input_file_name
    (example training set used is 6803PHX.nt)
  * use PROTEIN_Make_Markov.pl if you wish to translate the sequences and model protein instead
  * $input_file_name and $model_file_name should be modified appropriately
  
2) ..._Make_Markov program generates Markov Model of <$order> Order 
   which is outputted to a Storable binary object ~ FROM $model_file_name
  * if you wish, use Display_hash.pl on this Storable object to generate a human-readable output
     (see 1st_Order_proteinhash.txt, 3rd_Order_nucleotidehash.txt or 3rd_Order_proteinhash.txt as examples)
     
3) Feed your model and sequences-to-be-scored into NUCLEOTIDE/PROTEIN_Use_Markov.pl 
    ~ FROM $model_file_name and FROM $input_file_name respectively
    (example test files used are 6803Orfs.nt (nucleotide) and 6803prots.fa (protein))
  
4) Resultant scores, calculated as summed log-odds ratios, are outputted to a score file
    alongside the annotations for each of your sequences ~FROM $output_file_name
    

NOTES:

** Generally speaking, higher and more positive scores indicate closer relation to the training set

** If you choose a very high (3rd+) Markov Order or if your training set is not sufficiently large, 
    certain transition states may never be seen and, therefore, not be accounted for in the model. 
    This can especially be a problem when creating a model for protein sequences, since there are 
    so many more possible transition states. 
    
    (20 amino acids vs 4 nucleotides -- so for a 3rd order protein Markov Model there are:
      20 ^ 3 or 8000 possible  states --> 160000 possible state/transition combinations)
    
*** You would either need to lower your Markov ORDER, increase the size and/or robustness of your 
    training set, or distribute pseudo-counts to adjust for non-occurrences. 
    
*** I've currently and naively assigned a constant negative pseudoscore for those cases where 
    a transition was never counted ($dummy_pseudo_score -- currently set to -8). This will need 
    to be tweaked for greater accuracy/validity sometime the future.
