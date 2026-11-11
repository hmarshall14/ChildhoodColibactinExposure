Contains output from signature extraction and attribution.

`sigprofiler_sbs_patients`, `msigHDP_patients` provide outputs from signature extraction using these two modalities.

`mmsig_attribution`, `msa`_id_attribution`, `msa_sbs_attribution`, `musical_attribution`, `sigfit_attribution`, `sigprofiler_bootstrapped` contain output files 
from attribution of mutational signatures to each crypt.

`msa_phylogenetic_branches_attribution` contains the output of MSA attribution to branches.

`sigprofiler_sbs_coding_mutations` contains the output of sigprofiler extraction from coding mutations.

`input_matrices` contains all input signature matrices relevant for extraction and or attribution. It contains matrices derived from SigProfiler MatrixGenerator for both SBS and ID. For SBS there is a matrix of unique mutations per donor, a matrix of unique mutations per crypt, and a matrix of unique mutations per phylogenetic branch. For ID there is a matrix of unique mutations per donor. 


