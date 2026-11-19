# Rca3_R3 centromeric analysis

### Description
- the project includes file out_hor.rawTOP5b_rc.fa which contains centromeric repeat sequences in fasta format
- sequences has name like this R1L3::2867492-2867968::+ nHOR-1_2_3::rHOR-2_3_1 where R1...R5 is type of sequence
- R1, R2, R3 are active repeats while R4 and R5 are inactive
- the inactive R4 and R5 contains single nucleotide mutations
- there are two scripts gemini-code-1793024519427.r and gemini-code-1793024521722.r that attempted to analyze the sequences

### Goal of the project
- we need to do cluster anlysis and divide sequences into clusters based on their length and number of mutations
- the number of mutations should be measured as difference from consensus sequence
- the hypothesis is that we can divide the sequences into two clusters L1, L2, L3 should be similar to each other and in one cluster whereas L4 and L5 should be in second cluster
- do the cluster analysis and confirm or reject the hypothesis
- output of the script should be plot

### info
- use local virtual environment for downloading libraries
