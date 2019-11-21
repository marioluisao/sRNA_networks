# sRNA_networks
Inference of bacterial small non-coding RNA networks using the *Inferelator*

This repository contains all the input files, *Inferelator* code and the R scripts used for generating the results reported in the revised version of the paper:

Arrieta-Ortiz ML, Hafemeister C, Shuster B, Baliga NS, Bonneau R, Eichenberger P. 2019. Inference of bacterial small RNA regulatory networks  and integration with transcription factor driven regulatory networks. bioRxiv. doi: https://doi.org/10.1101/657478  
 
Below, we briefly describe the contents of this repository. 

**1.	Input_files:** it contains the files required to run the *Inferelator* in each the four analyzed bacterial species. Files for inferring the *Escherichia coli* network that includes eight sRNAs are in the Ecoli_8sRNAs folder. Similarly,  files for inferring the FsrA, PrrF and RsaE regulons for *Bacillus subtilis*, *Pseudomonas aeruginosa* and *Staphylococcus aureus* are in the corresponding folders. We included the input folders for the instances of shuffled *E. coli* sRNA priors and noisy *E. coli* sRNA priors (with false sRNA priors added). The ratio of true: false sRNA priors (1:1, 1:2, 1:5) is indicated by the n*i* term in last part of the noisy_sRNA_priors_1_\*.csv files. Although we did not infer sRNA regulons for *S. enterica*, we made available the files required to infer a transcriptional regulatory network for this species (sRNA priors are not included). 

Each folder contains an expression matrix (\*.RData), TFs and sRNA priors (also refer to as gold standard-*gs*) and a .csv file with the list of all potential regulators to be considered (TFs and sRNAs). The values in the prior matrix are in the {0, 1, -1} set, where 1 and -1 indicate activation and repression, respectively.

**2.	Inferelator_code:** this folder contains the code used to infer all the reported TF-controlled and sRNA-controlled networks. For more details about how to run the *Inferelator*, please check the README file in this folder. 

This code was previously used for generating the improved transcriptional regulatory network model of *B. subtilis* reported in: Arrieta-Ortiz ML, Hafemeister C, et al. 2015. An experimentally supported model of the *Bacillus subtilis* global transcriptional regulatory network. Mol. Syst. Biol.  11. 

This code is also available at: https://github.com/ChristophH/Inferelator

**3.	Inferelator_output_files:** it contains all the networks inferred with the *Inferelator* (but the ones with noisy *E. coli*  priors). Each output folder contains three files: regression coefficients (betas_\*.RData files), confidence scores (combined_conf\*.RData files) and the parameters and input file (params_and_input.RData).

**4.	Miscellaneous_scripts:** it contains multiple R scripts that can be used to generate input files and to analyze the output of the *Inferelator*. 


