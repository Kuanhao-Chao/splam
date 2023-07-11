
# Generalization Tests:

- all outputs from every step of the data processing will be saved in the corresponding folder of the step `{#}_output/`
- `{name}` refers to the name of the gene database

## Positive Dataset:


## Negative Dataset:

1. Generate and pre-process the data:
    
    You will need to randomly generate the negative splice junction dataset. This works by taking the existing protein-coding genes, then selecting the opposite strand to guarantee unique sequences, and creating pseudo-splice-junctions from GT-AG pairs found on this strand.
    
    1. To start, retrieve the protein-coding genes from all four genomes, making use of the `gffutils` library: 

            $ python 1_get_genes.py

        Inputs:
        - `{name}.gff` = annotations corresponding to a genome
        - `{name}_genomic.fa` = genome's fasta file 
        - `{name}_annotation_report.txt` = annotation report 
        *These were downloaded together from the NCBI Genome Database*

        Outputs:
        - `databases/{name}.db` = sqlite3-style databases parsed from the `.gff` annotation files
        - `{name}.bed` = protein-coding genes

    2. Then, generate the dataset of splice junctions, and process into a format readable by Splam. 
        
            $ python 2_extract_all_splam.py

        Inputs:
        - `{name}.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `donor.bed`, `acceptor.bed` = bed files referring to the specific 400nt donor and acceptor sequences
        - `donor_seq.fa`, `acceptor_seq.fa`= fasta files containing the 400nt sequences
        - `d_a.bed` = bed file referring to the intron coordinates of the splice junction (midpoint of the donor and acceptor coords)
        - `input_neg_random.fa` = fasta file containing the 800nt sequence that is given to Splam

    3. Now process the same splice junctions into a format readable by SpliceAI (it will take a subset of the junctions for efficiency).

            $ python 3_extract_all_spliceai.py
        
        Inputs:
        - `d_a.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

        Outputs:
        - `{name}_coords.bed` = bed file referring to the start and end positions of the *whole* SpliceAI input
        - `{name}_seq_noN.fa` = fasta file containing the SpliceAI input, with the full flanking sequence (coords refer to splice junction)
        - `{name}_seq_N.fa` = fasta file containing the SpliceAI input, with repeating N flanking sequence (coords refer to splice junction)


2. Run Splam and SpliceAI on the processed transcripts, and record the scores

    4. Run Splam.

            $ ./4_splam_runner.sh
        
        Inputs:
        - `input_neg_random.fa` 

        Outputs:
        - `{name}.score.bed` = bed file containing the Splam-scored splice junctions


    5. Run SpliceAI. Depending on your system, this may take several days, so you can run each dataset separately: 

            $ ./5_spliceai_prediction_wrapper.sh {name}
        
        Inputs:
        - `{name}_seq_noN.fa`
        - `{name}_seq_N.fa`

        Outputs:
        *There are 5 model output folders, each containing 4 folders with the database names*
        - `spliceai_all_seq.name.noN.{name}.tsv`, `spliceai_all_seq.name.N.{name}.tsv` = names and identifiers for the scored splice junctions
        - `spliceai_all_seq.score.a.noN.{name}.tsv`, `spliceai_all_seq.score.a.N.{name}.tsv` = acceptor site scores for every nt in sequence
        - `spliceai_all_seq.score.d.noN.{name}.tsv`, `spliceai_all_seq.score.d.N.{name}.tsv` = donor site scores for every nt in sequence
        - `spliceai_all_seq.score.n.noN.{name}.tsv`, `spliceai_all_seq.score.n.N.{name}.tsv` = neutral (neither) scores for every nt in sequence

    6. Post-process the Splam and SpliceAI scores into a single file for comparison:

            


## Plotting Result:

### Within a set (positive or negative)

1. Comparative score distributions on histogram

2. Intron lengths vs. scores scatterplot with marginal distributions

### Across both sets

3. ROC/PR curves

4. DT plots


