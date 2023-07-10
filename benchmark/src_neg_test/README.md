
# Generalization Tests:

## How to run:

** Note: all outputs from every step of the data processing will 

Positive Dataset:


Negative Dataset:

1. Generate and pre-process the data:
    
    You will need to randomly generate the negative splice junction dataset. This works by taking the existing protein-coding genes, then selecting the opposite strand to guarantee unique sequences, and creating pseudo-splice-junctions from GT-AG pairs found on this strand.
    
    (1) To start, retrieve the protein-coding genes from all four genomes, making use of the `gffutils` library: 

    `$ python 1_get_genes.py`

    Inputs:
        - `{name}.gff` = the annotations corresponding to a genome
        - `{name}_genomic.fa` = the genome's fasta file 
        - `{name}_annotation_report.txt` = the annotation report 
        ** These were downloaded together from the NCBI Genome Database

    Outputs:
        - `databases/{name}.db` = the sqlite3-style databases parsed from the `.gff` annotation files
        - `{name}.bed` = the protein-coding genes

    (2) Then, generate the dataset of splice junctions, and process into a format readable by Splam. 
    
    `$ python 2_extract_all_splam.py`

    Inputs:
        - `{name}.bed`
        - `{name}_genomic.fa`
        - `{name}_annotation_report.txt`

    Outputs:
        - 

    (3) Now process the same splice junctions into a format readable by SpliceAI (it will take a subset of the junctions for efficiency).

    `$ python 3_extract_all_spliceai.py`


2. Run the Splam and SpliceAI on the processed transcripts, and record the scores

    (1) Run Splam


    (2) Run SpliceAI. Depending on your system, this may take a few days, so you can run in separate blocks: 
    
    `$ ./5_spliceai_prediction_wrapper.sh`
    

    
    (3) 


Plotting Result: