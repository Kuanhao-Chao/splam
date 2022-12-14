# SpliceNN

SpliceNN is a deep learning-based splice junction predictor. It takes read alignments in BAM or CRAM format and predict highly accurate exon/intron boundaries.


## Data preprocessing

### Processing positive data
* Step 1 - extracting junctions in a BAM file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/1_GET_BAM_JUNCS
* Step 2 - extracting junctions in a GFF file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/2_GET_REF_JUNCS

### Processing negative data
* Step 1 - extracting random negative splice junctions in a FASTA file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/4_GET_NONCANONICAL_NEG_JUNCS

* Step 2 - extracting random canonical negative junctions in a FASTA file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/3_GET_CANONICAL_NEG_JUNCS

* Step 3 - extracting 1-alignment supported negative junctions from a BAM file:  https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/5_GET_BAM_NEG_JUNCS


### Models & Evaluation
* SpliceAI reimplementation: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_seq.

* This is the script that you can find SpliceNN's loss function: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_utils.py



* For the main model training script, it can be found at: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_train_model.py

* And the full implementation of SpliceNN can be found in the src/ directory: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src


* Visualization: We visualized various different model assessment metrics during the training process. The code used to create theses visualizations can be found here: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_visualization.py

## Downstream biological analysis assessment
* In total there are seven steps, and the script of each step can be found here: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/Experiment_BAM.

## Comparison to baseline model: SpliceAI
* https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/Experiment_SpliceAI

## Data preprocessing

We wrapped SpliceNN into C++ code. Here is the link to the github: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_SPLAM

## Data preprocessing

We implemented a C++ tool to assess the model by checking the matching intron chain because StringTie assembly process will clean up some spurious spliced alignments. Here is the link to the C++ implementation: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_intron_matcher