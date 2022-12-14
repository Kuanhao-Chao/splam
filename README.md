# SpliceNN

SpliceNN is a deep learning-based splice junction predictor. It takes read alignments in BAM or CRAM format and predict highly accurate exon/intron boundaries.

## 1. Data preprocessing

### Processing positive data
* Step 1 - extracting junctions in a BAM file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/1_GET_BAM_JUNCS
  * In total there are 3 steps.
* Step 2 - extracting junctions in a GFF file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/2_GET_REF_JUNCS
  * In total there are 5 steps.

### Processing negative data
* Step 1 - extracting random (noncanonical) negative splice junctions in a FASTA file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/4_GET_NONCANONICAL_NEG_JUNCS

* Step 2 - extracting random canonical negative junctions in a FASTA file: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/3_GET_CANONICAL_NEG_JUNCS

* Step 3 - extracting 1-alignment supported negative junctions from a BAM file:  https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/5_GET_BAM_NEG_JUNCS


## 2. SpliceAI Models (baseline)
We first reimplemented SpliceAI to observe its results before building our own model, SpliceNN.

* SpliceAI reimplementation: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_seq.
  * This is the fixed implementation of SpliceAI in keras Python2: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src_seq/TRAIN_model_spliceAI_keras.py
  * This is our new implementation of SpliceAI in Pytorch: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src_seq/TRAIN_model_spliceAI_pytorch.py
  * This is the main classes we defined for spliceAI: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src_seq/spliceai_pytorch.py


## 3. SpliceNN Models (our new proposed model)
We implemented our new model, SpliceNN, using deep residual CNN with cardinality `C` in Pytorch.

* This is the script that you can find SpliceNN's loss function and some assessment functions: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_utils.py

* For the main model training script, it can be found at: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_train_model.py

* For the main model classes we defined, it can be found at: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN.py

* And the full implementation of SpliceNN can be found in the src/ directory: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src


## 4. Assessment Metrics

* Visualization: We visualized various different model assessment metrics during the training process. The code used to create theses visualizations can be found here: https://github.com/Kuanhao-Chao/SpliceNN/blob/main/src/SpliceNN_visualization.py

## 5. Downstream biological analysis assessment
* In total there are seven steps, and the script of each step can be found here: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/Experiment_BAM.

## 6. SpliceAI and SpliceNN Comparison
* In total there are seven steps; https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src/Experiment_SpliceAI

## 7. SpliceNN in C++ code

* We wrapped SpliceNN into C++ code. Here is the link to the github: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_SPLAM

## 8. Intron Chain matcher in C++ code

* We implemented a C++ tool to assess the model by checking the matching intron chain because StringTie assembly process will clean up some spurious spliced alignments. Here is the link to the C++ implementation: https://github.com/Kuanhao-Chao/SpliceNN/tree/main/src_intron_matcher
