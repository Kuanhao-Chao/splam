.. _annotation-detailed-section:

Annotation file / assembeled transcripts evaluation (:code:`GFF`)
=================================================

We all know that annotation files are very noisy. There are mistakes in human annotation, and not to mention those species with less manual curation process. 

Image you are an user run RNA-Seq analysis. You run Fastqc, aligned reads onto genome, and assembled alignments into a transcripts. Now you get a :code:`GFF` or :code:`GTF` file generated from a transcriptome assembler. The next question you would like to ask is, what is the quality of those assembled transcripts?

Standard pipeline would be run :code:`gffcompare` to compare your assembled transcripts to the gold standard annotation file. 

**We provide another approach. You can run splam!**

The main purpose of running splam on annotation file is that it gives users an overview of what is the quality of all the introns (splice junctions) inside your annotation file or assembled transcripts. 

.. Since how accurate the splice junctions are and 
splam outputs the scores of every donor and acceptor sites. You can (1) have an overview of what's the quality of the assembled transcripts. Furthermore, you can (2) assess each transcript by checking how many bad splice junctions a transcripts have, which can be a good filtering criterion.

There are more potential usage!

Step 1: Preparing your input files
+++++++++++++++++++++++++++++++++++

There are three files you need to get prepared:

1. An annotation file in :code:`GFF` or :code:`GTF` format.
2. A reference genome in :code:`FASTA` format.
3. The splam model, which you can find here: `splam.pt <https://github.com/Kuanhao-Chao/splam/blob/main/model/splam_script.pt>`_

.. code-block:: bash


   cd test

   splam extract MANE.GRCh38.v1.1.subset.gff

   splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed

Step 2: Extracting introns in your annotation file
+++++++++++++++++++++++++++++++++++


Step 3: Scoring extracted introns
+++++++++++++++++++++++++++++++++++


Step 4: Visualization & reports
+++++++++++++++++++++++++++++++++++


|
|
