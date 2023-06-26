.. _annotation-detailed-section:

Annotation file / assembeled transcripts evaluation (:code:`GFF`)
=========================================================================


Image you are an user run RNA-Seq analysis. You ran quality control, aligned reads onto genome, and assembled alignments into a transcripts. Now you get a :code:`GFF` or :code:`GTF` file generated from a transcriptome assembler. The next question you would like to ask is, what is the quality of those assembled transcripts?

Standard pipeline would be run :code:`gffcompare` to compare your assembled transcripts to the reference annotation file, for example, `RefSeq <https://ftp.ncbi.nlm.nih.gov/refseq/>`_, `Gencode <https://www.gencodegenes.org>`_, or `CHESS <http://ccb.jhu.edu/chess/>`_. However, it is based on the assumption that the annotation file is correct, and we are comparing how close our assembled transcripts are to the annotations. We all know that annotation files are very noisy. There are mistakes in human annotation, and not to mention those species with less manual curation process. 

|
 
**We provide another approach to evaluate transcripts. You can run splam!**

|

You can run splam on **(1) annotation files** or **(2) assembled transcripts**. splam outputs the scores of every donor and acceptor sites, by using these scores, 

1. You can get an overview of what is the quality of all the introns (splice junctions) inside your annotation file.


2. You (I) have an overview of what's the quality of the assembled transcripts. Furthermore, you can (II) assess each transcript by checking how many bad splice junctions a transcripts have, which can be a good filtering criteria for assembled transcripts.

3. There are more potential usage!


.. Splice junctions play a critical rule of determine the 

.. The main purpose of running splam is that it 

.. 

.. _annotation-prepareintput:

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
+++++++++++++++++++++++++++++++++++++++++++++++++++++


Step 3: Scoring extracted introns
+++++++++++++++++++++++++++++++++++


Step 4: Visualization & reports
+++++++++++++++++++++++++++++++++++


|
|
