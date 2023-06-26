.. _alignment-detailed-section:

Alignment file evaluation & cleanup (:code:`GFF`)
=================================================

Image you are an user run RNA-Seq analysis. You ran quality control, aligned reads onto genome. Popular spliced aligners are :code:`HISAT2`, :code:`STAR`, and so on. 

If you checked the alignments on IGV, you will observe there are some read that are spliced aligned on the reference genome crossing different gene loci or at the intergenic regions. The question arised is that are these spliced alignment correct? Are these reads actually come from here? Are there any bias of the splice aligner that aligns reads to this position?

However, there are no tools that assess the outputs of the aligners. Therefore, 

|

**We propose a new step in standard RNA-Seq analysis pipeline. You can run splam!**

|


You can run splam on alignment files. 
splam outputs the scores of every donor and acceptor sites, by using these scores,



.. _alignment-prepareintput:

Step 1: Preparing your input files
+++++++++++++++++++++++++++++++++++

Step 2: Extracting splice junctions in your alignment file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Step 3: Scoring extracted splice junctions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Step 4: Cleaning up your alignment file
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Step 5: Visualization & reports
+++++++++++++++++++++++++++++++++++

|
|
