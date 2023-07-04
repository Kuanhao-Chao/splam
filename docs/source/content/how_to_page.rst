Q & A ...
==========

+++++

.. Q: What is splam?
.. -------------------------------------------


.. <div style="padding-left:20px">
    
.. dropdown:: Q: What is splam?
    :animate: fade-in-slide-down
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    splam means two things: **(1)** splam refers to the deep grouped residual convolutional neural network model that we designed to accurately predict splice junctions based solely on an input DNA sequence, and **(2)** it also stands for this software that and clean up alignment files and evaluate annotation files.

|


.. Q: Why do we need splam?
.. -------------------------------------------

.. dropdown:: Q: Why do we need splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left


    We are concerned about the way of training splice junction predictor simply replying on splice junctions in only canonical transcripts. Designing a splice site recognition method based only on one isoform per gene may result in mis-labeling alternative splice sites even when they are perfectly valid. Therefore, 

        * **we designed a biologically realistic model.** splam was trained on combined donor and acceptor pairs, with a focus on a narrow window of 400 base pairs surrounding each splice site. This approach is inspired by the understanding that the splicing process primarily relies on signals within this specific region.


    Furthermore, there are two applications of splam: 

    When inspecting an alignment file in IGV, it becomes apparent that some reads are spliced and aligned across different gene loci or intergenic regions. This raises the question, "Are these spliced alignments correct?" Therefore,

        * **we need a trustworthy way to evaluate all the spliced alignments in the alignment file.** splam learns splice junction patterns, and we have demonstrated that applying Splam to remove spurious spliced alignments improves transcript assembly! :ref:`alignment evaluation section <alignment-detailed-section>`.

    Additionally, we acknowledge that annotation files are not perfect, and there are more errors in the assembled transcripts. The current approach to assessing assembled transcripts involves comparing them with the annotation.

        * **we can utilize splam to score all introns in transcripts and provide a reference-free evalutation.**  :ref:`annotation evaluation section <annotation-detailed-section>`.



|

.. Q: What makes splam different from spliceAI?
.. -------------------------------------------

.. dropdown:: Q: What makes splam different from spliceAI?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left


    SPLAM and SpliceAI are both frameworks used for predicting splice junctions in DNA sequences, but they have some key differences.


    #. **Input constraints:**
 
       * **splam**: It follows the design principle of using biologically realistic input constraints. It uses a window limited to 200 base pairs on each side of the donor and acceptor sites, totaling 800 base pairs. Furthermore, we pair each donor and acceptor
       .. figure::  ../image/splam_input.png
            :align:   center
            :scale:   40 %
     
       * **SpliceAI**: The previous state-of-the-art CNN-based system, SpliceAI, relies on a window of 10,000 base pairs flanking each splice site to obtain maximal accuracy. However, this window size is much larger than what the splicing machinery in cells can recognize.


    #. **Training data**
    
       * **splam** was trained using a high-quality dataset of human donor and acceptor sites. Check out the :ref:`data curation section <data_curation>`.
    
       * **SpliceAI** was trained with canonical transcripts only, and it does not consider alternative splicing.



| 

.. Q: What is the model architecture of splam?
.. -----------------------------------------


.. dropdown:: Q: What is the model architecture of splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`model architecture section <model_architecture>`.

|

.. Q: What is the model architecture of splam?
.. -----------------------------------------


.. dropdown:: Q: What is the model architecture of splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`model architecture section <model_architecture>`.

| 

.. Q: How is splam trained?
.. --------------------------------

.. dropdown:: Q: How is splam trained?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Check out the :ref:`splam training and testing section <splam_train_test>`.

| 

.. Q: Which mode should I run splam, :code:`cpu`, :code:`cuda`, or :code:`mps`?
.. -------------------------------------------------------------------------------

.. dropdown:: Q: Which mode should I run splam, :code:`cpu`, :code:`cuda`, or :code:`mps`?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

| 

.. Q: How do I interpret splam scores?
.. -------------------------------------

.. dropdown:: Q: How do I interpret splam scores?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

|

.. .. Q: What is canonical transcripts? 
.. .. ------------------------------------------

.. .. dropdown:: Q: What is canonical transcripts? 
..     :animate: fade-in-slide-down
..     :container: + shadow
..     :title: bg-light font-weight-bolder
..     :body: bg-light text-left


.. |

.. .. Q: What is alternative splicing?
.. .. ------------------------------------------

.. .. dropdown:: Q: What is alternative splicing?
..     :animate: fade-in-slide-down
..     :container: + shadow
..     :title: bg-light font-weight-bolder
..     :body: bg-light text-left


|
|
