Q & A ...
==========

+++++

Q: What is splam?
-------------------------------------------

.. <div style="padding-left:20px">

    
.. dropdown:: Q: What is splam?
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    splam means two things: (1) splam is a deep residual convolutional neural networks that accurately predict splice junctions based solely on an input DNA sequence, and (2) it also stands for this software that evaluates the annotation files and clean up the alignment files.

|


Q: Why do we need splam?
-------------------------------------------

.. dropdown:: Ans:
    :animate: fade-in-slide-down
    :container: + shadow
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    Given an alignment file, 

    1. We need a good way to evaluate the annotation file and the alignment file.

    2. SpliceAI is biased toward canonical transcript.
    
    3. We designed a biologically realistic model. 



|

Q: What makes splam different from spliceAI?
-------------------------------------------

.. dropdown:: Ans:
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
    
       * **splam** was trained using a high-quality dataset of human donor and acceptor sites. We curated
    
       * **SpliceAI** was trained with canonical transcripts only, and it does not consider alternative splicing.



| 

Q: What is the model architecture of splam?
-----------------------------------------

.. raw:: html

    <details>
    <summary style="font-size:25px; font-weight: bold; padding-left:50px">Ans:</summary>
    <p style="padding-left:50px"></p>
    </details>

| 

Q: How is splam trained?
--------------------------------

.. raw:: html

    <details>
    <summary style="font-size:25px; font-weight: bold; padding-left:50px">Ans:</summary>
    <pre style="padding-left:50px">lots_of_code = "this text block"</pre>
    </details>

| 

Q: Which mode should I run splam, :code:`cpu`, :code:`cuda`, or :code:`mps`?
-------------------------------------------------------------------------------

.. raw:: html

    <details>
    <summary style="font-size:25px; font-weight: bold; padding-left:50px">Ans:</summary>
    <pre style="padding-left:50px">lots_of_code = "this text block"</pre>
    </details>

| 

Q: How do I interpret splam scores?
-------------------------------------

.. raw:: html

    <details>
    <summary style="font-size:25px; font-weight: bold; padding-left:50px">Ans:</summary>
    <pre style="padding-left:50px">lots_of_code = "this text block"</pre>
    </details>

|

Q: What is canonical transcripts? 
------------------------------------------


|

Q: What is alternative splicing?
------------------------------------------

|
|

.. raw:: html

    <style>
    .text-header {
        font-size: 30px;
    }
    </style>