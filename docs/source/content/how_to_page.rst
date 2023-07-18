.. raw:: html

    <script type="text/javascript">
        var observer = new MutationObserver(function(mutations) {
            var dark = document.body.dataset.theme == 'dark';

            if (document.body.dataset.theme == 'auto') {
                dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
            }
            
            console.log(dark);
            document.getElementsByClassName('header-image')[0].src = dark ? '../_images/jhu-logo-white.png' : "../_images/jhu-logo-dark.png";
            document.getElementsByClassName('sidebar_ccb')[0].src = dark ? '../_images/JHU_ccb-white.png' : "../_images/JHU_ccb-dark.png";
            document.getElementsByClassName('sidebar_wse')[0].src = dark ? '../_images/JHU_wse-white.png' : "../_images/JHU_wse-dark.png";



            for (let i=0; i < document.getElementsByClassName('summary-title').length; i++) {
                console.log("document.getElementsByClassName('sidebar_wse')[0].src: ", document.getElementsByClassName('summary-title')[i]);

                if (dark) {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-dark font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-dark text-left docutils";
                } else {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-light font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-light text-left docutils";
                }
            }

        })
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>
    <link rel="preload" href="../_images/jhu-logo-dark.png" as="image">



.. raw:: html
    
    <script type="text/javascript">
        var block_to_insert ;
        var container_block ;
        
        block_to_insert = document.createElement( 'div' );
        block_to_insert.innerHTML = '<img alt="My Logo" style="width:80%;  margin:10px; padding-top:30px" class="logo sidebar_ccb align-center" src="../_images/JHU_ccb-dark.png"><img alt="My Logo" class="logo sidebar_wse align-center" style="width:80%;  margin:10px" src="../_images/JHU_wse-dark.png">' ;
        
        container_block = document.getElementsByClassName( 'sidebar-sticky' )[0];
        console.log("container_block: ", container_block);
        container_block.appendChild( block_to_insert );
    </script>


|


Q & A ...
==========

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

    Given an input of length 800bp, Splam outputs a Tensor with dimensions 3 * 800. The first channel represents the "acceptor scores", the second channel represents the "donor scores", and the third channel represents the "non-splice site scores". Each score is between 0 and 1, representing Splam's confidence in a given site being a splice site. A score closer to one indicates a higher level of confidence in its classification.

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
|
|
|

.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center
