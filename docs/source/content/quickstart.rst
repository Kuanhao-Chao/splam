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

            console.log("document.getElementsByClassName('sidebar_wse')[0].src: ", document.getElementsByClassName('sidebar_wse')[0].src);
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

Quick Start Guide
=================

This page provides simple quick-start information for using splam with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load splam.

|

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

There are two use case scenarios of Splam. The first one is :ref:`running with an alignment file <splam_bam_quick>`, and second one is :ref:`running with an annotation file <splam_gff_quick>`. Both can be done in three lines of code. Following are the examples:

|

.. _splam_bam_quick:
Cleaning up alignment files  (:code:`BAM`)
-------------------------------------------

The most minimal example gets the job done in three lines of code. More details below.

.. code-block:: bash

    $ cd test

    $ splam extract -P SRR1352129_chr9_sub.bam -o tmp_out_alignment

    $ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_alignment tmp_out_alignment/junction.bed

    $ splam clean -o tmp_out_alignment

| 

.. _splam_gff_quick:
Evaluation annotation files / assembled transcripts (:code:`GFF`)
----------------------------------------------------------------------

.. code-block:: bash

    $ cd test

    $ splam extract refseq_40_GRCh38.p14_chr_fixed.gff -o tmp_out_annotation

    $ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

    $ splam clean -o tmp_out_annotation

|


For more detailed analysis steps, please check :

.. seealso::
    
    * :ref:`alignment-detailed-section`

    * :ref:`annotation-detailed-section`


|
|
|
|
|

.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center

