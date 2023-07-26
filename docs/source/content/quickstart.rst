.. raw:: html

    <script type="text/javascript">

        let mutation_fuc = function(mutations) {
            var dark = document.body.dataset.theme == 'dark';

            if (document.body.dataset.theme == 'auto') {
                dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
            }
            
            document.getElementsByClassName('sidebar_ccb')[0].src = dark ? '../_static/JHU_ccb-white.png' : "../_static/JHU_ccb-dark.png";
            document.getElementsByClassName('sidebar_wse')[0].src = dark ? '../_static/JHU_wse-white.png' : "../_static/JHU_wse-dark.png";



            for (let i=0; i < document.getElementsByClassName('summary-title').length; i++) {
                console.log(">> document.getElementsByClassName('summary-title')[i]: ", document.getElementsByClassName('summary-title')[i]);

                if (dark) {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-dark font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-dark text-left docutils";
                } else {
                    document.getElementsByClassName('summary-title')[i].classList = "summary-title card-header bg-light font-weight-bolder";
                    document.getElementsByClassName('summary-content')[i].classList = "summary-content card-body bg-light text-left docutils";
                }
            }

        }
        document.addEventListener("DOMContentLoaded", mutation_fuc);
        var observer = new MutationObserver(mutation_fuc)
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>
    <link rel="preload" href="../_images/jhu-logo-dark.png" as="image">


|

Quick Start Guide
=================

This page provides simple quick-start information for using Splam with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load Splam.

|

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

There are 2 use case scenarios of Splam. The first one is :ref:`running with an alignment file <splam_bam_quick>`, and second one is :ref:`running with an annotation file <splam_gff_quick>`. Both can be done in three lines of code. We demonstrate below:

|

.. _splam_bam_quick:
Cleaning up alignment files  (:code:`BAM`)
-------------------------------------------

The most minimal example gets the job done in three lines of code:

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


.. image:: ../_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ../_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center