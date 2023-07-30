Quick Start Guide
=================

This page provides simple quick-start information for using Splam with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load Splam.

|

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

There are two use case scenarios of Splam. The first one is :ref:`running with an alignment file <splam_bam_quick>`, and second one is :ref:`running with an annotation file <splam_gff_quick>`. Both can be done in three lines of code. 

Before you get started, make sure you have already cloned the :ref:`Splam GitHub repository <install_from_source>`. We demonstrate below:


|

.. _splam_bam_quick:
Cleaning up alignment files  (:code:`BAM`)
-------------------------------------------

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

    $ splam extract MANE.GRCh38.v1.1.subset.gff -o tmp_out_annotation

    $ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

    $ splam clean -o tmp_out_annotation

|

Try Splam on Google Colab
+++++++++++++++++++++++++++++++++++

We created some reproducible and easy-to-run Splam examples on Google Colab. It's a good starting point, so go ahead and check them out!


.. image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb



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