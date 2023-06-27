Quick Start Guide
=================

This page provides simple quick-start information for using splam with :code:`BAM` and :code:`GFF` files. Please read the :ref:`alignment-detailed-section` or :ref:`annotation-detailed-section` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load splam.

|

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

The most minimal example gets the job done in three lines of code for two use case scenarios of splam. More details below:

|

Alignment file evalutation (:code:`BAM`)
-------------------------------------------

The most minimal example gets the job done in three lines of code. More details below.

.. code-block:: bash

   cd test

   splam extract -P -o tmp_out SRR1352129_chr9_sub.bam 

   splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed

   splam clean -o tmp_out

| 

Annotation file evalutation (:code:`GFF`)
-------------------------------------------

.. code-block:: bash


   cd test

   splam extract MANE.GRCh38.v1.1.subset.gff

   splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed

|


For more detailed analysis steps, please check :

* :ref:`alignment-detailed-section`

* :ref:`annotation-detailed-section`

|
|