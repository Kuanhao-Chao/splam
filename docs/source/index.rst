.. raw:: html

    <script type="text/javascript">
        var observer = new MutationObserver(function(mutations) {
            var dark = document.body.dataset.theme == 'dark';

            if (document.body.dataset.theme == 'auto') {
                dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
            }
            
            console.log(dark);
            document.getElementsByClassName('header-image')[0].src = dark ? './_images/jhu-logo-white.png' : "./_images/jhu-logo-dark.png";
            document.getElementsByClassName('sidebar_ccb')[0].src = dark ? './_images/JHU_ccb-white.png' : "./_images/JHU_ccb-dark.png";
            document.getElementsByClassName('sidebar_wse')[0].src = dark ? './_images/JHU_wse-white.png' : "./_images/JHU_wse-dark.png";

            console.log("document.getElementsByClassName('sidebar_wse')[0].src: ", document.getElementsByClassName('sidebar_wse')[0].src);
        })
        observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
        console.log(document.body);
    </script>
    <link rel="preload" href="./_images/jhu-logo-dark.png" as="image">



.. raw:: html
    
    <script type="text/javascript">
        var block_to_insert ;
        var container_block ;
        
        block_to_insert = document.createElement( 'div' );
        block_to_insert.innerHTML = '<img alt="My Logo" style="width:80%;  margin:10px; padding-top:30px" class="logo sidebar_ccb align-center" src="./_images/JHU_ccb-dark.png"><img alt="My Logo" class="logo sidebar_wse align-center" style="width:80%;  margin:10px" src="./_images/JHU_wse-dark.png">' ;
        
        container_block = document.getElementsByClassName( 'sidebar-sticky' )[0];
        console.log("container_block: ", container_block);
        container_block.appendChild( block_to_insert );
    </script>


|


splam's tutorial
*************************

.. _splam_logo:
.. figure::  ./image/logo.png
   :align:   center
   :scale:   20 %

|

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg

.. image:: https://img.shields.io/badge/version-v.0.2.6-blue

.. image:: https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download

.. image:: https://img.shields.io/badge/platform-macOS_/Linux_/Windows-green.svg


Why splam?
==================

splam is a splice junction recognition model based on a deep grouped residual convolutional neural network that offers fast and precise assessment of splice junctions. 

There are two primary applications of splam:

.. splam is useful if you want to :

1. improve your **alignmnet file**. splam evaluates the quality of splice alignments and removes those that contain spurious splice junctions. This removal process significantly enhances the quality of the downstream transcriptome assembly [:ref:`Link <alignment-detailed-section>`].

2. evaluate the quality of introns in your **annotation file or assembled transcripts** [:ref:`Link <annotation-detailed-section>`].



splam is free, it's open source, it's a light weight deep learning model, and it's in Python!

.. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.



|

Main features
=============
* **Biologically inspired training process**: splam was trained on combined donor and acceptor pairs, with a focus on a narrow window of 400 base pairs surrounding each splice site. This approach is inspired by the understanding that the splicing process predominantly relies on signals within this specific region.
.. * **Visualization**: splam produces high-quality figures ready for your publication.
* **Python + C++ integration**: We have taken care of all the engineer work for you! splam is easy to install and runs efficiently due to its underlying C++ implementation. You can install and run splam with just one simple command!
* **Run splam in three steps**: With just three lines of code, you can obtain a new alignment file that is cleaned and sorted.
* **Pytorch implementation**: splam is implemented and trained using the popular Pytorch framework.

|

What splam **doesn't** do
==================================

One feature that splam does not have is the ability to scan through the genome. Many splice site transcriptome tools take the DNA sequence and predict the splice sites within it, such as `SpliceAI <https://github.com/Illumina/SpliceAI>`_. However, their training step only considers splice sites in the canonical transcripts while disregarding the isoforms. This raises the question of whether their deep learning model learns the correct splice junction pattern. Therefore, splam zooms in on predicting at the "splice junction level" rather than the "transcript-level." Splam considers paired splice sites within limited windows of 400bp and provides information about the quality of the splice junctions.

|

.. User Manual
.. ===========
.. If you are already familiar with splam and want to have a quick look at function signatures, please refer to `sangeranalyseR user manual <https://bioconductor.org/packages/devel/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_

.. |

User support
============
Please go through the :ref:`Documentation` below first. If you have questions about using the package, a bug report, or a feature request, please use the GitHub issue tracker here:

https://github.com/Kuanhao-Chao/splam/issues


|

Key contributors
================

splam deep residual convolutional neural network was trained using the PyTorch framework by Kuan-Hao Chao. Kuan-Hao Chao also implemented the package that applies splam to evaluate annotation files and clean up alignment files.

|

Table of content
==================

.. toctree::
   :maxdepth: 2

   content/installation
   content/quickstart
   content/alignment_evaluation
   content/annotation_evaluation
   content/behind_scenes
   content/how_to_page
   content/function_manual
   content/license
   content/contact



|
|
|
|

.. image:: ./image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center


.. raw:: html

    <footer align="center" style="margin-top:-5px">&copy; Copyright 2023, Kuan-Hao Chao</footer> 
