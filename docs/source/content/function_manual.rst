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

User Manual 
=======================

splam
---------------------------------

.. code-block:: text

   usage: splam [-h] [-v] [-c] {extract,score,clean} ...

   splice junction predictor to improve alignment files (BAM / CRAM)

   optional arguments:
   -h, --help            show this help message and exit
   -v, --version
   -c, --citation

   Commands:
   {extract,score,clean}
      extract             Extracting all splice junctions from an alignment or annotation file
      score               Scoring all splice junctions
      clean               Cleaning up spurious splice alignment


|

splam extract
-----------------------------------

.. code-block:: text

   usage: splam extract [-h] [-V] [-P] [-n] [-f FILE_FORMAT] [-o DIR] [-M DIST] [-g GAP] INPUT

   positional arguments:
   INPUT                 target alignment file in BAM format or annotation file in GFF format.

   optional arguments:
   -h, --help            show this help message and exit
   -V, --verbose         running splam in verbose mode.
   -P, --paired          bundling alignments in "paired-end" mode.
   -n, --write-junctions-only
                           writing out splice junction bed file only without other temporary files.
   -f FILE_FORMAT, --file-format FILE_FORMAT
                           the file type for SPLAM to process. It can only be "BAM", "GFF", or "GTF". The default value is "BAM".
   -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
   -M DIST, --max-splice DIST
                           maximum splice junction length
   -g GAP, --bundle-gap GAP
                           minimum gap between bundles

|

splam score 
-----------------------------------

.. code-block:: text

   usage: splam score [-h] [-V] [-o DIR] [-b BATCH] [-d pytorch_dev] -G REF.fasta -m MODEL.pt junction_BED

   positional arguments:
   junction_BED          target splice junctions in bed files.

   optional arguments:
   -h, --help            show this help message and exit
   -V, --verbose
   -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
   -b BATCH, --batch-size BATCH
                           the number of samples that will be propagated through the network. By default, the batch size is set to 10.
   -d pytorch_dev, --device pytorch_dev
                           the computing device that is used to perform computations on tensors and execute operations in the PyTorch framework. By
                           default, this parameter is detectd automatically.
   -G REF.fasta, --reference-genome REF.fasta
                           The path to the reference genome.
   -m MODEL.pt, --model MODEL.pt
                           the path to the SPLAM! model

|                     

splam clean 
-----------------------------------

.. code-block:: text

   usage: splam clean [-h] [-@ threads] [-t threshold] -o DIR

   optional arguments:
   -h, --help            show this help message and exit
   -@ threads, --threads threads
                           Set number of sorting, compression and merging threads. By default, operation is single-threaded.
   -t threshold, --threshold threshold
                           The cutoff threshold for identifying spurious splice junctions.
   -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed".



|
|
|
|
|

.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center
