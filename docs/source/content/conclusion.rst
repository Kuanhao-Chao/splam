.. raw:: html

    <script type="text/javascript">
        var observer = new MutationObserver(function(mutations) {
            const dark = document.body.dataset.theme == 'dark';
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


Conclusion
==========

splam is a new deep-learning-based splice junction recognizer that helps us to better understand splice junctions. 


We propose a novel splice junction recognition model called SPLAM, which is a deep convolutional neural
network being able to train on laptop computers. Our approach incorporates several innovative training
concepts. First, we leverage expression information from thousands of samples across different human
tissues in the GTEX project to enhance the quality of splice junctions. Furthermore, we propose a distinctive
method for selecting negative splice junction-like sequences, with a specific focus on those supported by
soley one single splice alignment situated on the opposite strand of the protein coding genes. 


|
|
|
|

.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center


.. raw:: html

    <footer align="center" style="margin-top:-5px">&copy; Copyright 2023, Kuan-Hao Chao</footer> 