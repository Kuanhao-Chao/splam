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

Installation
===============

System requirements
-------------------

.. admonition:: Software dependency

   * Python >= 3.6.0
   * pytorch >= 1.12.0
   * pybedtools >= 0.9.0
   * gffutils >= 0.10.0
   * pybind11 >= 2.10.0

|

There are three ways that you can install splam:



Install through pip
-------------------------

Splam is on `PyPi 3.12 <https://pypi.org/project/splam/>`_ now. Check out all the release `here <https://pypi.org/manage/project/splam/releases/>`_. Pip automatically resolves and installs any dependencies required by Splam.

.. code-block:: bash
   
   $ pip install splam

|

Install from source
-------------------------

You can also install Splam from source. Check out the latest version on `GitHub <https://github.com/Kuanhao-Chao/splam>`_
!

.. code-block:: bash

   $ git clone https://github.com/Kuanhao-Chao/splam --recursive

   $ cd splam/src/

   $ python setup.py install

|

Install through conda
-------------------------------

Install Splam through conda is the easiest way to go

.. code-block:: bash
   
   $ conda install -c bioconda liftoff
   
   TBC

|


Now, you are ready to go !
--------------------------
Please continue to the :ref:`Quick Start Guide`.



|
|
|
|
|

.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image
   :align: center
