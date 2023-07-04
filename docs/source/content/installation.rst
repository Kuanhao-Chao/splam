.. raw:: html

  <script type="text/javascript">
    var observer = new MutationObserver(function(mutations) {
        const dark = document.body.dataset.theme == 'dark';
        console.log(dark);
        document.getElementsByClassName('mainlogo')[0].src = dark ? '../_images/jhu-logo-white.png' : "../_images/jhu-logo-dark.png";
        console.log(document.getElementsByClassName('mainlogo')[0].src);
    })
    observer.observe(document.body, {attributes: true, attributeFilter: ['data-theme']});
    console.log(document.body);
  </script>
  <link rel="preload" href="../_images/jhu-logo-dark.png" as="image">


.. image:: ../image/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, mainlogo
   :align: center

|

Installation
===============

System requirements
-------------------

   * Python >= 3.6.0
   * pytorch >= 1.12.0
   * pybedtools >= 0.9.0
   * gffutils >= 0.10.0
   * pybind11 >= 2.10.0

|

There are three ways that you can install splam:

|

Install through conda
-------------------------------

Install splam through conda is the easiest way to go

.. code-block:: bash
   
   $ conda install -c bioconda liftoff
   
   TBC

|

Install through pip
-------------------------

splam is on `PyPi 3.12 <https://pypi.org/project/splam/>`_ now. Check out all the release `here <https://pypi.org/manage/project/splam/releases/>`_.

.. code-block:: bash
   
   $ pip install splam

|

Install from source
-------------------------

You can also install splam from source

.. code-block:: bash

   $ git clone https://github.com/Kuanhao-Chao/splam --recursive

   $ cd splam/src/

   $ python setup.py install

|


Now, you are ready to go !
--------------------------
Please continue to the :ref:`Quick Start Guide`.
