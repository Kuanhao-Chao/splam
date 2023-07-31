.. _installation:

Installation
===============

.. _sys-reqs:

System requirements
-------------------

.. admonition:: Software dependency

   * Python >= 3.6.0
   * pytorch >= 1.12.0
   * pybedtools >= 0.9.0
   * gffutils >= 0.10.0
   * pybind11 >= 2.10.0

|

There are three ways that you can install Splam:


.. _install-through-pip:

Install through pip
-------------------------

Splam is on `PyPi 3.12 <https://pypi.org/project/splam/>`_ now. Check out all the release `here <https://pypi.org/manage/project/splam/releases/>`_. Pip automatically resolves and installs any dependencies required by Splam.

.. code-block:: bash
   
   $ pip install splam

|

.. _install-from-source:

Install from source
-------------------------

You can also install Splam from source. Check out the latest version on `GitHub <https://github.com/Kuanhao-Chao/splam>`_
!

.. code-block:: bash

   $ git clone https://github.com/Kuanhao-Chao/splam --recursive

   $ cd splam/src/

   $ python setup.py install

|

.. _install-through-conda: 

Install through conda
-------------------------------

Installing Splam through conda is the easiest way to go:

.. code-block:: bash
   
   TBC

   $ conda install -c bioconda liftoff

|

.. _check-splam-installation:

Check Splam installation
-------------------------------------

Run the following command to make sure Splam is properly installed:

.. code-block:: bash
   
   $ splam -h


.. dropdown:: Terminal output
    :animate: fade-in-slide-down
    :title: bg-light font-weight-bolder
    :body: bg-light text-left

    .. code-block::

        ====================================================================
        An accurate spliced alignment pruner and spliced junction predictor.
        ====================================================================


        ███████╗██████╗ ██╗      █████╗ ███╗   ███╗
        ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║
        ███████╗██████╔╝██║     ███████║██╔████╔██║
        ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║
        ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║
        ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝

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

.. _installation-complete:

Now, you are ready to go !
--------------------------
Please continue to the :ref:`Quick Start Guide`.



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