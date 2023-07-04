
<p align="center">
  <img src="./logo.png" alt="Italian Trulli" style="width:70%">
</p>


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.0.2.6-blue) [![GitHub Downloads](https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download)](https://github.com/Kuanhao-Chao/splam/releases) ![os](https://img.shields.io/badge/platform-macOS_/Linux_/Windows-green.svg) [![](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb)

<!-- <a href="https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a> -->


<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/minimap2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/minimap2) -->
<!-- [![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat)](https://pypi.python.org/pypi/mappy) -->



<!-- [![Travis build status](https://travis-ci.org/roblanf/sangeranalyseR.svg?branch=master)](https://travis-ci.org/roblanf/sangeranalyseR)  -->
<!-- ![documentation build status](https://readthedocs.org/projects/pip/badge/) -->

splam is a splice junction recognition model based on a deep grouped residual convolutional neural network that offers fast and precise assessment of splice junctions. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.

There are two main use case scenarios:

1. improving your **alignmnet file**. splam evaluates the quality of splice alignments and removes those that contain spurious splice junctions. This removal process significantly enhances the quality of the downstream transcriptome assembly [[Link](https://ccb.jhu.edu/splam/content/alignment_evaluation.html#alignment-detailed-section)].

2. evaluating the quality of introns in your **annotation file or assembled transcripts** [[Link](https://ccb.jhu.edu/splam/content/annotation_evaluation.html#annotation-detailed-section)].

<br>


## <a name="documentation"></a>Documentation
 📒 The full documentation is on [ccb](http://ccb.jhu.edu/splam/)


<br>

## <a name="installation"></a>Installation

You can install splam with conda package manager. This is the easiest approach.
``` bash
$ conda install -c bioconda liftoff
```


splam is on [PyPi](https://pypi.org/) now. Check out all the release [here](https://pypi.org/manage/project/splam/releases/).
```bash
$ pip install splam
```

You can also install splam from source
```bash
$ git clone https://github.com/Kuanhao-Chao/splam --recursive

$ cd splam/src/

$ python setup.py install
```

<br>

## <a name="quick_start"></a>Quick Start

The simplest example in just three lines of code!

Check this example on google Colab [![](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb)


```
$ cd test

$ splam extract -P -o tmp_out SRR1352129_chr9_sub.bam

$ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed

$ splam clean -o tmp_out
```

<br>


## <a name="publication"></a>Publications
TBC