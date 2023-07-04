
<p align="center">
  <img src="./logo.png" alt="Italian Trulli" style="width:70%">
</p>



[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.0.2.6-blue) [![GitHub Downloads](https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download)](https://github.com/Kuanhao-Chao/splam/releases)![os](https://img.shields.io/badge/platform-macOS_/Linux_/Windows-green.svg) [![Build Status](https://github.com/lh3/minimap2/actions/workflows/ci.yaml/badge.svg)](https://github.com/lh3/minimap2/actions)
<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/minimap2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/minimap2) -->
<!-- [![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat)](https://pypi.python.org/pypi/mappy) -->



<!-- [![Travis build status](https://travis-ci.org/roblanf/sangeranalyseR.svg?branch=master)](https://travis-ci.org/roblanf/sangeranalyseR)  -->
<!-- ![documentation build status](https://readthedocs.org/projects/pip/badge/) -->

splam is a splice junction recognition model based on a deep grouped residual convolutional neural network that offers fast and precise assessment of splice junctions. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.

There are two main use case scenarios:

1. improving your alignmnet file. splam evaluates the quality of splice alignments and removes those that contain spurious splice junctions. This removal process significantly enhances the quality of the downstream transcriptome assembly [ ![Link](https://ccb.jhu.edu/splam/content/alignment_evaluation.html#alignment-detailed-section)].

2. evaluating the quality of introns in your annotation file or assembled transcripts [ ![Link](https://ccb.jhu.edu/splam/content/annotation_evaluation.html#annotation-detailed-section)].

<br>


<!-- # Table of Contents
- [Table of Contents](#table-of-contents)
- [User's Guide](#users-guide)
  - [Installation](#installation)
- [Model Architecture](#model-architecture) -->


## <a name="installation"></a>Installation

You can install splam with conda package manager. This is the easiest approach.
```
TBC
```

You can also install splam from source
```
git clone https://github.com/Kuanhao-Chao/splam --recursive
cd splam/src/
python setup.py install
```

or with pip
```
pip install splam
```

<br>



## <a name="publication"></a>Publications
**Here is the link to the Colab directory: https://colab.research.google.com/drive/1y5q_OjZcjVRwZB78cZ9xLCwRODOiG_eB**
