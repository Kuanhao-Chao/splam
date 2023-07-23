
<p align="center">
  <img src="./logo.png" alt="Italian Trulli" style="width:85%">
</p>


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.1.0.0-blue) [![GitHub Downloads](https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download)](https://github.com/Kuanhao-Chao/splam/releases) ![os](https://img.shields.io/badge/platform-macOS_/Linux_/Windows-green.svg) [![](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb)

<!-- <a href="https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a> -->


<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/minimap2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/minimap2) -->
<!-- [![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat)](https://pypi.python.org/pypi/mappy) -->



<!-- [![Travis build status](https://travis-ci.org/roblanf/sangeranalyseR.svg?branch=master)](https://travis-ci.org/roblanf/sangeranalyseR)  -->
<!-- ![documentation build status](https://readthedocs.org/projects/pip/badge/) -->

Splam is a splice junction recognition model based on a deep grouped residual convolutional neural network that offers fast and precise assessment of splice junctions. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.

There are two main use case scenarios:

1. Improving your **alignment file**. Splam evaluates the quality of splice alignments and removes those that contain spurious splice junctions. This removal process significantly enhances the quality of the downstream transcriptome assembly [[Link](https://ccb.jhu.edu/splam/content/alignment_evaluation.html#alignment-detailed-section)].

2. Evaluating the quality of introns in your **annotation file or assembled transcripts** [[Link](https://ccb.jhu.edu/splam/content/annotation_evaluation.html#annotation-detailed-section)].

<br>


## <a name="documentation"></a>Documentation<a class="headerlink" href="#documentation" title="Permalink to this heading">#</a>
📒 The full user manual is available **[here](http://ccb.jhu.edu/splam/)**
<section id="table-of-content" class="">
<h2>Table of content<a class="headerlink" href="#table-of-content" title="Permalink to this heading">#</a></h2>
<div class="toctree-wrapper compound">
<ul>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html">Installation</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#system-requirements">System requirements</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#install-through-pip">Install through pip</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#install-from-source">Install from source</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#install-through-conda">Install through conda</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#check-splam-installation">Check Splam installation</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/installation.html#now-you-are-ready-to-go">Now, you are ready to go !</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/quickstart.html">Quick Start Guide</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/quickstart.html#super-quick-start-3-lines-of-code">Super-Quick Start (3 lines of code)</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html">Alignment file evaluation &amp; cleanup (<code class="code docutils literal notranslate"><span class="pre">BAM</span></code>)</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#introduction">Introduction</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#workflow-overview">Workflow overview</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-1-preparing-your-input-files">Step 1: Preparing your input files</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-2-extracting-splice-junctions-in-your-alignment-file">Step 2: Extracting splice junctions in your alignment file</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-3-scoring-extracted-splice-junctions">Step 3: Scoring extracted splice junctions</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-4-cleaning-up-your-alignment-file">Step 4: Cleaning up your alignment file</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-5-igv-visualization">Step 5: IGV visualization</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#step-6-assembling-alignments-into-transcripts">Step 6: Assembling alignments into transcripts</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/alignment_evaluation.html#what-s-next">What's next?</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html">Annotation file / assembeled transcripts evaluation (<code class="code docutils literal notranslate"><span class="pre">GFF</span></code>)</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html#step-1-preparing-your-input-files">Step 1: Preparing your input files</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html#step-2-extracting-introns-in-your-annotation-file">Step 2: Extracting introns in your annotation file</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html#step-3-scoring-extracted-introns">Step 3: Scoring extracted introns</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html#step-4-evaluating-isoforms-by-splam-scores">Step 4: Evaluating isoforms by Splam scores</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/annotation_evaluation.html#what-s-next">What's next?</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/behind_scenes.html">Behind the scenes of Splam</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/behind_scenes.html#data-curation">Data curation</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/behind_scenes.html#model-architecture">Model architecture</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/behind_scenes.html#splam-training-testing">Splam training &amp; testing</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/behind_scenes.html#reference">Reference</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/generalization.html">Splam generalizes on no-human species</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/generalization.html#test">Test</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/how_to_page.html">Q &amp; A ...</a></li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/function_manual.html">User Manual</a><ul>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/function_manual.html#splam">splam</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/function_manual.html#splam-extract">splam extract</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/function_manual.html#splam-score">splam score</a></li>
<li class="toctree-l2"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/function_manual.html#splam-clean">splam clean</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/license.html">License</a></li>
<li class="toctree-l1"><a class="reference internal" href="http://ccb.jhu.edu/splam/content/contact.html">Contact</a></li>
</ul>
</div>
<div class="line-block">
<div class="line"><br></div>
<div class="line"><br></div>
<div class="line"><br></div>
<div class="line"><br></div>
<div class="line"><br></div>
</div>
<img alt="My Logo" class="logo header-image only-light align-center" src="_images/jhu-logo-dark.png">
<img alt="My Logo" class="logo header-image only-dark align-center" src="_images/jhu-logo-white.png">
</section>

<br>

## <a name="installation"></a>Installation<a class="headerlink" href="#installation" title="Permalink to this heading">#</a>

You can install Splam with conda package manager. This is the easiest approach.
``` bash
$ conda install -c bioconda splam
```


Splam is on [PyPi](https://pypi.org/) now. Check out all the release [here](https://pypi.org/manage/project/splam/releases/).
```bash
$ pip install splam
```

You can also install Splam from source
```bash
$ git clone https://github.com/Kuanhao-Chao/splam --recursive

$ cd splam/src/

$ python setup.py install
```

<br>

## <a name="quick_start"></a>Quick Start<a class="headerlink" href="#quick-start" title="Permalink to this heading">#</a>

The simplest example uses just three lines of code!

Check this example on Google Colab [![](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kuanhao-Chao/splam/blob/main/notebook/splam_example.ipynb)


### Running Splam to clean up alignment files ('BAM')
```
$ cd test

$ splam extract -P SRR1352129_chr9_sub.bam -o tmp_out_alignment

$ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_alignment tmp_out_alignment/junction.bed

$ splam clean -o tmp_out_alignment
```

### Running Splam to evaluate annotation files / assembled transcripts ('GFF')

```
$ cd test

$ splam extract refseq_40_GRCh38.p14_chr_fixed.gff -o tmp_out_annotation

$ splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

$ splam clean -o tmp_out_annotation
```

<br>


## <a name="publication"></a>Publications<a class="headerlink" href="#publication" title="Permalink to this heading">#</a>
TBC
