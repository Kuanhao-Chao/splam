
<p align="center">
  <img src="./logo.png" alt="Italian Trulli" style="width:70%">
</p>



[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![version](https://img.shields.io/badge/version-v.0.2.6-blue) [![GitHub Downloads](https://img.shields.io/github/downloads/Kuanhao-Chao/splam/total.svg?style=social&logo=github&label=Download)](https://github.com/Kuanhao-Chao/splam/releases)![os](https://img.shields.io/badge/platform-macOS_/Linux_/Windows-green.svg) [![Build Status](https://github.com/lh3/minimap2/actions/workflows/ci.yaml/badge.svg)](https://github.com/lh3/minimap2/actions)
<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/minimap2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/minimap2) -->
<!-- [![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat)](https://pypi.python.org/pypi/mappy) -->



<!-- [![Travis build status](https://travis-ci.org/roblanf/sangeranalyseR.svg?branch=master)](https://travis-ci.org/roblanf/sangeranalyseR)  -->
<!-- ![documentation build status](https://readthedocs.org/projects/pip/badge/) -->

splam is a splice junction recognition model based on a deep residual convolutional neural network. It was trained on donor and acceptor pairs combined and focuses on a narrow window of 400 basepairs surrounding each splice site, inspired by the understanding that the splicing process primarily depends on signals within this specific region.

There are two main use case scenarios:

1. Evaluating all splice sites in an annotation file or assembled transcripts [[Link](#annotation_splam)].
2. Evaluating all splice junctions in an alignment file and remove spliced alignment containing spurious splice sites [[Link](#alignment_splam)]. Removing spurious splice alignments surprisingly improves transcriptome assembly.

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



## <a name="getting_started"></a>Getting started

### <a name="annotation_splam"></a>I. Evaluating splice junctions in an annotation file.

The first use case scenario is to evaluate all splice junctions in an annotation file or transcripts assembled by assemblers such as StringTie and Scallop. The acceptable file formats for this analysis are `GFF` and `GTF`. Users can execute this mode in two steps: `extract` and `score`.

#### <a name="annotation_splam_extract"></a> I.I. Extracting splice junctions
In this example, given a `GFF` file, you first run 

```
cd test
splam extract MANE.GRCh38.v1.1.subset.gff
```
splam iterates through the GFF file, extracts all introns in transcripts, and writes their coordinates into a `BED` file. The BED file consists of six columns: `CHROM`, `START`, `END`, `JUNC_NAME`, `INTRON_NUM`, and `STRAND`. Here are a few entries from the BED file:

```
chr9    4849549 4860125 JUNC00000007    3       +
chr9    5923308 5924658 JUNC00000008    6       -
chr9    5924844 5929044 JUNC00000009    8       -
```


#### <a name="annotation_splam_score"></a> I.II. Scoring splice junctions
In this step, the goal is to score all the extracted splice junctions. To accomplish this, you will need three essential files. Firstly, you should have the BED file that was generated in the previous step. Additionally, you will require two additional files: (1) the reference genome, which shares coordinates with the junction BED file [[file link](https://github.com/Kuanhao-Chao/splam/blob/main/test/chr9_subset.fa)], and (2) the splam model, accessible at `model/splam_script.pt` [[file link](https://github.com/Kuanhao-Chao/splam/blob/main/model/splam_script.pt)]. Once you have these files in place, you can proceed to score all the splice junctions by executing the following command:


By default, splam automatically detects your environment and runs in `cuda` mode if CUDA is available. However, if your computer is running macOS, splam will check if `mps` mode is available. If neither `cuda` nor `mps` are available, splam will run in `cpu` mode. You can manually specify the mode using the `-d / --device` argument.

Additionally, you can adjust the batch size using the `-b / --batch-size` argument. We recommend setting a small batch size (default is 10) when running splam in `cpu` mode.


```
splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed
```

After this step, a new BED file is produced, featuring eight columns. Two extra columns, namely `DONOR_SCORE` and `ACCEPTOR_SCORE`, are appended to the file. It is worth noting that any unstranded introns are excluded from the output. (p.s. they might be from unstranded transcripts assembled by StringTie).

```
chr9    4849549 4860125 JUNC00000007    3       +       0.7723698       0.5370769
chr9    5923308 5924658 JUNC00000008    6       -       0.9999831       0.9999958
chr9    5924844 5929044 JUNC00000009    8       -       0.9999883       0.9999949
```


### <a name="alignment_splam"></a> II. Improving spliced alignment
The second use case scenario is to remove spurious spliced alignments from an alignment file in order to enhance the overall alignment quality. The process is illustrated in the following figure:

![splam improve alignment](img/improve_alignment.png)


We propose using splam as a new quality trimming step in the RNA-Seq data analysis pipeline. This step is positioned after splice alignment and before transcriptome assembly. In the figure, the quality trimming step is highlighted by the red box, and a detailed breakdown of the trimming process is presented in Figure b.

For this analysis, the acceptable file format is `BAM`. The input to splam should be a sorted BAM file, and it generates a new, cleaned, and sorted BAM file as output. Users can execute this mode in three steps: `extract`, `score`, and `clean`.


#### <a name="alignment_splam_extract"></a> II.I. Extracting splice junctions
Like in the [previous section](#annotation_splam_extract), the initial step involves extracting all splice junctions from an alignment file. The input alignment file must be in `BAM` format, and it is crucial that the <b>input file is sorted</b>.

One important argument is `-P / --paired`.  This argument should be selected based on the RNA sequencing read type. There are two types of RNA sequencing read types: single-read and paired-end sequencing. For a more detailed explanation, you can refer to [this page](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html).


By default, splam processes alignments without pairing them. If your RNA-Seq sample is single-read, there is no need to set this argument. However, if your RNA-Seq sample is from paired-end sequencing, it is highly recommended to run splam with the `-P / --paired` argument. Otherwise, if an alignment is removed, the flag of its mate will not be unpaired. It is worth noting that it takes longer to pair alignments in the BAM file, but it produces more accurate flags.

The primary outputs for this step is a `BED` file containing the coordinates of each junction and some temporary files. If you only want to extract splice junctions from the BAM file without running the subsequent cleaning step, you can use the `-n / --write-junctions-only` argument to skip writing out temporary files.

```
cd test 

splam extract -P SRR1352129_chr9_sub.bam
```

#### <a name="alignment_splam_score"></a>II.II. Scoring splice junctions
It is the same as [this section](#annotation_splam_score). In this step, the input is the BED file generated from the previous stage, and the output is a new BED file that includes the scores for donor and acceptor sites.

```
splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed
```

#### <a name="alignment_splam_clean"></a>II.III. Cleaning up alignment file
After scoring every splice junction in your alignment file, the final step of this analysis is to remove alignments with low-quality splice junctions and update 'NH' tag and flags for multi-mapped reads. You can pass the directory path to splam using the `clean` mode, which will output a new cleaned and sorted BAM file. The implementation of this step utilizes the core functions of `samtools sort` and `samtools merge`. If you want to run this step with multiple threads, you can set the `-@ / --threads` argument accordingly.

```
splam clean -o tmp_out -@ 5
```

<br>


## <a name="manual"></a>splam manual

There are three mode of splam, which are `extract`, `score`, and `clean`. 

```
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
```


### <a name="junction_extract"></a>Extracting splice junctions

```
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
```

### <a name="score_splam"></a>Scoring splice junctions 

```
usage: splam score [-h] [-V] [-o DIR] -G REF.fasta -m MODEL.pt junction_BED

positional arguments:
  junction_BED          target splice junctions in bed files.

optional arguments:
  -h, --help            show this help message and exit
  -V, --verbose
  -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
  -G REF.fasta, --reference-genome REF.fasta
                        The path to the reference genome.
  -m MODEL.pt, --model MODEL.pt
                        the path to the SPLAM! model
```

## <a name="clean_splam"></a>Cleaning up BAM file 

```
usage: splam clean [-h] [-@ threads] [-t threshold] -o DIR

optional arguments:
  -h, --help            show this help message and exit
  -@ threads, --threads threads
                        Set number of sorting, compression and merging threads. By default, operation is single-threaded.
  -t threshold, --threshold threshold
                        The cutoff threshold for identifying spurious splice junctions.
  -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed".
```

<br>



## <a name="example"></a>Example
Here is a simple reproducible example with toy data:

```
cd test

splam extract -P -o tmp_out SRR1352129_chr9_sub.bam 

splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed

splam clean -o tmp_out
```

<br>


## <a name="m_architecture"></a>Model Architecture
![model_architecture](img/splam_architecture.png)


<br>



## <a name="publication"></a>Publications
**Here is the link to the Colab directory: https://colab.research.google.com/drive/1y5q_OjZcjVRwZB78cZ9xLCwRODOiG_eB**
