<!-- <h1 align="center">splam</h1> --> 
![Splam Logo](./logo.png) 

splam is a deep learning-based splice junction recognition model. It was trained on MANE and alternate

It has two main use case scenarios:

1. The first use case scenario is to evaluate all the splice sites in the annotation file.
2. The second one is to evaluate splice junctions in the alignment file and remove spurious spliced alignment. 

It takes read alignments in BAM or CRAM format and predict highly accurate exon/intron boundaries.


<!-- # Table of Contents
- [Table of Contents](#table-of-contents)
- [User's Guide](#users-guide)
  - [Installation](#installation)
- [Model Architecture](#model-architecture) -->

# <a name="getting_started"></a>Getting started

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



## <a name="usage"></a>Usage

### <a name="annotation_splam"></a>Evaluating splice junctions in an annotation file.




### <a name="alignment_splam"></a>Evaluating splice junctions in an alignment file.





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
    extract             Extracting all splice junctions from a BAM file
    score               Scoring all splice junctions
    clean               Cleaning up spurious splice alignment
```


## <a name="junction_extract"></a>Extracting splice junctions

In this step, you input a sorted BAM file. SPLAM extracts all splice junctions from the alignments and writes out a BED file with the coordinates of every junction. You can run it with or without alignment pairing by setting the argument '-P'. The latter option unpairs the alignment's mate if it gets removed.

If you simply want to extract splice junctions in the BAM file without running the following steps, you can run with the `--write-junctions-only` argument to skip writing out temporary files.

```
usage: splam extract [-h] [-V] [-P] [-n] [-o DIR] [-M DIST] [-g GAP] BAM_INPUT

positional arguments:
  BAM_INPUT             target alignment file in BAM format.

optional arguments:
  -h, --help            show this help message and exit
  -V, --verbose         running splam in verbose mode
  -P, --paired          paired-end bundling alignments
  -n, --write-junctions-only
                        only write out splice junction bed file without other temporary files.
  -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
  -M DIST, --max-splice DIST
                        maximum splice junction length
  -g GAP, --bundle-gap GAP
                        minimum gap between bundles
```

## <a name="score_splam"></a>Scoring splice junctions 

After extracting splice junctions, you can run splam with the `score` mode to score each splice junction in the BED file generated from the previous step. Stranded splice junctions will be scored using splam's deep convolutional neural network in batches, and the result of this process will be a new BED file containing the donor and acceptor scores.

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

Finally, you can pass the directory path to splam using the `clean` mode, which will output a new cleaned and sorted BAM file. Alignments with low-quality splice junctions are removed and flags are updated.

```
usage: splam clean [-h] -o DIR

optional arguments:
  -h, --help            show this help message and exit
  -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
```


## <a name="example"></a>Example

Here is a simple reproducible example with toy data:

```
cd test

splam extract -P SRR1352129_chr9_sub.bam

splam score -G chr9.fa -m ../model/splam_script.pt -o tmp_splam_out tmp_splam_out/junction.bed

splam clean -o tmp_splam_out
```



# <a name="m_architecture"></a>Model Architecture
![model_architecture](./splam_architecture.png)


# <a name="installation"></a>Publications


**Here is the link to the Colab directory: https://colab.research.google.com/drive/1y5q_OjZcjVRwZB78cZ9xLCwRODOiG_eB**