<!-- <h1 align="center">splam</h1> --> 
![Splam Logo](./logo.png) 


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

#### <a name="annotation_splam_extract"></a> Extracting splice junctions
In this example, given a `GFF` file, you first run 

```
splam extract input.gff
```
splam iterates through the GFF file, extracts all introns in transcripts, and writes their coordinates into a `BED` file. The BED file consists of six columns: `CHROM`, `START`, `END`, `JUNC_NAME`, `INTRON_NUM`, and `STRAND`. Here are a few entries from the BED file:

```
chr9    14940   15080   JUNC00000001    251     -
chr9    14940   15040   JUNC00000002    4       -
chr9    14954   15080   JUNC00000003    3       -
```


#### <a name="annotation_splam_score"></a>Scoring splice junctions
In this step, the goal is to score all the extracted splice junctions. To accomplish this, you will need three essential files. Firstly, you should have the BED file that was generated in the previous step. Additionally, you will require two additional files: (1) the reference genome, which shares coordinates with the junction BED file [[example link](https://github.com/Kuanhao-Chao/splam/blob/main/test/chr9_subset.fa)], and (2) the splam model, accessible at [`model/splam_script.pt`](https://github.com/Kuanhao-Chao/splam/blob/main/model/splam_script.pt). Once you have these files in place, you can proceed to score all the splice junctions by executing the following command:

```
splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed
```

After this step, a new BED file is produced, featuring eight columns. Two extra columns, namely `DONOR_SCORE` and `ACCEPTOR_SCORE`, are appended to the file. It is worth noting that any unstranded introns are excluded from the output. (p.s. they might be from unstranded transcripts assembled by StringTie).

```
chr9    14940   15080   JUNC_0  1       -       0.9999496       0.99997985
chr9    14940   15040   JUNC_1  1       -       0.9967939       0.9991217
chr9    14954   15080   JUNC_2  1       -       0.9030796       0.9573919
chr9    15347   16025   JUNC_4  1       -       5.2404507e-09   1.0171946e-08
```


### <a name="alignment_splam"></a> II. Improving spliced alignment

#### <a name="alignment_splam_extract"></a> Extracting splice junctions
Like in the [previous section](#annotation_splam_extract), the initial step involves extracting all splice junctions from an alignment file. The input alignment file must be in `BAM` format, and it is crucial that the input file is sorted.

One important argument is `-P / --paired`.  This argument should be selected based on the RNA sequencing read type. There are two types of RNA sequencing read types: single-read and paired-end sequencing. For a more detailed explanation, you can refer to [this page](https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html).


By default, splam processes alignments without pairing them. If your RNA-Seq sample is single-read, there is no need to set this argument. However, if your RNA-Seq sample is from paired-end sequencing, it is highly recommended to run splam with the `-P / --paired` argument. Otherwise, if an alignment is removed, the flag of its mate will not be unpaired. It is worth noting that it takes longer to pair alignments in the BAM file, but it produces more accurate flags.

The primary outputs for this step is a `BED` file containing the coordinates of each junction and some temporary files. If you only want to extract splice junctions from the BAM file without running the subsequent cleaning step, you can use the `-n / --write-junctions-only` argument to skip writing out temporary files.

```
splam extract -P -o tmp_out SRR1352129_chr9_sub.bam
```

#### <a name="alignment_splam_score"></a>Scoring splice junctions
It is the same as [this section](#annotation_splam_score). 
```
splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out tmp_out/junction.bed
```

#### <a name="alignment_splam_clean"></a>Cleaning up alignment file
After scoring every splice junction in your alignment file, the final step of this analysis is to remove alignments with low-quality splice junctions and update 'NH' tag and flags for multi-mapped reads. You can pass the directory path to splam using the `clean` mode, which will output a new cleaned and sorted BAM file. 

```
splam clean -o tmp_out
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
    extract             Extracting all splice junctions from a BAM file
    score               Scoring all splice junctions
    clean               Cleaning up spurious splice alignment
```


### <a name="junction_extract"></a>Extracting splice junctions

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
usage: splam clean [-h] -o DIR

optional arguments:
  -h, --help            show this help message and exit
  -o DIR, --outdir DIR  the directory where the output file is written to. Default output filename is "junction_score.bed"
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
![model_architecture](./splam_architecture.png)


<br>


## <a name="publication"></a>Publications


**Here is the link to the Colab directory: https://colab.research.google.com/drive/1y5q_OjZcjVRwZB78cZ9xLCwRODOiG_eB**