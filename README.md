<!-- <h1 align="center">splam</h1> --> 
![Splam Logo](./logo.png) 

splam is a deep learning-based splice junction predictor. It takes read alignments in BAM or CRAM format and predict highly accurate exon/intron boundaries.



<!-- # Table of Contents
- [Table of Contents](#table-of-contents)
- [User's Guide](#users-guide)
  - [Installation](#installation)
- [Model Architecture](#model-architecture) -->

# <a name="uguide"></a>User's Guide 

## <a name="installation"></a>Installation

```
git clone https://github.com/Kuanhao-Chao/splam --recursive
cd splam/src/
python setup.py install
```

## <a name="example"></a>Example

```
splam extract -P SRR1352129_chr9_sub.bam

splam score -G chr9.fa -m ../model/splam_script.pt -o tmp_splam_out tmp_splam_out/junction.bed

splam clean -o tmp_splam_out%
```

# <a name="getting_started"></a>Getting started

## <a name="junction_extract"></a>Extracting splice junctions

## <a name="junction_extract"></a>Scoring splice junctions 

## <a name="junction_extract"></a>Cleaning up BAM file 


# <a name="m_architecture"></a>Model Architecture
![model_architecture](./splam_architecture.pngig)


# <a name="installation"></a>Publications


**Here is the link to the Colab directory: https://colab.research.google.com/drive/1y5q_OjZcjVRwZB78cZ9xLCwRODOiG_eB**