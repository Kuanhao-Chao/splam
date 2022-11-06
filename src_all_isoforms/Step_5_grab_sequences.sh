#!/bin/bash

source constants.py

CLr=$((CL_max/2))
CLl=$(($CLr+1))
# First nucleotide not included by BEDtools

# cat $splice_table | awk -v CLl=$CLl -v CLr=$CLr '{if(($5-CLl) >= 0 ){print $3"\t"($5-CLl)"\t"($6+CLr)}}' > temp.bed
cat $splice_table | awk -v CLl=$CLl -v CLr=$CLr '{print $3"\t"($5-CLl)"\t"($6+CLr)}' > temp.bed

bedtools getfasta -bed temp.bed -fi $ref_genome -fo $sequence -tab