CL_max=10000
# Maximum nucleotide context length (CL_max/2 on either side of the 
# position of interest)
# CL_max should be an even number

SL=5000
# Sequence length of SpliceAIs (SL+CL will be the input length and
# SL will be the output length)

splice_table='dataset.txt'
ref_genome='/Users/chaokuan-hao/Documents/Projects/SpliceNN/Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa'
# Input details

data_dir='./'
sequence='sequence.txt'
# Output details