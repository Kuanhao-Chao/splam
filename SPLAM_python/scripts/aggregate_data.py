# aggregate the outputs from SPLAM to a single .csv file

import os
import csv
from progress.bar import Bar
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
#from Bio import SeqIO
from pyfaidx import Fasta 
#from itertools import groupby


# use the score file as the input file, then check with the original GRCh28.p14 fasta file for the donor and acceptor dimers
# cross-reference with the INPUT.fa file to see why some of the entries get deleted -> write those to a new file?
# if you set first line to track_name=junctions, you can visualize better using igv

def run_aggregator(db_name):

    # collect the data to be aggregated from multiple sources
    input_dir = './outputs/' + db_name + '/'
    score_file = input_dir + db_name + '.score.bed'
    df = collect_data(score_file)

    # write output to .csv file
    output_file = input_dir + db_name + ''
    write_output(df, output_file)


def collect_data(score_file):
    # read score data from score file
    df = pd.read_csv(score_file, sep='\t', header=0, \
        names=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score','strand', 'donorScore', 'acceptorScore'])

    # add new columns for the dimers
    df['donorDimer'] = ''
    df['acceptorDimer'] = ''

    # ping the GRCh38.p14 fasta file using pyfaidx to obtain the sequences as a dictionary
    genes = Fasta('../Dataset/GRCh38_p14_genomic_ucsc.fa')
    print(f'Found {len(genes.keys())} unique genes in GRCh38.p14 dataset.')

    # iterate over every coordinate and extract corresponding dimer from fasta
    pbar = Bar('[Info] Getting dimers...', max=len(df))
    for index, row in df.iterrows():
        # obtain the start and end sequences from the 
        chromosome = row['chrom'].strip('>')
        row['chrom'] = chromosome # fixing some naming error in the input file
        donor_start = int(row['chromStart(donor)'])
        acceptor_end = int(row['chromEnd(acceptor)'])

        # extract the dimers from the fasta file
        donor_dimer = str(genes[chromosome][donor_start:donor_start+2]).upper()
        acceptor_dimer = str(genes[chromosome][acceptor_end-2:acceptor_end]).upper()

        # insert the dimers into the pandas dataframe
        df.at[index, 'donorDimer'] = donor_dimer
        df.at[index, 'acceptorDimer'] = acceptor_dimer

        pbar.next()

    pbar.finish()
    print(df.head())
   
    return df

def write_output(df, output_file):
    
    df.to_csv(output_file, index=False)

    pass

# to elucidate the reason why there is a discrepancy between 
def compare_error(df, input_fa_file):
    pass


def stats(df):
    pass


def visualize():
    pass


if __name__ == '__main__':
    dbs = ['chess3', 'gencode_all', 'MANE', 'refseq']
    nums = [0] # CHANGEME

    for num in nums:
        run_aggregator(dbs[num])