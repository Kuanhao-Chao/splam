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
    output_file = input_dir + db_name + '.full_data.csv'

    # make csv file if it does not exist, write to pd dataframe
    if not os.path.exists(output_file):
        df = collect_data(score_file)
        write_output(df, output_file)
    else:
        df = pd.read_csv(output_file)

    # compare the scores to dimer
    avg_df = get_average(df)

    # visualize results
    visualize(avg_df)


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
    # convert df to csv file 
    df.to_csv(output_file, index=False)
    print(f'Full data csv file saved to {output_file}')

# to elucidate the reason why there are genes missing from score file between input.fa file
def compare_error(df, input_fa_file):
    
    pass


# to check average score given flanking dimer sequence
def get_average(df):
    # get the counts of both dimers
    donor_counts = df['donorDimer'].value_counts()
    acceptor_counts = df['acceptorDimer'].value_counts()

    # combine both series into a df
    full_df = pd.concat([donor_counts, acceptor_counts], axis=1).fillna(0).astype('int64')

    # get the average scores for each dimer
    full_df['donorAvgScore'] = 0.0
    full_df['acceptorAvgScore'] = 0.0

    for target_dimer in full_df.index:
        # filter out the target dimers, then get average score
        donor_filt = df[df['donorDimer'] == target_dimer]
        donor_avg = donor_filt['donorScore'].mean()
        acceptor_filt = df[df['acceptorDimer'] == target_dimer]
        acceptor_avg = acceptor_filt['acceptorScore'].mean()

        # save to corresponding index in df
        full_df.at[target_dimer, 'donorAvgScore'] = donor_avg
        full_df.at[target_dimer, 'acceptorAvgScore'] = acceptor_avg

    print(full_df)
    return full_df


def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext' format
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path


def visualize(df):
    # do the plot
    ax = sns.scatterplot(data=df, x='donorDimer', y='donorAvgScore')
    plt.xscale('log')

    # set plot title and axis labels
    plt.title('Correlation between Donor Dimer Frequency and Score')
    plt.xlabel('Dimer Frequency')
    plt.ylabel('Dimer Score')

    # display the plot
    plt.show()

    # save the plot
    fig_path = './outputs/FIGURES/acceptor_dimer_freq_vs_score.png'
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    ax.figure.savefig(handle_duplicate_names(fig_path))


if __name__ == '__main__':
    dbs = ['chess3', 'gencode_all', 'MANE', 'refseq']
    nums = [0] # CHANGEME

    for num in nums:
        run_aggregator(dbs[num])