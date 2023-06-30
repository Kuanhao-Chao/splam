# aggregate the outputs from SPLAM to a single .csv file

import os
import csv
from progress.bar import Bar
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text
#from Bio import SeqIO
from pyfaidx import Fasta 
#from itertools import groupby


# use the score file as the input file, then check with the original GRCh28.p14 fasta file for the donor and acceptor dimers
# cross-reference with the INPUT.fa file to see why some of the entries get deleted -> write those to a new file?
# if you set first line to track_name=junctions, you can visualize better using igv

def run_aggregator(db_name):
    BASE_LEN = 800

    # collect the data to be aggregated from multiple sources
    score_file = f'../SPLAM/5_output/{BASE_LEN}bp/{db_name}/{db_name}.score.bed'
    output_file = f'../data/{db_name}.splam_data.csv'

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # make csv file if it does not exist, write to pd dataframe
    if not os.path.exists(output_file):
        df = collect_data(score_file, db_name)
        write_output(df, output_file)
    else:
        df = pd.read_csv(output_file)


    # CHANGEME select which figures to make
    make_fig1 = False
    make_fig2 = True

    if make_fig1:
        # compare the scores to dimer
        avg_df = get_average(df)

        # visualize results (fig1)
        da = ['donor', 'acceptor']
        log = [True, False]
        for args in [(avg_df, db_name, i, j) for i in da for j in log]:
            visualize_f1(*args)
    
    if make_fig2:
        # visualize results (fig2)
        visualize_f2(df, db_name)


def collect_data(score_file, db_name):
    # read score data from score file
    df = pd.read_csv(score_file, sep='\t', header=None, \
        names=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score', 'strand', 'donorScore', 'acceptorScore'])

    # add new columns for the dimers
    df['donorDimer'] = ''
    df['acceptorDimer'] = ''

    # ping the fasta file using pyfaidx to obtain the sequences as a dictionary
    genes = Fasta('../../SPLAM_python/extraction/primates/' + db_name + '_genomic.fa')
    print(f'Found {len(genes.keys())} unique genes in {db_name} dataset.')

    # iterate over every coordinate and extract corresponding dimer from fasta
    pbar = Bar('[Info] Getting dimers...', max=len(df))
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for index, row in df.iterrows():
        # obtain the start and end sequences from the fasta file
        chromosome = row['chrom'].strip('>')
        row['chrom'] = chromosome # fixing some naming error in the input file
        donor_start = int(row['chromStart(donor)'])
        acceptor_end = int(row['chromEnd(acceptor)'])

        # extract the dimers from the fasta file
        donor_dimer = str(genes[chromosome][donor_start:donor_start+2]).upper()
        acceptor_dimer = str(genes[chromosome][acceptor_end-2:acceptor_end]).upper()

        # convert to reverse complement if on reverse strand
        if (row['strand'] == '-'):
            reverse = ''.join([donor_dimer, acceptor_dimer])[::-1]
            r_complement = ''.join(complement[base] for base in reverse)
            donor_dimer = r_complement[:2]
            acceptor_dimer = r_complement[2:]

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

    # handle NaN float values
    full_df = full_df.fillna(0)

    print(full_df)
    return full_df

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path


def visualize_f1(df, db, d_a = 'donor', log_xscale = True):
    # do the plot
    sns.set_palette('deep')
    sns.set(font_scale=0.8)
    ax = sns.scatterplot(data=df, x=f'{d_a}Dimer', y=f'{d_a}AvgScore', s=15)
    if (log_xscale):
        plt.xscale('log')

    # add labels to the points
    texts = [ax.text(row[f'{d_a}Dimer'], row[f'{d_a}AvgScore'], str(i), size=5, ha='center', va='center') for i, row in df.iterrows()]
    adjust_text(texts, precision=0.5)

    # set plot title and axis labels
    plt.title(f'Correlation between {d_a.title()} Dimer Counts and Score for {db} Database')
    plt.xlabel('Dimer Counts')
    plt.ylabel('Dimer Score')

    # display the plot
    plt.show()

    # save the plot
    fig_path = handle_duplicate_names(f'../figures/fig1/{d_a}_dimer_freq_vs_score_{db}_log{log_xscale}.png')
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    ax.figure.savefig(fig_path, dpi=300)
    print(f'Saved figure to {fig_path}.')


# to elucidate the reason why there are genes missing from score file between input.fa file
def compare_error(df_score, input_fa_file, output_file):
    # extract rows from input fasta file into dataframe
    columns = ['chromosome', 'donorPos', 'acceptorPos', 'strand', 'score', 'seq']
    data = []
    linenums = []
    with open(input_fa_file, 'r') as fasta:
        linenum = 1
        for line in fasta:
            if line.startswith('>'):
                linenums.append(linenum)
                linenum += 1
                # get attributes from header
                header = line.split(';')
                chromosome = header[0][1:]
                donor_pos = int(header[1])
                acceptor_pos = int(header[2])
                strand = header[3]
                score = header[4]
                
            else:
                data.append((chromosome, donor_pos, acceptor_pos, strand, score, line))
    df_inp = pd.DataFrame(data, columns=columns)

    print(f'Input df length: {len(df_inp)}\nScore df length: {len(df_score)}')
    print(f'Lines with leading > {len(linenums)}')

# plot intron length vs. score
def visualize_f2(df, db): 
    # calculate the intron length for each entry and format into new df
    newdf = pd.DataFrame()
    newdf['intron_length'] = df['chromEnd(acceptor)'] - df['chromStart(donor)'] + 1
    newdf = pd.concat([newdf, newdf], axis=0)
    d_df = pd.DataFrame()
    d_df['score'] = df['donorScore']
    d_df['type'] = 'donor'
    a_df = pd.DataFrame()
    a_df['score'] = df['acceptorScore']
    a_df['type'] = 'acceptor'
    scoredf = pd.concat([d_df, a_df], axis=0)
    newdf = pd.concat([newdf, scoredf], axis=1)

    # get rows where intron length is <400
    # newdf = newdf[newdf['intron_length'] < 400]
    # print(len(newdf))

    # plt.figure(figsize=(12,8))
    # sns.scatterplot(x='intron_length', y='donorScore', data=df, color='blue', label='Donor Score', **plot_params)
    # sns.scatterplot(x='intron_length', y='acceptorScore', data=df, color='orange', label='Acceptor Score', **plot_params)
    
    # plot the scatterplot with marginal distributions
    sns.set(font_scale=0.8)
    plot_params = {'s': 3, 'alpha': 0.25}
    sns.jointplot(data=newdf, x='intron_length', y='score', hue='type', height=8, ratio=4, **plot_params)
    
    # set up the plot
    plt.xlabel('Intron Length')
    plt.ylabel('Score')
    plt.title(f'Intron Lengths vs. Score for {db} Database')
    plt.xscale('log')
    plt.legend(fontsize=8)

    # save the plot
    fig_path = handle_duplicate_names(f'../figures/fig2/intron_length_vs_scores_{db}.png')
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=500)
    print(f'Saved figure to {fig_path}.')

    # show the plot
    plt.show()


if __name__ == '__main__':
    #dbs = ['chess3', 'gencode_all', 'MANE', 'refseq']
    dbs = ['Mmul_10', 'NHGRI_mPanTro3', 'GRCm39', 'TAIR10']
    nums = [0,1,2,3] # CHANGEME

    for num in nums:
        run_aggregator(dbs[num])