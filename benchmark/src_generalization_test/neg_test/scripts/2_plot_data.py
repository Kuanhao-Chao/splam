import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
from progress.bar import Bar
import itertools

def run_plotter(TYPE, SPLICEAI_VERSION, db):

    # define identifiers for this run
    print('*'*170)
    print(f'Parsing for type {TYPE}, SpliceAI version {SPLICEAI_VERSION}, database {db}')

    # input filepaths
    dirpath = f'../SpliceAI/4_output/spliceai_result_{SPLICEAI_VERSION}/{db}/'
    d_score_file = f'{dirpath}spliceai_all_seq.score.d.{TYPE}.{db}.tsv'
    a_score_file = f'{dirpath}spliceai_all_seq.score.a.{TYPE}.{db}.tsv'
    name_file = f'{dirpath}spliceai_all_seq.name.{TYPE}.{db}.tsv'
    splam_data = f'../data/{db}.splam_data.csv'

    # output filepaths
    id = f'.v{SPLICEAI_VERSION}.{TYPE}.{db}'
    pack = [SPLICEAI_VERSION, TYPE, db]
    csvpath = f'../data/comparison/full_data{id}.csv'
    fig1_path = f'..figures/fig1/scores_distribution{id}.png'
    fig2_path = f'../figures/fig2/filt_scores_distribution_comparison{id}.png'
    
    # collect the full data into a dataframe
    if not os.path.exists(csvpath):
        score_df, donor_df, acceptor_df = collect_data(d_score_file, a_score_file)
        full_df = index_compare(splam_data, name_file, score_df)
        df = write_df(full_df, csvpath)
    else:
        df = pd.read_csv(csvpath)

    # make figures 
    #make_fig1(df, fig1_path)
    make_fig2(df, fig2_path, pack)

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path

'''Parse the SpliceAI donor and acceptor files into a single dataframe'''
def collect_data(d_file, a_file):

    # read donor and acceptor spliceai files into list
    print('Reading donor score file...')
    with open(d_file, 'r') as d_read:
        d_lines = d_read.read().splitlines()
    print('Reading acceptor score file...')
    with open(a_file, 'r') as a_read:
        a_lines = a_read.read().splitlines()

    # collect 400 scores on donor side, then donor position score for each intron
    pbar = Bar('Collecting donor scores...', max=len(d_lines))
    donors = []
    donor_scores = []
    for line in d_lines:
        values = line.strip().split()
        donor_row = values[:400]
        donors.append(donor_row)
        donor_scores.append(donor_row[199])
        pbar.next()
    pbar.finish()

    # collect 400 scores on acceptor side, then acceptor position score for each intron
    pbar = Bar('Collecting acceptor scores...', max=len(a_lines))
    acceptors = []
    acceptor_scores = []
    for line in a_lines:
        values = line.split()
        acceptor_row = values[-400:]
        acceptors.append(acceptor_row)
        try:
            acceptor_scores.append(acceptor_row[200])
        except: 
            acceptor_scores.append(acceptor_row[199]) # for rare cutoffs in chrom length
        pbar.next()
    pbar.finish()

    # compile scores into a single dataframe
    score_df = pd.DataFrame()
    score_df['donor_score'] = pd.Series(donor_scores, dtype='float64')
    score_df['acceptor_score'] = pd.Series(acceptor_scores, dtype='float64')
    donor_df = pd.DataFrame(donors, dtype='float64')
    acceptor_df = pd.DataFrame(acceptors, dtype='float64')

    # print results of score 
    print(f'Preview SpliceAI scores:\n{score_df}')
    print(f'Averages:\n{score_df.mean()}')
    print('Finished parsing')

    return score_df, donor_df, acceptor_df

'''Compare the indices between the SPLAM data file and SpliceAI output and get only the matches'''
def index_compare(full_data_file, name_file, score_df):

    # read files and create dfs (df1 = SPLAM file; df2 = SpliceAI file w. scores)
    splam_full_df = pd.read_csv(full_data_file)
    spliceai_name_df = pd.read_csv(name_file, delimiter=' ')

    df1 = splam_full_df
    df2 = pd.concat([spliceai_name_df, score_df], axis=1)
    df2.columns = ['id', 'start', 'end', 'strand', 'd_score_spliceai', 'a_score_spliceai']

    # get the identifiers (seqid, start, end) as a list of tuples -> show uniques for each 
    print('Filtering full data...')
    df1_rows = df1.iloc[:,:3].apply(tuple,axis=1)
    df2_rows = df2.iloc[:,:3].apply(tuple,axis=1)
    print(f'\t\tRows\tUniques\nSPLAM:\t\t{len(df1_rows)}\t{len(df1_rows.drop_duplicates())}\nSpliceAI:\t{len(df2_rows)}\t{len(df2_rows.drop_duplicates())}')
    
    # merge the two dataframes by ID, sort result, and drop duplicates 
    # NOTE: this will still keep some duplicate results as SpliceAI may give different scores to same sequence 
    filtered_df = df1.merge(df2, how='inner', left_on=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'strand'], 
                            right_on=['id', 'start', 'end', 'strand'], sort=True).drop_duplicates()

    print(f'Merged rows: {len(filtered_df)}')
    #filtered_df.to_csv('./output/test.csv')

    print(f'Preview filtered SPLAM values:\n{filtered_df}')
    print('Finished filtering')

    return filtered_df

'''Combine the full data with SpliceAI scores and save to csv file'''
def write_df(full_df, csvpath):

    # drop unneccesary columns and fix indices
    full_df.drop(['id', 'start', 'end'], inplace=True, axis=1)
    full_df.reset_index(drop=True, inplace=True)

    # reorder and rename the final df columns
    full_df = full_df.reindex(columns=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score', 'strand', 'donorDimer', 'acceptorDimer', 
                                       'donorScore', 'acceptorScore', 'd_score_spliceai', 'a_score_spliceai'])
    full_df.columns = ['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                       'd_score_splam', 'a_score_splam', 'd_score_spliceai', 'a_score_spliceai']

    print(f'Preview final saved df:\n{full_df}')

    # save to file
    os.makedirs(os.path.dirname(csvpath), exist_ok=True)
    full_df.to_csv(csvpath, index=False)
    print(f'Saved csv file to {csvpath}.')
    print('Finished writing')

    return full_df

def make_fig1(df, fig_path):
    # plotting the score distribution for SpliceAI output as violin plot

    reshaped_df1 = pd.DataFrame(columns=['Type', 'Score'])
    reshaped_df1['Score'] = df['d_score_spliceai']
    reshaped_df1['Type'] = 'donor'
    reshaped_df2 = pd.DataFrame(columns=['Type', 'Score'])
    reshaped_df2['Score'] = df['a_score_spliceai']
    reshaped_df2['Type'] = 'acceptor'
    reshaped_df = pd.concat([reshaped_df1, reshaped_df2], axis=0)

    sns.set(font_scale=0.8)
    sns.violinplot(data=reshaped_df, x='Score', y='Type', palette='viridis')
    plt.tight_layout()
    plt.show()

    # save the plot
    fig_path = handle_duplicate_names(fig_path)
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=500)
    print(f'Saved figure to {fig_path}.')

def make_fig2(df, fig_path, id):
    # plotting the comparative score distributions between SpliceAI and SPLAM as a split violin plot

    d_melt = pd.melt(df, value_vars=['d_score_spliceai', 'd_score_splam'], var_name='Method', value_name='Score').replace('d_score_spliceai', 'SpliceAI').replace('d_score_splam', 'SPLAM')
    a_melt = pd.melt(df, value_vars=['a_score_spliceai', 'a_score_splam'], var_name='Method', value_name='Score').replace('a_score_spliceai', 'SpliceAI').replace('a_score_splam', 'SPLAM')
    d_melt['Type'] = 'Donor'
    a_melt['Type'] = 'Acceptor'
    reshaped_df = pd.concat([d_melt, a_melt], axis=0)
    print(reshaped_df)

    # filter out super low values
    reshaped_df = reshaped_df[reshaped_df['Score'] > 0.0005]
    print(reshaped_df['Method'].value_counts())

    sns.set(font_scale=0.8)
    plt.figure(figsize=(6,8))
    sns.violinplot(data=reshaped_df, x='Type', y='Score', hue='Method', hue_order=['SPLAM', 'SpliceAI'],
                    split=True, inner='quart', linewidth=0.9)
    plt.title(f'Score Distributions between SPLAM and SpliceAI for {id[2]} Dataset')
    plt.xlabel('Site')
    plt.tight_layout()
    plt.show()

    # save the plot
    fig_path = handle_duplicate_names(fig_path)
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=500)
    print(f'Saved figure to {fig_path}.')
    pass

if __name__ == '__main__':
    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    type = 'N'
    versions = [1,2,3,4,5]
    db_nums = [0,1,2,3]

    for ver, num in itertools.product(versions, db_nums):
        run_plotter(type, ver, databases[num])