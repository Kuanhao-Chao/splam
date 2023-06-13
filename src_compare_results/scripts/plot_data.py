import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
from progress.bar import Bar

def run_plotter(argv, db):

    # define identifiers for this run
    TYPE = 'noN'
    SPLICEAI_VERSION = 2
    print(f'Parsing for type {TYPE}, SpliceAI version {SPLICEAI_VERSION}, database {db}')

    # input filepaths
    dirpath = f'./spliceai_result_{SPLICEAI_VERSION}/{db}/'
    d_score_file = f'{dirpath}spliceai_all_seq.score.d.{TYPE}.{db}.tsv'
    a_score_file = f'{dirpath}spliceai_all_seq.score.a.{TYPE}.{db}.tsv'
    name_file = f'{dirpath}spliceai_all_seq.name.{TYPE}.{db}.tsv'
    full_data_score_file = f'./data/{db}.full_data.csv'
    # n_score_file = f'{dirpath}spliceai_all_seq.score.n.{TYPE}.{db}.tsv'

    # output filepaths
    id = f'.v{SPLICEAI_VERSION}.{TYPE}.{db}'
    pack = [SPLICEAI_VERSION, TYPE, db]
    csvpath = f'./output/comparison/full_data{id}.csv'
    fig1_path = f'./figures/fig1/scores_distribution{id}.png'
    fig2_path = f'./figures/fig2/scores_distribution_comparison{id}.png'
    
    # collect the full data into a dataframe
    index_compare(full_data_score_file, name_file)


    # if not os.path.exists(csvpath):
    #     data_df, name_df, indices = index_compare(full_data_score_file, name_file)
    #     score_df, donor_df, acceptor_df = collect_data(d_score_file, a_score_file)
    #     df = combine_and_write(score_df, name_df, data_df, indices, csvpath)
    # else:
    #     df = pd.read_csv(csvpath)

    # make figures 
    #make_fig1(df, fig1_path)
    #make_fig2(df, fig2_path, pack)

def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path

'''Compare the indices between the SPLAM data file and SpliceAI output and get only the matches'''
def index_compare(full_data_file, name_file, score_df):

    # read files and create dfs (df1 = SPLAM file; df2 = SpliceAI file w. scores)
    splam_full_df = pd.read_csv(full_data_file)
    spliceai_name_df = pd.read_csv(name_file, delimiter=' ', names=['id', 'start', 'end', 'strand'])

    df1 = splam_full_df
    df2 = pd.concat([spliceai_name_df, score_df], axis=1)

    # starts1 = df1.iloc[:,1:3].apply(tuple, axis=1)
    # print(starts1)
    # set1 = set(sorted(starts1))
    # starts2 = df2.iloc[:,1:3].apply(tuple, axis=1)
    # print(starts2)
    # set2 = set(sorted(starts2))
    # combined = set1.union(set2)
    # print(len(set1), len(set2), len(combined))
    # print(df1.shape, df2.shape)

    # get the identifiers (seqid, start, end) as a list of tuples to search for matches
    print('Filtering full data...')
    df1_rows = df1.iloc[:,:3].apply(tuple,axis=1)
    df2_rows = df2.iloc[:,:3].apply(tuple,axis=1)
    print(f'\t\tRows\t\tUniques\nSPLAM:\t\t{len(df1_rows)}\t{len(df1_rows.drop_duplicates())}\nSpliceAI:\t{len(df2_rows)}\t{len(df2_rows.drop_duplicates())}')
    
    filtered_df = df1.merge(df2, how='inner', left_on=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'strand'], right_on=['id', 'start', 'end', 'strand'], sort=True).drop_duplicates()
    indices = df2_rows.drop_duplicates().index 

    print(f'Shape of merged: {filtered_df.shape}, size of indices: {indices.shape}')
    filtered_df.to_csv('./output/test.csv')
    # get the uniques from each row
    # df1_rows = df1_rows.drop_duplicates()
    # df2_rows = df2_rows.drop_duplicates()
  
    # get the matches (filtered_df = rows of df1 which match; indices = indices of SpliceAI files which are unique and are used here)
    # filtered_df = df1[df1_rows.isin(df2_rows)].sort_values(by=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)'])
    # indices = df2_rows.index

    

    print(f'Preview filtered SPLAM values:\n{filtered_df}')
    print('Finished filtering')

    return filtered_df, df2, indices

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

'''Combine the full data with SpliceAI scores and save to csv file'''
def combine_and_write(score_df, name_df, full_df, indices, csvpath):
    
    # combine the spliceai scores, then get only the unique rows corresponding to full data
    comb_df = pd.concat([name_df, score_df], axis=1)
    comb_df.columns = ['name', 'start', 'end', 'strand', 'donor_score', 'acceptor_score']
    comb_df = comb_df.iloc[indices].sort_values(by=['name', 'start', 'end']).reset_index(drop=True)
    full_df = full_df.reset_index(drop=True)

    print(f'Length of\nspliceai: {len(comb_df)}\nfull splam: {len(full_df)}')
    print(comb_df)
    print(full_df)

    # maybe make a method to validate the matches between the two dataframes? (chrom, start, end, strand?)
    

    # rearrange, rename, and combine the final df columns
    new_headers = ['seqid', 'start', 'end', 'name', 'exp_score', 'strand', \
                   'd_dimer', 'a_dimer', 'd_score_splam', 'a_score_splam']
    full_df = full_df.reindex(columns=['chrom', 'chromStart(donor)', 'chromEnd(acceptor)', 'name', 'score', \
                               'strand', 'donorDimer', 'acceptorDimer', 'donorScore', 'acceptorScore'])
    full_df.columns = new_headers
    full_df['d_score_spliceai'] = comb_df['donor_score']
    full_df['a_score_spliceai'] = comb_df['acceptor_score']

    print(f'Preview final saved df:\n{full_df}')

    # save to file
    os.makedirs(os.path.dirname(csvpath), exist_ok=True)
    full_df.to_csv(csvpath, index=False)
    print(f'Saved csv file to {csvpath}.')
    print('Finished combining into csv')

    return full_df



def check_match(id1, id2):
    pass


def make_fig1(df, fig_path):

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

    d_melt = pd.melt(df, value_vars=['d_score_spliceai', 'd_score_splam'], var_name='Method', value_name='Score').replace('d_score_spliceai', 'SpliceAI').replace('d_score_splam', 'SPLAM')
    a_melt = pd.melt(df, value_vars=['a_score_spliceai', 'a_score_splam'], var_name='Method', value_name='Score').replace('a_score_spliceai', 'SpliceAI').replace('a_score_splam', 'SPLAM')
    d_melt['Type'] = 'Donor'
    a_melt['Type'] = 'Acceptor'
    reshaped_df = pd.concat([d_melt, a_melt], axis=0)
    print(reshaped_df)

    sns.set(font_scale=0.8)
    plt.figure(figsize=(6,8))
    #sns.violinplot(data=reshaped_df, x='Score', y='Type', palette='viridis')
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
    nums = [3]

    for num in nums:
        run_plotter(sys.argv[1:], databases[num])