import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
import ast
from progress.bar import Bar

def run_plotter(type, db):

    # obtain the aggregated data from all 5 versions of spliceai and write to file
    agg_ofp = f'./output/aggregate/agg_full_data.{type}.{db}.csv'
    if not os.path.exists(agg_ofp):
        df = aggregate_data(type, db, agg_ofp)
    else:
        #df = pd.read_csv(ofp, converters={'d_score_spliceai': ast.literal_eval, 'a_score_spliceai': ast.literal_eval})
        df = pd.read_csv(agg_ofp)
        df = parse_list(df)
    
    # calculate the averages of all score columns
    avg_ofp = f'./output/aggregate/avg_data.{type}.{db}.csv'
    df = get_averages(df, avg_ofp)

    # visualize the results
    fig3_path = f'./figures/fig3/agg_scores_distribution.{type}.{db}.png'

    #make_fig3(df, fig3_path, db)


def aggregate_data(type, db, ofp):

    # define the df which will collect all scores together
    aggregate_df = pd.DataFrame(columns=['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                                            'd_score_splam', 'a_score_splam', 'd_score_spliceai', 'a_score_spliceai', 'spliceai_version'])

    # aggregate all versions of spliceai
    for version_num in range(1, 6):
    
        # read the input file
        ifp = f'./output/comparison/full_data.v{version_num}.{type}.{db}.csv'
        full_df = pd.read_csv(ifp)

        # add column indicating which version of the model data is from
        full_df['spliceai_version'] = version_num

        # concat onto aggregate df
        aggregate_df = pd.concat([aggregate_df, full_df], axis=0)

        print(aggregate_df)

    # group the dataframe by common columns and aggregate the spliceai scores into a list
    # NOTE: since there are duplicate sequences with different spliceai scores for a given file, there can be more than 5 scores in a merged list
    merged_df = aggregate_df.groupby(['seqid', 'start', 'end', 'name', 'expected_score', 'strand', 'd_dimer', 'a_dimer', 
                                      'd_score_splam', 'a_score_splam']).agg(lambda x: x.tolist()).reset_index()

    print(f'Preview merged_df:\n{merged_df}')
    
    # save to csv
    os.makedirs(os.path.dirname(ofp), exist_ok=True)
    merged_df.to_csv(ofp, index=False)
    print(f'Saved aggregate csv file to {ofp}.')

    return merged_df


def parse_list(df):

    # revive the lists from string representations into np arrays
    df['d_score_spliceai'] = df['d_score_spliceai'].apply(lambda x: np.array(x.strip('[]').split(','), dtype='float64'))
    df['a_score_spliceai'] = df['a_score_spliceai'].apply(lambda x: np.array(x.strip('[]').split(','), dtype='float64'))
    df['spliceai_version'] = df['spliceai_version'].apply(lambda x: np.array(x.strip('[]').split(','), dtype='float64'))

    return df


def get_averages(df, ofp):

    # take the averages of the spliceai score columns
    df['d_score_spliceai'] = df['d_score_spliceai'].apply(np.average)
    df['a_score_spliceai'] = df['a_score_spliceai'].apply(np.average)

    print(f'Preview df with averages:\n{df}')
    
    # save to csv
    os.makedirs(os.path.dirname(ofp), exist_ok=True)
    df.to_csv(ofp, index=False)
    print(f'Saved average csv file to {ofp}.')

    return df


def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext'
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path


def make_fig3(df, fig_path, db):
    # plotting the comparative score distributions between SPLAM and SpliceAI as a score distribution, for averaged SpliceAI across 5 models

    d_melt = pd.melt(df, value_vars=['d_score_spliceai', 'd_score_splam'], var_name='Method', value_name='Score').replace('d_score_spliceai', 'SpliceAI').replace('d_score_splam', 'SPLAM')
    a_melt = pd.melt(df, value_vars=['a_score_spliceai', 'a_score_splam'], var_name='Method', value_name='Score').replace('a_score_spliceai', 'SpliceAI').replace('a_score_splam', 'SPLAM')
    d_melt['Type'] = 'Donor'
    a_melt['Type'] = 'Acceptor'
    reshaped_df = pd.concat([d_melt, a_melt], axis=0)
    print(reshaped_df)

    sns.set(font_scale=0.8)
    plt.figure(figsize=(6,8))
    sns.violinplot(data=reshaped_df, x='Type', y='Score', hue='Method', hue_order=['SPLAM', 'SpliceAI'],
                    split=True, inner='quart', linewidth=0.9)
    plt.title(f'Score Distributions between SPLAM and Averaged SpliceAI Models for {db} Dataset')
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

    type = 'noN'
    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    db_nums = [0,1,2,3]

    for num in db_nums:
        run_plotter(type, databases[num])