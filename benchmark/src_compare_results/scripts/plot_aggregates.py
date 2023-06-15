import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
import ast
from progress.bar import Bar

def run_plotter(type, db):

    # obtain the aggregated and averaged data from all 5 versions of spliceai and write to file
    agg_ofp = f'./output/aggregate/agg_full_data.{type}.{db}.csv'
    avg_ofp = f'./output/aggregate/avg_data.{type}.{db}.csv'
    if not os.path.exists(avg_ofp):
        agg_df = aggregate_data(type, db, agg_ofp)
        df = get_averages(agg_df, avg_ofp)
    else:
        agg_df = pd.read_csv(agg_ofp)
        agg_df = revive_list(agg_df, ['d_score_spliceai', 'a_score_spliceai', 'spliceai_version'], ['float64', 'float64', 'int'])
        df = pd.read_csv(avg_ofp)
        df = revive_list(df, ['spliceai_version'], ['int'])

    # visualize the results
    id = [type, db]

    #make_fig3(df, id)
    #make_fig4(agg_df, id)
    make_fig5(agg_df, id)


def revive_list(df, col_names, dtypes):
    
    for name, dt in zip(col_names, dtypes):
        # revive the lists in column from string representations into np arrays
        df[name] = df[name].apply(lambda x: np.array(x.strip('[]').split(','), dtype=dt))

    return df


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


def save_fig(figpath):
    fig_path = handle_duplicate_names(figpath)
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=500)
    print(f'Saved figure to {fig_path}.')


def make_fig3(df, id):
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
    plt.title(f'Score Distributions between SPLAM and Averaged SpliceAI Models for {id[1]} Dataset')
    plt.xlabel('Site')
    plt.tight_layout()
    plt.show()

    save_fig(f'./figures/fig3/agg_scores_distribution.{id[0]}.{id[1]}.png')


def make_fig4(agg_df, id):
    # plotting SpliceAI scores as a comparative histogram for different versions of the model

    df = agg_df.iloc[:,10:]
    df = df.apply(pd.Series.explode).astype('float64').reset_index() # explodes the lists into more rows
    df['spliceai_version'] = df['spliceai_version'].astype('int')
    df.drop('index', axis=1, inplace=True)
    print(df)
    
    f, axes = plt.subplots(1,2, figsize=(12,6), sharey=True)
    sns.set(font_scale=0.8)
    f.suptitle(f'Score Distributions between SPLAM and Averaged SpliceAI Models for {id[1]} Dataset')
    axes[0].set_title('Donor')
    axes[1].set_title('Acceptor')
    #sns.kdeplot(data=df, x='d_score_spliceai', hue='spliceai_version', fill=True, alpha=0.05, ax=axes[0], palette='deep')
    sns.histplot(data=df, x='d_score_spliceai', hue='spliceai_version', kde=True, ax=axes[0], palette='deep')
    sns.move_legend(axes[0], 'upper left')
    #sns.kdeplot(data=df, x='a_score_spliceai', hue='spliceai_version', fill=True, alpha=0.05, ax=axes[1], palette='deep')
    sns.histplot(data=df, x='a_score_spliceai', hue='spliceai_version', kde=True, ax=axes[1], palette='deep')
    sns.move_legend(axes[1], 'upper left')
    plt.show()

    save_fig(f'./figures/fig4/spliceai_versions_distribution.{id[0]}.{id[1]}.png')


def make_fig5(agg_df, id):
    # plot SpliceAI scores across different models against intron length

    df = agg_df[['start', 'end', 'd_score_spliceai', 'a_score_spliceai', 'spliceai_version']]
    df = df.apply(pd.Series.explode).astype('float64').reset_index() # explodes the lists into more rows
    df['spliceai_version'] = df['spliceai_version'].astype('int')
    df['intron_length'] = df['end'] - df['start']
    df['intron_length'] = df['intron_length'].astype('int')
    df.drop(['index', 'start', 'end'], axis=1, inplace=True)
    print(df)

    sns.set(font_scale=0.8)
    plot_params = {'s': 3, 'alpha': 0.3}
    # f, axes = plt.subplots(1,2, figsize=(12,6), sharey=True)
    # f.suptitle(f'Intron Length vs. Score across SpliceAI Models for {id[1]} Dataset')
    # axes[0].set_title('Donor')
    # axes[1].set_title('Acceptor')
    # sns.scatterplot(data=df, x='intron_length', y='d_score_spliceai', hue='spliceai_version', ax=axes[0], palette='deep', **plot_params)
    # sns.scatterplot(data=df, x='intron_length', y='a_score_spliceai', hue='spliceai_version', ax=axes[1], palette='deep', **plot_params)
    # axes[0].set_xscale('log')
    # axes[1].set_xscale('log')

    sns.relplot(data=df, x='intron_length', y='a_score_spliceai', hue='spliceai_version', col='spliceai_version', palette='deep', **plot_params)
    plt.xscale('log')

    plt.show()

    save_fig(f'./figures/fig5/acceptor_len_vs_score_all_sai_vers.{id[0]}.{id[1]}.png')

if __name__ == '__main__':

    type = 'noN'
    databases = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    db_nums = [3]

    for num in db_nums:
        run_plotter(type, databases[num])