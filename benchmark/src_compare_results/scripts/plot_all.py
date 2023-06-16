# run in /benchmark/src_compare_results/

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

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


def make_fig7(type):
    # defining paper figure for comparing SPLAM and SpliceAI in histogram across all species for donor and acceptor
    
    # data gathering
    species = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    sites = ['d_splam', 'd_spliceai', 'a_splam', 'a_spliceai']
    cols = [x+'_'+y for x in species for y in sites]
    all_species_df = pd.DataFrame()

    for idx, db in zip(range(0, len(cols), 4), species):
        avg_ofp = f'./output/aggregate/avg_data.{type}.{db}.csv'
        df = pd.read_csv(avg_ofp)

        t_df = pd.DataFrame(columns=cols[idx:idx+4])
        t_df.iloc[:,0] = df['d_score_splam']
        t_df.iloc[:,1] = df['d_score_spliceai']
        t_df.iloc[:,2] = df['a_score_splam']
        t_df.iloc[:,3] = df['a_score_spliceai']
        
        all_species_df = pd.concat([all_species_df, t_df], axis=1)

    print(all_species_df)


    # plot setup
    sns.set(font_scale=0.8)
    sns.set_palette('deep')
    f, axs = plt.subplots(2, 4, figsize=(16, 10), sharey=True)
    plt.subplots_adjust(hspace=None,wspace=None)
    f.suptitle('Comparison of SPLAM performance across non-human species', y=0.95)
    f.supylabel('Site')
    f.supxlabel('Score')
    plt.rcParams['font.family'] = 'Arial'

    # iterative plotting for each column of dataframe
    axs = axs.flatten(order='F')
    for i in range(0,16,2):
        data = all_species_df.iloc[:,i:i+2]
        ax = axs[i//2]

        sns.kdeplot(data=data, ax=ax, clip=(0.0, 1.0), fill=True)

        if i % 2 == 0 and i < 8:
            axs[i].set_title(species[i//2])

    axs[0].set_ylabel('Donor')
    axs[1].set_ylabel('Acceptor')

    plt.show()

    save_fig(f'./figures/fig7/all_scores_compare.{type}.png')

if __name__ == '__main__':

    type = 'noN'

    make_fig7(type)