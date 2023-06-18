# run in /benchmark/src_compare_results/

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import scipy

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
    f, axs = plt.subplots(2, 4, figsize=(16, 10), sharey=True, sharex=True)
    plt.subplots_adjust(hspace=None,wspace=None)
    #f.suptitle('Comparison of SPLAM performance across non-human species', y=0.95, size=12)
    f.supylabel('Density', size=17, weight=400)
    f.supxlabel('Score', size=17, weight=400)
    #plt.rcParams['font.family'] = 'Arial'

    # calculate quartiles
    quartiles = all_species_df.quantile([0.25, 0.50, 0.75], axis=0)
    print(quartiles)
    cols = ['b','b','b', 'r','r','r']

    # iterative plotting for each column of dataframe
    axs = axs.flatten(order='F')
    for i in range(0,16,2):
        data = all_species_df.iloc[:,i:i+2]
        ax = axs[i//2]

        sns.kdeplot(data=data, ax=ax, clip=(0.0, 1.0), fill=True, alpha=0.5).legend_.remove()
        ax.tick_params(axis='x', labelsize=10)

        # plot quarts
        quarts = quartiles.iloc[:,i:i+2].values.flatten(order='F')
        for xc, c in zip(quarts, cols):
            ax.axvline(x=xc, c=c, alpha=0.5)
    
        if i % 2 == 0 and i < 8:
            # x-axis species label
            axs[i].set_title(species[i//2], size=14, weight=500)

    # y-axis site label
    axs[0].set_ylabel('Donor Site', size=14, weight=500)
    axs[1].set_ylabel('Acceptor Site', size=14, weight=500)
    axs[0].tick_params(axis='y', labelsize=12, width=12)
    axs[1].tick_params(axis='y', labelsize=12, width=12)

    

    # move the legend out of the plot
    plt.legend(['SpliceAI-10k-noN', 'SPLAM'], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    
    # display and save
    plt.tight_layout()
    plt.show()

    save_fig(f'./figures/fig7/all_scores_compare.{type}.png')

if __name__ == '__main__':

    type = 'noN'

    make_fig7(type)