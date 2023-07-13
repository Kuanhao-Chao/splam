import pandas as pd
import os


def main(db):
    old_file = f'./1_output/old_method/{db}_introns.bed'
    new_file = f'./1_output/{db}_introns.bed'

    old_df = pd.read_csv(old_file, delimiter='\t', header=None, usecols=range(6),
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])
    new_df = pd.read_csv(new_file, delimiter='\t', header=None, usecols=range(6),
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])

    print(len(old_df), len(new_df))

    cols = ['seqid', 'start', 'end', 'strand']
    old_df.drop_duplicates(subset=cols, inplace=True)
    new_df.drop_duplicates(subset=cols, inplace=True)


    print(len(old_df), len(new_df))

    df_diff = pd.merge(old_df, new_df, on=cols, indicator=True, how='outer')
    old_only = df_diff.query("_merge == 'left_only'").drop('_merge', axis=1)
    new_only = df_diff.query("_merge == 'right_only'").drop('_merge', axis=1)

    print(len(old_only), len(new_only))

    os.makedirs('./1_output/compare/', exist_ok=True)
    old_only.to_csv(f'./1_output/compare/{db}_old_only.bed', sep='\t', header=None, index=False)
    new_only.to_csv(f'./1_output/compare/{db}_new_only.bed', sep='\t', header=None, index=False)

if __name__ == '__main__':
    
    for db in ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']:
        main(db)


        ######################################################################################################
        # Verdict: gffutils method is less accurate (picks up less introns) than my previous method
        # (verified using igv) 
        ######################################################################################################