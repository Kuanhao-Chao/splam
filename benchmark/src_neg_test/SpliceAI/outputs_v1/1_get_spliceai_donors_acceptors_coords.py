import os
import sys
import pandas as pd

def get_chrom_size(path):
    chrs = {}
    with open(path, 'r') as file:
        # read the file line by line
        for line in file:  
            if line.startswith('#'):
                continue
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    
    return chrs

def main(argv):
    db = argv[0]
    bl = argv[1]

    chrs = get_chrom_size(f'../../SPLAM_python/extraction/primates/{db}_assembly_report.txt')

    df = pd.read_csv(f'../SPLAM/2_output/{bl}bp/{db}/d_a.bed', delimiter='\t', header=None, usecols=range(6),
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])
    
    # add the flanking sequence to either end of the transcript
    df['start'] -= 5200
    df['end'] += 5200

    # perform validation to shorten sequence to ends of the chromosome if over length
    cuts_file = './1_output/{db}_cuts.csv'
    os.makedirs(os.path.dirname(cuts_file), exist_ok=True)
    with open(cuts_file, 'w') as f:
        for index, row in df.iterrows():
            chrom = row['seqid']
            if (row['start'] < 0):
                df.at[index, 'start'] = 0
                # df.drop(index, axis=0, inplace=True)
                # continue
            if (row['end'] > chrs[chrom]):
                df.at[index, 'end'] = chrs[chrom]
                # df.drop(index, axis=0, inplace=True)
                # continue
        

    # obtain a random sample (reproducible) for further analysis
    n = 50000
    df = df.sample(n, random_state=3217)
    
    save_path = f'./1_output/{db}_spliceai.bed'
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    df.to_csv(save_path, sep='\t', header=None, index=0)

if __name__ == "__main__":

    if os.getcwd() != 'SpliceAI':
        os.chdir('/home/smao10/SPLAM/benchmark/src_neg_test/SpliceAI/')

    main(sys.argv[1:])