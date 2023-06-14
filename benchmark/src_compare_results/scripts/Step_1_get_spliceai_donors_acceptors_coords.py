import os
import sys
import pandas as pd

def get_chrom_size(path):
    chrs = {}
    with open(path, 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    
    return chrs

def main(argv):

    chrs = get_chrom_size(f'../Dataset/{argv[1]}_assembly_report.txt')

    df = pd.read_csv(f'../{argv[0]}/{argv[1]}_parsed.bed', delimiter='\t', header=None, \
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])
    
    # add the flanking sequence to either end of the transcript
    df['start'] -= 5200
    df['end'] += 5200

    # perform validation to shorten sequence to ends of the chromosome if over
    for index, row in df.iterrows():
        chrom = row['seqid']
        if (row['start'] < 0):
            df.at[index, 'start'] = 0
        if (row['end'] > chrs[chrom]):
            df.at[index, 'end'] = chrs[chrom]
    
    save_path = f'../output/{argv[0]}/{argv[1]}_spliceai.bed'
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    df.to_csv(save_path, sep='\t', header=None, index=0)

if __name__ == "__main__":
    main(sys.argv[1:])