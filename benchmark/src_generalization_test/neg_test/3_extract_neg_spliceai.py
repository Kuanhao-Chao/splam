import os
import pandas as pd
from pyfaidx import Fasta
from progress.bar import Bar 

'''Get the length of chromosomes'''
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

'''Extract the negative SpliceAI-formatted sequences'''
def extract(df, seqs, chrs):

    # add the flanking sequence to either end of the transcript
    df['flank_start'] = df['start'] - 5200
    df['flank_end'] = df['end'] + 5200

    pbar = Bar('Extracting negative sequences...', max=len(df))
    for idx, row in df.iterrows():

        # get sequence and size
        chrom = row['seqid']
        seq = seqs[chrom]
        chrom_size = chrs[chrom]

        # sanity checks pt. 1
        prestart = row['start']
        preend = row['end']
        prelength = preend - prestart
        extract = seq[prestart:preend]
        if row['strand'] == '+':
            assert(extract[:2] == 'GT')
            assert(extract[-2:] == 'AG')
        elif row['strand'] == '-':
            assert(extract[:2] == 'CT')
            assert(extract[-2:] == 'AC')

        # init pads
        left_pad = 0
        right_pad = 0
        start = row['flank_start']
        end = row['flank_end']
        strand = row['strand']

        # check for any sequences out of bounds
        if start < 0:
            left_pad = 0 - start
            df.at[idx, 'flank_start'] = 0 # modifies the df in-place
            start = 0 # modifies local variable

        if end > chrom_size:
            right_pad = end - chrom_size
            df.at[idx, 'flank_end'] = chrom_size # modifies the df in-place
            end = chrom_size # modifies local variable

        # extract the noN sequence
        if strand == '+':
            noN_sequence = (left_pad*'N') + str(seq[start:end]) + (right_pad*'N')
        elif strand == '-':
            noN_sequence = (right_pad*'N') + str(seq[start:end].reverse.complement) + (left_pad*'N')

        # modify to get the N sequence
        N_sequence = (5000*'N') + noN_sequence[5000:-5000] + (5000*'N')

        # add sequences to df
        df.at[idx, 'noN_seq'] = noN_sequence
        df.at[idx, 'N_seq'] = N_sequence

        # sanity checks pt. 2
        assert(noN_sequence[5200:5202] == 'GT')
        assert(noN_sequence[-5202:-5200] == 'AG')
        assert(N_sequence[5200:5202] == 'GT')
        assert(N_sequence[-5202:-5200] == 'AG')    
        assert(len(noN_sequence) == prelength + 10400)
        assert(len(N_sequence) == prelength + 10400)
        assert(len(noN_sequence) == len(N_sequence))

        pbar.next()
    pbar.finish()

    # obtain a random sample (reproducible) for further analysis
    n = 25000
    df = df.sample(n, random_state=3217)

    return df

def write_results(df, fw_coords, fw_noN_seq, fw_N_seq):

    print('Writing results...')

    # write coordinate BED file - these coords refer to the full sequence, not but including the flank padding
    bed_columns = ['seqid', 'flank_start', 'flank_end', 'name', 'score', 'strand']
    df[bed_columns].to_csv(fw_coords, sep='\t', header=None, index=0) 

    for idx, row in df.iterrows():
        # write header to both FASTA files
        # NOTE: this is not FASTA formatted, but makes it more SpliceAI compatible (0-indexed start, 1-indexed end)
        header = f'{row["seqid"]};{row["start"]};{row["end"]};{row["strand"]}\n'
        fw_noN_seq.write(header)
        fw_N_seq.write(header)

        # write noN sequence to FASTA file
        fw_noN_seq.write(f'{row["noN_seq"]}\n')

        # write N sequence to FASTA file
        fw_N_seq.write(f'{row["N_seq"]}\n')


def main(db):

    print(f'Parsing for {db} dataset:')

    # input files
    bed_file = f'./2_output/{db}/d_a.bed'
    fasta_file = f'../data/{db}_genomic.fa'
    assembly_file = f'../data/{db}_assembly_report.txt'
    
    # output files
    output_dir = f'./3_output/{db}/'
    os.makedirs(output_dir, exist_ok=True)
    fw_coords = open(f'{output_dir}coords.bed', 'w')
    fw_noN_seq = open(f'{output_dir}seq_noN.fa', 'w')
    fw_N_seq = open(f'{output_dir}seq_N.fa', 'w')

    # process inputs
    df = pd.read_csv(bed_file, delimiter='\t', header=None, usecols=range(6),
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])
    seqs = Fasta(fasta_file, sequence_always_upper=True)
    chrs = get_chrom_size(assembly_file)

    # extract the sequences
    df = extract(df, seqs, chrs)

    # write results to files
    write_results(df, fw_coords, fw_noN_seq, fw_N_seq)
    print(f'Result in {output_dir}')

    fw_coords.close()
    fw_noN_seq.close()
    fw_N_seq.close()

if __name__ == "__main__":

    if os.getcwd() != 'neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_generalization_test/neg_test/')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [0,1,2,3] #CHANGEME

    for idx in idxs:
        main(datasets[idx])