import os
from progress.bar import Bar
from pyfaidx import Fasta
import pandas as pd

SEQ_LENGTH = "800"
QUARTER_SEQ_LEN = int(SEQ_LENGTH) // 4
 
'''Obtain chromosome size'''
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


'''Extract the positive set sequences'''
def extract(df, seqs, chrs):
    
    pbar = Bar('Extracting positive sequences (Splam)...', max=len(df))
    for idx, row in df.iterrows():

        # initialize vars
        chrom = row['seqid']
        start = row['start']
        end = row['end']
        strand = row['strand']
        seq = seqs[chrom]
        chrom_size = chrs[chrom]
        splice_junc_len = end - start

        # filter out erroneous introns (can't have donor and acceptor sites overlap)
        if splice_junc_len < 4:
            print(f'Impossibly short intron at {chrom}:{start}-{end}({strand}). Skipping.')
            df.drop(idx, axis=0, inplace=True)
            continue

        # generate the 200nt flanking sequence each side (decreased if fully overlaps)
        flanking_size = QUARTER_SEQ_LEN
        center_pad = 0
        if splice_junc_len < QUARTER_SEQ_LEN:
            flanking_size = splice_junc_len
            center_pad = QUARTER_SEQ_LEN - flanking_size

        first_s = start - QUARTER_SEQ_LEN
        first_e = start + flanking_size
        second_s = end - flanking_size
        second_e = end + QUARTER_SEQ_LEN
        
        # create padding for boundary cases
        left_pad = 0
        right_pad = 0
        if first_s < 0:
            left_pad = 0 - first_s
            first_s = 0
        if second_e > chrom_size:
            right_pad = second_e - chrom_size
            second_e = chrom_size

        # sanity checks
        try: 
            assert(first_e > first_s)
            assert(second_e > second_s)
            assert(second_e > first_s)
        except:
            print(f'Mislabeled intron at {chrom}:{start}-{end}({strand}). Skipping.')
            df.drop(idx, axis=0, inplace=True)
            continue
        if strand == '+':
            df.at[idx, 'donor_dimer'] = str(seq[start:start+2])
            df.at[idx, 'acceptor_dimer'] = str(seq[end-2:end])
        elif strand == '-':
            df.at[idx, 'donor_dimer'] = str(seq[end-2:end].reverse.complement)
            df.at[idx, 'acceptor_dimer'] = str(seq[start:start+2].reverse.complement)
        else: 
            print('No strand information. Skipping...')
            df.drop(idx, axis=0, inplace=True)
            continue
        try:
            assert('N' not in df.at[idx, 'donor_dimer'])
            assert('N' not in df.at[idx, 'acceptor_dimer'])
        except: 
            print(f'Dimer contains N at {chrom}:{start}-{end}({strand}). Skipping.')
            df.drop(idx, axis=0, inplace=True)
            continue
        
        # get the donor and acceptor sequences
        if strand == '+':  
            donor_seq = seq[first_s:first_e]
            acceptor_seq = seq[second_s:second_e]
            
            # add padding
            donor_seq = (left_pad*'N') + str(donor_seq) + (center_pad*'N')
            acceptor_seq = (center_pad*'N') + str(acceptor_seq) + (right_pad*'N')

            # save pos in df
            df.at[idx, 'donor_s'] = first_s
            df.at[idx, 'donor_e'] = first_e
            df.at[idx, 'acceptor_s'] = second_s
            df.at[idx, 'acceptor_e'] = second_e

        elif strand == '-':
            donor_seq = seq[second_s:second_e].reverse.complement
            acceptor_seq = seq[first_s:first_e].reverse.complement

            # add padding
            donor_seq = (right_pad*'N') + str(donor_seq) + (center_pad*'N')
            acceptor_seq = (center_pad*'N') + str(acceptor_seq) + (left_pad*'N')

            # save pos in df
            df.at[idx, 'donor_s'] = second_s
            df.at[idx, 'donor_e'] = second_e
            df.at[idx, 'acceptor_s'] = first_s
            df.at[idx, 'acceptor_e'] = first_e

        df.at[idx, 'donor_seq'] = donor_seq
        df.at[idx, 'acceptor_seq'] = acceptor_seq

        # sanity checks
        try:
            assert(len(donor_seq) == 400)
            assert(len(acceptor_seq) == 400)
            assert(donor_seq[200:202] == df.at[idx, 'donor_dimer'])
            assert(acceptor_seq[198:200] == df.at[idx, 'acceptor_dimer'])
        except AssertionError as e:
            print('Debug values: ', start, end, strand, flanking_size, left_pad, right_pad, center_pad, 
                  donor_seq[200:202], acceptor_seq[198:200], df.at[idx, 'donor_dimer'], df.at[idx, 'acceptor_dimer'])
            raise(e)

        pbar.next()
    pbar.finish()
        
    return df

'''Write all to files'''
def write_results(df, fw_donor, fw_acceptor, fw_donor_seq, fw_acceptor_seq, fw_da, fw_input):

    print('Writing results...')
    
    # write the BED files
    da_columns = ['seqid', 'start', 'end', 'name', 'score', 'strand']
    df[da_columns].to_csv(fw_da, sep='\t', header=None, index=0) 

    donor_columns = ['seqid', 'donor_s', 'donor_e', 'name', 'score', 'strand']
    df[donor_columns].to_csv(fw_donor, sep='\t', header=None, index=0) 

    acceptor_columns = ['seqid', 'acceptor_s', 'acceptor_e', 'name', 'score', 'strand']
    df[acceptor_columns].to_csv(fw_acceptor, sep='\t', header=None, index=0)         

    # write the FASTA files
    for idx, row in df.iterrows():

        chrom = row['seqid']
        strand = row['strand']
        donor_s = row['donor_s']
        donor_e = row['donor_e']
        donor_seq = row['donor_seq']
        acceptor_s = row['acceptor_s']
        acceptor_e = row['acceptor_e']
        acceptor_seq = row['acceptor_seq']
        start = row['start']
        end = row['end']
        input_seq = donor_seq + acceptor_seq

        fw_donor_seq.write(f'>{chrom}:{donor_s+1}-{donor_e}({strand})\n{donor_seq}\n')
        fw_acceptor_seq.write(f'>{chrom}:{acceptor_s+1}-{acceptor_e}({strand})\n{acceptor_seq}\n')

        # NOTE: the first position is 0-indexed, second position is 1-indeed for SPLAM compatibility
        fw_input.write(f'>{chrom};{start};{end};{strand};1\n{input_seq}\n')

'''Show dimer frequencies'''
def display_stats(df):
    
    donor_cts = df['donor_dimer'].value_counts()
    acceptor_cts = df['acceptor_dimer'].value_counts()
    
    count_df = pd.DataFrame({'Donor Dimer Counts': donor_cts, 'Acceptor Dimer Counts': acceptor_cts}).fillna(0).astype('int')

    print(count_df)


def main(db):

    # make the necessary dirs
    print(f'Parsing for {db} dataset:')

    # inputs 
    bed_file = f'./1_output/{db}_introns.bed'
    fasta_file = f'../data/{db}_genomic.fa'
    assembly_file = f'../data/{db}_assembly_report.txt'
    
    # outputs
    output_dir = f'./2_output/{db}/'
    os.makedirs(output_dir, exist_ok=True)
    fw_donor = open(f'{output_dir}donor.bed', 'w')
    fw_acceptor = open(f'{output_dir}acceptor.bed', 'w')
    fw_donor_seq = open(f'{output_dir}donor_seq.fa', 'w')
    fw_acceptor_seq = open(f'{output_dir}acceptor_seq.fa', 'w')
    fw_da = open(f'{output_dir}d_a.bed', 'w')
    fw_input = open(f'{output_dir}input_pos.fa', 'w')
    
    # process inputs
    print('Processing inputs...')
    df = pd.read_csv(bed_file, delimiter='\t', header=None, usecols=range(6),
                     names=['seqid', 'start', 'end', 'name', 'score', 'strand'])
    seqs = Fasta(fasta_file, sequence_always_upper=True)
    chrs = get_chrom_size(assembly_file)

    # remove any duplicate isoforms in the bed file
    len_before = len(df)
    df.drop_duplicates(subset=['seqid', 'start', 'end', 'strand'], inplace=True)
    print(f'{len_before} introns before\n{len(df)} introns after\n{len_before-len(df)} duplicate introns removed')

    # extract the sequences
    df = extract(df, seqs, chrs)

    # write results to files
    write_results(df, fw_donor, fw_acceptor, fw_donor_seq, fw_acceptor_seq, fw_da, fw_input)
    print(f'Result in {output_dir}')

    # display dimer statistics
    display_stats(df)

    fw_donor.close()
    fw_acceptor.close()
    fw_donor_seq.close()
    fw_acceptor_seq.close()
    fw_da.close()
    fw_input.close()


if __name__ == "__main__":

    if os.getcwd() != 'pos_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_generalization_test/pos_test/')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [0,1,2,3] #CHANGEME

    for idx in idxs:
        main(datasets[idx])