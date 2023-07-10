import random
import os
from progress.bar import Bar
from pyfaidx import Fasta

SEQ_LENGTH = "800"
QUARTER_SEQ_LEN = int(SEQ_LENGTH) // 4
EACH_JUNC_PER_LOCUS = 20 #15 for mammalian
MIN_JUNC = 200
MAX_JUNC = 20000
OUTPUT_DIR = f'./2_output/{SEQ_LENGTH}bp/'

'''Make the directory(ies) needed for the given path'''
dir = lambda path : os.makedirs(os.path.dirname(path), exist_ok=True)

 
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


'''Obtain the negative sequence'''
def task(chromosome, sequence, start, end, name, strand, size, fw_donor, fw_acceptor, fw_da):

    junc_count = 0
    # does EACH_JUNC_PER_LOCUS number of seed generations
    for idx in range(EACH_JUNC_PER_LOCUS):

        # random seed to start searching 
        base_num = random.randint(start, end)
        # base length of junction within range MIN_JUNC and MAX_JUNC
        JUNC_BASE = random.randint(MIN_JUNC, MAX_JUNC)
        

        ################################
        # Finding the first and second dimers
        ################################
        first_idx = base_num
        second_idx = JUNC_BASE
        
        has_first = False
        has_second = False

        # select the opposite strand
        strand_select = "."
        if strand == "+":
            strand_select = "-"
        elif strand == "-":
            strand_select = "+"
        else:
            #print('No strand.')
            continue

        if strand_select == "-":
            # I need to find junctions on the reverse strand (CT-AC pairs).
            first_dimer = 'CT' # acceptor
            second_dimer = 'AC' # donor

        elif strand_select == "+":
            # I need to find junctions on the forward strand (GT-AG pairs).
            first_dimer = 'GT' # donor
            second_dimer = 'AG' # acceptor

        
        ################################
        # Finding the first dimer
        ################################
        has_first = False
        while has_first == False:
            # (1) Found the first (leftmost) dimer
            cur_dimer = sequence[first_idx:first_idx+2]
            #print(f'Donor dimer: {cur_dimer}')
            if (cur_dimer == first_dimer):
                has_first = True
                break

            # (2) Dimer is out of range
            if first_idx+2 > end:
                #print('\tDonor out of range.')
                break

            first_idx += 1

        if not has_first:
            continue
        

        ################################
        # Finding the second dimer
        ################################
        has_second = False
        while has_second == False:
            # (1) Found the second (rightmost) dimer
            cur_dimer = sequence[first_idx+second_idx-2:first_idx+second_idx]
            #print(f'Acceptor dimer: {cur_dimer}')
            if (cur_dimer == second_dimer):
                has_second = True
                break

            # (2) Dimer is out of range.
            if first_idx+second_idx > end:
                #print('\tAcceptor out of range.')
                break

            second_idx += 1

        if not has_second:
            continue


        ####################
        # Found both dimers
        ####################
        if has_first and has_second:

            # print(f'Donor   : {sequence[first_idx:first_idx+2]}')
            # print(f'Acceptor: {sequence[first_idx+second_idx-2:first_idx+second_idx]}')

            # generate the 200bp flanking sequence each side (decreased if fully overlaps)
            splice_junc_len = second_idx
            flanking_size = QUARTER_SEQ_LEN

            if splice_junc_len < QUARTER_SEQ_LEN:
                flanking_size = splice_junc_len

            first_pos = first_idx
            second_pos = first_idx+second_idx

            first_s = first_pos - QUARTER_SEQ_LEN
            first_e = first_pos + flanking_size
            second_s = second_pos - flanking_size
            second_e = second_pos + QUARTER_SEQ_LEN

            ###################################
            # Check if the dimers are in range
            ###################################
            if second_e >= size or first_s < 0:
                continue

            # sanity checks
            assert(first_e > first_s)
            assert(second_e > second_s)
            assert(second_e > first_s)
            if strand_select == '+':
                assert(sequence[first_pos:first_pos+2] == 'GT')
                assert(sequence[second_pos-2:second_pos] == 'AG')
            if strand_select == '-':
                assert(sequence[first_pos:first_pos+2] == 'CT')
                assert(sequence[second_pos-2:second_pos] == 'AC')

            # write to file
            fw_da.write(f'{chromosome}\t{first_idx}\t{first_idx+second_idx}\t{name}__{junc_count}\t0\t{strand_select}\n')
            
            if (strand_select == "+"):
                fw_donor.write(f'{chromosome}\t{first_s}\t{first_e}\t{name}__{junc_count}_donor\t0\t{strand_select}\n')
                fw_acceptor.write(f'{chromosome}\t{second_s}\t{second_e}\t{name}__{junc_count}_acceptor\t0\t{strand_select}\n')

            elif (strand_select == "-"):
                fw_acceptor.write(f'{chromosome}\t{first_s}\t{first_e}\t{name}__{junc_count}_donor\t0\t{strand_select}\n')
                fw_donor.write(f'{chromosome}\t{second_s}\t{second_e}\t{name}__{junc_count}_acceptor\t0\t{strand_select}\n')
            
            junc_count += 1


def main(db):

    # make the necessary dirs
    print(f'Parsing for {db} dataset')
    dir(OUTPUT_DIR+db+'/')

    # inputs 
    fasta_file = f'../SPLAM_python/extraction/primates/{db}_genomic.fa'
    bed_file = f'./1_output/{db}_genes.bed'
    
    # outputs
    fw_donor = open(f'{OUTPUT_DIR}{db}/donor.bed', 'w')
    fw_acceptor = open(f'{OUTPUT_DIR}{db}/acceptor.bed', 'w')
    fw_da = open(f'{OUTPUT_DIR}{db}/d_a.bed', 'w')

    with open(bed_file) as f:

        record_dict = Fasta(fasta_file, sequence_always_upper=True)
        #print(record_dict)
        size_dict = get_chrom_size(f'../SPLAM_python/extraction/primates/{db}_assembly_report.txt')

        lines = f.read().splitlines()
        pbar = Bar('Generating negative samples... ', max=len(lines))
        for line in lines:
            eles = line.split("\t")
            #print(eles)

            chr = eles[0]
            start = int(eles[1])
            end = int(eles[2])
            name = eles[3]
            strand = eles[5]
            sequence = record_dict[chr]
            size = size_dict[chr]
            #chromosome = sequence.long_name

            # Extract individual parts of the FASTA record
        
            #print(f'Searching in: {chromosome}:{start}-{end};{strand}')
            task(chr, sequence, start, end, name, strand, size, fw_donor, fw_acceptor, fw_da)

            pbar.next()
        pbar.finish()
    
    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()


if __name__ == "__main__":

    if os.getcwd() != 'src_neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_neg_test')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [0,1,2] #CHANGEME

    for idx in idxs:
        main(datasets[idx])