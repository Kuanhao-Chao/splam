import random
import os
from progress.bar import Bar
from pyfaidx import Fasta

SEQ_LENGTH = "800"
QUARTER_SEQ_LEN = int(SEQ_LENGTH) // 4
EACH_JUNC_PER_LOCUS = 20
MIN_JUNC = 200
MAX_JUNC = 20000
 
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
def extract(chromosome, sequence, start, end, name, strand, chrom_size, 
         fw_donor, fw_acceptor, fw_donor_seq, fw_acceptor_seq, fw_da, fw_input):
# NOTE: start is 0-indexed, end is 1-indexed (BED format)

    junc_count = 0
    # does EACH_JUNC_PER_LOCUS number of seed generations
    for attempt in range(EACH_JUNC_PER_LOCUS):
        
        #####################
        # Search Setup
        #####################
        # random seed to start searching 
        base_num = random.randint(start, end)
        # base length of junction within range MIN_JUNC and MAX_JUNC
        JUNC_BASE = random.randint(MIN_JUNC, MAX_JUNC)
        
        # establish the first and second dimer start positions
        first_idx = base_num
        second_idx = JUNC_BASE

        # select the opposite strand
        if strand == "+":
            strand_select = "-"

            # I need to find junctions on the reverse strand (CT-AC pairs).
            first_dimer = 'CT' # acceptor
            second_dimer = 'AC' # donor
        elif strand == "-":
            strand_select = "+"

            # I need to find junctions on the forward strand (GT-AG pairs).
            first_dimer = 'GT' # donor
            second_dimer = 'AG' # acceptor
        else:
            #print('No strand.')
            continue
            
        ################################
        # Finding the first dimer
        ################################
        has_first = False
        while first_idx+2 < end:
            
            # Found the first (left) dimer
            cur_dimer = sequence[first_idx:first_idx+2]
            #print(f'Donor dimer: {cur_dimer}')
            if cur_dimer == first_dimer:
                has_first = True
                break

            first_idx += 1

        # out of bounds: fail attempt
        if not has_first:
            continue

        ################################
        # Finding the second dimer
        ################################
        has_second = False
        while first_idx+second_idx < end:
            
            # Found the second (right) dimer
            cur_dimer = sequence[first_idx+second_idx-2:first_idx+second_idx]
            #print(f'Acceptor dimer: {cur_dimer}')
            if cur_dimer == second_dimer:
                has_second = True
                break

            second_idx += 1

        # out of bounds: fail attempt
        if not has_second:
            continue

        ######################
        # Found both dimers
        ######################
        if has_first and has_second:

            # generate the 200nt flanking sequence each side (decreased if fully overlaps)
            flanking_size = QUARTER_SEQ_LEN
            center_pad = 0
            if second_idx < QUARTER_SEQ_LEN:
                flanking_size = second_idx
                center_pad = QUARTER_SEQ_LEN - flanking_size

            first_pos = first_idx
            second_pos = first_idx + second_idx

            first_s = first_pos - QUARTER_SEQ_LEN
            first_e = first_pos + flanking_size
            second_s = second_pos - flanking_size
            second_e = second_pos + QUARTER_SEQ_LEN

            ###################################
            # Prevent dropping of boundary cases
            ###################################
            left_pad = 0
            right_pad = 0
            if first_s < 0:
                left_pad = 0 - first_s
                first_s = 0
            if second_e > chrom_size:
                right_pad = second_e - chrom_size
                second_e = chrom_size

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
           
            #######################################
            # Extract donor and acceptor sequences
            #######################################
            # NOTE: ends are 1-indexed, starts are 0-indexed, but converted to 1-indexed in FASTA format

            # get the donor and acceptor sequences
            if strand_select == '+':
                donor_seq = sequence[first_s:first_e]
                acceptor_seq = sequence[second_s:second_e]
                
                # add padding
                donor_seq = (left_pad*'N') + str(donor_seq) + (center_pad*'N')
                acceptor_seq = (center_pad*'N') + str(acceptor_seq) + (right_pad*'N')

            elif strand_select == '-':
                donor_seq = sequence[second_s:second_e].reverse.complement
                acceptor_seq = sequence[first_s:first_e].reverse.complement

                # add padding
                donor_seq = (right_pad*'N') + str(donor_seq) + (center_pad*'N')
                acceptor_seq = (center_pad*'N') + str(acceptor_seq) + (left_pad*'N')

            # full input sequence for SPLAM
            input_sequence = donor_seq + acceptor_seq
            
            # sanity checks
            assert(len(donor_seq) == 400)
            assert(len(acceptor_seq) == 400)
            assert(donor_seq[200:202] == 'GT')
            assert(acceptor_seq[198:200] == 'AG')

            ####################
            # Write to files
            ####################
            fw_da.write(f'{chromosome}\t{first_pos}\t{second_pos}\t{name}__{junc_count}\t0\t{strand_select}\n')
            
            if strand_select == "+":
                fw_donor.write(f'{chromosome}\t{first_s}\t{first_e}\t{name}__{junc_count}_donor\t0\t{strand_select}\n')
                fw_acceptor.write(f'{chromosome}\t{second_s}\t{second_e}\t{name}__{junc_count}_acceptor\t0\t{strand_select}\n')

                fw_donor_seq.write(f'>{chromosome}:{first_s+1}-{first_e}({strand_select})\n{donor_seq}\n')
                fw_acceptor_seq.write(f'>{chromosome}:{second_s+1}-{second_e}({strand_select})\n{acceptor_seq}\n')

            elif strand_select == "-":
                fw_acceptor.write(f'{chromosome}\t{first_s}\t{first_e}\t{name}__{junc_count}_donor\t0\t{strand_select}\n')
                fw_donor.write(f'{chromosome}\t{second_s}\t{second_e}\t{name}__{junc_count}_acceptor\t0\t{strand_select}\n')

                fw_acceptor_seq.write(f'>{chromosome}:{first_s+1}-{first_e}({strand_select})\n{acceptor_seq}\n')
                fw_donor_seq.write(f'>{chromosome}:{second_s+1}-{second_e}({strand_select})\n{donor_seq}\n')

            # NOTE: the first position is 0-indexed, second position is 1-indexed for SPLAM compatibility
            fw_input.write(f'>{chromosome};{first_pos};{second_pos};{strand_select};0\n{input_sequence}\n')

            junc_count += 1


def main(db):

    # make the necessary dirs
    print(f'Parsing for {db} dataset')

    # inputs 
    bed_file = f'./1_output/{db}_genes.bed'
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
    fw_input = open(f'{output_dir}input_neg_random.fa', 'w')

    with open(bed_file) as f:

        record_dict = Fasta(fasta_file, sequence_always_upper=True)
        size_dict = get_chrom_size(assembly_file)

        lines = f.read().splitlines()

        pbar = Bar('Generating negative samples... ', max=len(lines))
        for line in lines:
            eles = line.split("\t")

            chr = eles[0]
            start = int(eles[1])
            end = int(eles[2])
            name = eles[3]
            strand = eles[5]
            sequence = record_dict[chr]
            chrom_size = size_dict[chr]

            # Extract individual parts of the FASTA record
            extract(chr, sequence, start, end, name, strand, chrom_size, 
                 fw_donor, fw_acceptor, fw_donor_seq, fw_acceptor_seq, fw_da, fw_input)

            pbar.next()
        pbar.finish()
    
    fw_donor.close()
    fw_acceptor.close()
    fw_donor_seq.close()
    fw_acceptor_seq.close()
    fw_da.close()
    fw_input.close()


if __name__ == "__main__":

    if os.getcwd() != 'src_neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_neg_test/')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    idxs = [1,2,3] #CHANGEME

    for idx in idxs:
        main(datasets[idx])