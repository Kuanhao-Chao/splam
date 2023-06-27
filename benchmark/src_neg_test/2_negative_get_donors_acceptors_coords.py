from Bio import SeqIO
import random
import os
from progress.bar import Bar
import time
import pyfaidx

SEQ_LENGTH = "800"
QUARTER_SEQ_LEN = int(SEQ_LENGTH) // 4
EACH_JUNC_PER_LOCUS = 5
MIN_JUNC = 200
MAX_JUNC = 20000
THRESHOLD = "100"   
OUTPUT_DIR = f'./2_output/{SEQ_LENGTH}bp/'

'''Make the directory(ies) needed for the given path'''
dir = lambda path : os.makedirs(os.path.dirname(path), exist_ok=True)

'''Obtain the negative sequence'''
def task(chromosome, sequence, start, end, strand, fw_donor, fw_acceptor, fw_da):

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
        donor_idx = 0
        acceptor_idx = JUNC_BASE
        
        has_donor = False
        has_acceptor = False

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
            donor_dimer = 'CT'
            acceptor_dimer = 'AC'

        elif strand_select == "+":
            # I need to find junctions on the forward strand (GT-AG pairs).
            donor_dimer = 'GT'
            acceptor_dimer = 'AG'

        
        ################################
        # Finding the donor dimer
        ################################
        has_donor = False
        while has_donor == False:
            # (1) Found the donor dimer
            cur_dimer = sequence[base_num+donor_idx:base_num+donor_idx+2]
            #print(f'Donor dimer: {cur_dimer}')
            if (cur_dimer == donor_dimer):
                has_donor = True

            # (2) Donor is out of range
            if base_num+donor_idx+2 > end:
                #print('\tDonor out of range.')
                break

            donor_idx += 1

        if not has_donor:
            continue
        

        ################################
        # Finding the acceptor dimer
        ################################
        has_acceptor = False
        while has_acceptor == False:
            # (1) Found the acceptor dimer
            cur_dimer = sequence[base_num+donor_idx+acceptor_idx-2:base_num+donor_idx+acceptor_idx]
            #print(f'Acceptor dimer: {cur_dimer}')
            if (cur_dimer == acceptor_dimer):
                has_acceptor = True

            # (2) Acceptor is out of range.
            if base_num+donor_idx+acceptor_idx > end:
                #print('\tAcceptor out of range.')
                break

            acceptor_idx += 1

        if not has_acceptor:
            continue


        ########################################
        # Found both donor and acceptor dimers
        ########################################
        if has_donor and has_acceptor:

            # print(f'Donor   : {sequence[base_num+donor_idx:base_num+donor_idx+2]}')
            # print(f'Acceptor: {sequence[base_num+donor_idx+acceptor_idx-2:base_num+donor_idx+acceptor_idx]}')

            # generate the 200bp flanking sequence each side
            splice_junc_len = acceptor_idx + 1
            flanking_size = QUARTER_SEQ_LEN

            if splice_junc_len < QUARTER_SEQ_LEN:
                flanking_size = splice_junc_len

            donor = base_num+donor_idx
            acceptor = base_num+donor_idx+acceptor_idx

            if (strand_select == "+"):
                donor_s = donor - QUARTER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + QUARTER_SEQ_LEN

            elif (strand_select == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + QUARTER_SEQ_LEN
                acceptor_s = acceptor - QUARTER_SEQ_LEN
                acceptor_e = acceptor + flanking_size

            ######################################################
            # Check if the donors and acceptors are in range.
            ######################################################
            if acceptor_s >= end:
                continue

            fw_da.write(chromosome + "\t" + str(base_num+donor_idx) + "\t" + str(base_num+donor_idx+acceptor_idx+1) + "\t" + "JUNC_" + str(junc_count) + "\t1\t"+strand_select+"\n")
            
            if (strand_select == "+"):
                fw_donor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_" + str(junc_count) + "\t1\t"+strand_select+"\n")
                fw_acceptor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_" + str(junc_count) + "\t1\t"+strand_select+"\n")

            elif (strand_select == "-"):
                fw_acceptor.write(chromosome + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + "JUNC_" + str(junc_count) + "\t1\t"+strand_select+"\n")
                fw_donor.write(chromosome + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + "JUNC_" + str(junc_count) + "\t1\t"+strand_select+"\n")
            
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

    with open(fasta_file, 'r') as handle, open(bed_file) as f:

        record_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        #print(record_dict)

        lines = f.read().splitlines()
        pbar = Bar('Generating negative samples... ', max=len(lines))
        for line in lines:
            
            s = time.time()

            eles = line.split("\t")

            chr = eles[0]
            start = int(eles[1])
            end = int(eles[2])
            strand = eles[5]
            record = record_dict[chr]

            # Extract individual parts of the FASTA record
            chromosome = record.description
            sequence = record.seq
            e = time.time()
            sequence = str(sequence).upper()

            

            #print(f'Searching in: {chromosome}:{start}-{end};{strand}')
            s1 = time.time()
            task(chromosome, sequence, start, end, strand, fw_donor, fw_acceptor, fw_da)
            s2 = time.time()

            print(f'Compare: {e-s}, task {s2-s1}')

            pbar.next()
        pbar.finish()
    
    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":

    if os.getcwd() != 'src_neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_neg_test')

    datasets = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10.1']
    idxs = [0] #CHANGEME

    for idx in idxs:
        main(datasets[idx])