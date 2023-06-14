import pandas as pd
import os 
import argparse
import sys
import pybedtools

from prediction import *
from helper import *

threshold = "100"
SEQ_LEN="800"
QUARTER_SEQ_LEN = int(SEQ_LEN) // 4
HALF_SEQ_LEN = int(SEQ_LEN) // 2

def parse_args(args):
    parser = argparse.ArgumentParser(description='score splice junctions from a BED file')
    parser.add_argument("junction_BED", help="target splice junctions in bed files.")

    parser.add_argument(
        '-o', '--output', default="score.bed", metavar='score.bed',
        help='write output to file in same format as input; by default, output is written to "score.bed"'
    )

    parser.add_argument(
        '-G', '--reference-genome',  metavar='reference_genome.fasta',
        help='The path to the reference genome.'
    )

    parser.add_argument(
        '-m', '--model', metavar='<model.pt>', 
        required=True, help='the path to the SPLAM! model'
    )

    args_r = parser.parse_args(args)

    return args_r

def create_donor_acceptor_bed(junction_bed, junction_dir, chrs):
    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    donor_bed = junction_dir+"/juncs/donor.bed"
    acceptor_bed = junction_dir+"/juncs/acceptor.bed"
    os.makedirs(junction_dir+"/juncs/", exist_ok=True)
    fw_donor = open(donor_bed, "w")
    fw_acceptor = open(acceptor_bed, "w")
    
    d_a_bed = junction_dir+"/juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()

    with open(junction_bed, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            if len(eles) == 1:
                continue
            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]

            
            if (strand == "+"):
                donor = int(eles[1])
                acceptor = int(eles[2])
                splice_junc_len = acceptor - donor
            elif (strand == "-"):
                acceptor = int(eles[1])
                donor = int(eles[2])
                splice_junc_len = donor - acceptor

            flanking_size = QUARTER_SEQ_LEN
            if splice_junc_len < QUARTER_SEQ_LEN:
                flanking_size = splice_junc_len

            if (strand == "+"):
                donor_s = donor - QUARTER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + QUARTER_SEQ_LEN

            elif (strand == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + QUARTER_SEQ_LEN
                acceptor_s = acceptor - QUARTER_SEQ_LEN
                acceptor_e = acceptor + flanking_size


            if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
                continue
            if donor_s < 0 or acceptor_s < 0:
                continue
            new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
            if new_junc in JUNCS:
                continue
            else:
                JUNCS.add(new_junc)
                fw_donor.write(chr + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
                fw_acceptor.write(chr + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")

                if (strand == "+"):
                    fw_da.write(chr + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
                elif (strand == "-"):
                    fw_da.write(chr + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()
    return donor_bed, acceptor_bed

def write_donor_acceptor_fasta(bed_file, reference_genome):

    # Create a BedTool object from the BED file
    bed_f = pybedtools.BedTool(bed_file)

    # Get the sequences from the BED intervals
    sequences = bed_f.sequence(fi=reference_genome, s=True)

    fasta_f = os.path.splitext(bed_file)[0]+ ".fasta"
    print("fasta_f: ", fasta_f)

    b = sequences.save_seqs(fasta_f)
    return fasta_f

def concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta):
    output_file = os.path.dirname(donor_fasta) + "/junction.fa"
    fw = open(output_file, "w")
    fr_donor = open(donor_fasta, "r")
    fr_acceptor = open(acceptor_fasta, "r")

    # parsing donor and acceptor fa files
    lines_d = fr_donor.read().splitlines()
    lines_a = fr_acceptor.read().splitlines()
    line_num = len(lines_d)

    # initializing stats
    canonical_d_count = 0 # GT
    noncanonical_d_count = 0
    canonical_a_count = 0 # AG
    noncanonical_a_count = 0
    donors = {}
    acceptors = {}
    num_skipped = 0
    for idx in range(0, line_num, 2):
        # PARSE FIRST LINE
        # >chr1:10000-20000(+)
        chr_name = lines_d[idx]
        strand = lines_d[idx][-2]
        chromosome = lines_d[idx].split(":")[0]

        d_splits = lines_d[idx].split(":")[1].split("(")
        d_start, d_end = d_splits[0].split("-")
        d_strand = d_splits[1][0]

        a_splits = lines_a[idx].split(":")[1].split("(")
        a_start, a_end = a_splits[0].split("-")
        a_strand = a_splits[1][0]

        # print("donor   : ", d_start, d_end)
        # print("d_strand: ", d_strand)

        # print("Acceptor   : ", a_start, a_end)
        # print("a_strand: ", a_strand)

        if strand == "+":
            donor_pos = int(d_start) + QUARTER_SEQ_LEN
            acceptor_pos = int(a_end) - QUARTER_SEQ_LEN
        elif strand == "-":
            donor_pos = int(d_end) - QUARTER_SEQ_LEN
            acceptor_pos = int(a_start) + QUARTER_SEQ_LEN

        # PARSE SECOND LINE
        idx += 1

        seq_d = lines_d[idx]
        seq_a = lines_a[idx]
        len_d = len(seq_d)
        len_a = len(seq_a)

        if len_d != len_a:
            print(f'Unequal lengths: seq_d {len_d}, seq_a {len_a}')

        if len_d == HALF_SEQ_LEN and len_a == HALF_SEQ_LEN:
            # combine normally
            x = seq_d + seq_a
        else:
            # pad with repeating N
            x = seq_d + (HALF_SEQ_LEN - len_d) * 'N' + (HALF_SEQ_LEN - len_a) * 'N' + seq_a

        x = x.upper()
        if (len(x) != int(SEQ_LEN)): 
            print("x: ", len(x))
       
        # skip sequence if there are Ns in the sequence
        if x[QUARTER_SEQ_LEN] == "N" or x[QUARTER_SEQ_LEN+1] == "N" or x[QUARTER_SEQ_LEN*3-2] == "N" or x[QUARTER_SEQ_LEN*3-1] == "N":
            num_skipped += 1
            continue
        
        # write the final fasta entry as two lines
        fw.write(chromosome + ";" + str(donor_pos) +";"+ str(acceptor_pos) + ";" + d_strand + ";1\n")
        fw.write(x + "\n")

        # get stats on the dimers 
        donor_dimer = x[QUARTER_SEQ_LEN:QUARTER_SEQ_LEN+2]
        acceptor_dimer = x[QUARTER_SEQ_LEN*3-2:QUARTER_SEQ_LEN*3]

        if donor_dimer not in donors.keys():
            donors[donor_dimer] = 1
        else:
            donors[donor_dimer] += 1

        if acceptor_dimer not in acceptors.keys():
            acceptors[acceptor_dimer] = 1
        else:
            acceptors[acceptor_dimer] += 1

        if (donor_dimer == "GT"):
            canonical_d_count += 1
        else:
            noncanonical_d_count += 1
            
        if (acceptor_dimer == "AG"):
            canonical_a_count += 1
        else:
            noncanonical_a_count += 1

    # output stats
    print("Number of skips due to N in dimer: ", num_skipped)

    print("Canonical donor count: ", canonical_d_count)
    print("Noncanonical donor count: ", noncanonical_d_count)
    print("Canonical acceptor count: ", canonical_a_count)
    print("Noncanonical acceptor count: ", noncanonical_a_count)
    for key, value in donors.items():
        print("Donor   : ", key, " (", value, ")")
    for key, value in acceptors.items():
        print("Acceptor: ", key, " (", value, ")")

    fw.close()
    fr_acceptor.close()
    fr_donor.close()

    return output_file

def main(argv=None):
    args = parse_args(argv)
    print(args)
    print(args.output)
    print(args.reference_genome)

    chrs = get_chrom_size('GRCh38.p14_assembly_report.txt')

    junction_bed = args.junction_BED
    junction_dir = os.path.dirname(junction_bed)
    junction_score_bed = junction_dir + "/junction_score.bed"
    reference_genome = args.reference_genome
    splam_model = args.model

    #################################
    # Step 1: creating donor acceptor bed file.
    #################################
    donor_bed, acceptor_bed = create_donor_acceptor_bed(junction_bed, junction_dir, chrs)

    #################################
    # Step 2: write donor acceptor fasta file.
    #################################
    donor_fasta = write_donor_acceptor_fasta(donor_bed, reference_genome)
    acceptor_fasta = write_donor_acceptor_fasta(acceptor_bed, reference_genome)

    #################################
    # Step 3: concatenate donor and acceptor into a fasta
    #################################
    junction_fasta = concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta)
    
    #################################
    # Step 4: splam score junctions
    #################################
    out_score_f = os.path.dirname(junction_fasta) + "/junction_score.bed"
    junction_fasta = splam_prediction(junction_fasta, junction_score_bed, splam_model)
    print("junction_fasta: ", junction_fasta)

if __name__ == "__main__":
    main()