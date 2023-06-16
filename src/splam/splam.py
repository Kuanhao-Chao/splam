import pandas as pd
import os 
import argparse
import sys

from splam import prediction, config, parse, chr_size

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

def main(argv=None):
    args = parse_args(argv)
    print(args)
    print(args.output)
    print(args.reference_genome)

    chrs = chr_size.chrs
    junction_bed = args.junction_BED
    junction_dir = os.path.dirname(junction_bed)
    junction_score_bed = junction_dir + "/junction_score.bed"
    reference_genome = args.reference_genome
    splam_model = args.model

    #################################
    # Step 1: creating donor acceptor bed file.
    #################################
    donor_bed, acceptor_bed = parse.create_donor_acceptor_bed(junction_bed, junction_dir, chrs)

    #################################
    # Step 2: write donor acceptor fasta file.
    #################################
    donor_fasta = parse.write_donor_acceptor_fasta(donor_bed, reference_genome)
    acceptor_fasta = parse.write_donor_acceptor_fasta(acceptor_bed, reference_genome)

    #################################
    # Step 3: concatenate donor and acceptor into a fasta
    #################################
    junction_fasta = parse.concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta)
    
    #################################
    # Step 4: splam score junctions
    #################################
    out_score_f = os.path.dirname(junction_fasta) + "/junction_score.bed"
    junction_fasta = prediction.splam_prediction(junction_fasta, junction_score_bed, splam_model)
    print("junction_fasta: ", junction_fasta)

if __name__ == "__main__":
    main()