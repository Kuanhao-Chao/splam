import pandas as pd
import os 
import argparse
import sys

from splam import prediction, config, parse, chr_size
import splam_extract

def parse_args(args):
    parser = argparse.ArgumentParser(description='score splice junctions from a BED file')
    parser.add_argument("junction_BED", help="target splice junctions in bed files.")

    parser.add_argument('-v', '--verbose',
                    action='store_true')  # on/off flag

    parser.add_argument(
        '-o', '--outdir', default="tmp", metavar='score.bed',
        help='write output to file in same format as input; by default, output is written to "score.bed"'
    )

    parser.add_argument(
        '-G', '--reference-genome',  metavar='reference_genome.fasta',
        required=True, help='The path to the reference genome.'
    )

    parser.add_argument(
        '-m', '--model', metavar='model.pt', 
        required=True, help='the path to the SPLAM! model'
    )

    args_r = parser.parse_args(args)

    return args_r

def main(argv=None):
    print("Hello world!")

    splam_extract.splam_extract(["splam-extract", "-o", "SRR1352129_chr9_sub", "--paired", "../../Dataset/SRR1352129_chr9_sub/SRR1352129_chr9_sub.bam"])
    
    # chrs = chr_size.chrs
    # args = parse_args(argv)
    # verbose =args.verbose
    # outdir = args.outdir

    # junction_bed = args.junction_BED
    # junction_score_bed = os.path.join(outdir, "junction_score.bed")
    # reference_genome = args.reference_genome
    # splam_model = args.model

    # #################################
    # # Step 1: creating donor acceptor bed file.
    # #################################
    # donor_bed, acceptor_bed = parse.create_donor_acceptor_bed(junction_bed, outdir, chrs)

    # #################################
    # # Step 2: write donor acceptor fasta file.
    # #################################
    # donor_fasta = parse.write_donor_acceptor_fasta(donor_bed, reference_genome)
    # acceptor_fasta = parse.write_donor_acceptor_fasta(acceptor_bed, reference_genome)

    # #################################
    # # Step 3: concatenate donor and acceptor into a fasta
    # #################################
    # junction_fasta = parse.concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta, verbose)
    
    # #################################
    # # Step 4: splam score junctions
    # #################################
    # out_score_f = os.path.dirname(junction_fasta) + "/junction_score.bed"
    # junction_fasta = prediction.splam_prediction(junction_fasta, junction_score_bed, splam_model)

if __name__ == "__main__":
    main()