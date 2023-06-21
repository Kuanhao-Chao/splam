import pandas as pd
import os 
import argparse
import sys

from splam import prediction, config, parse, chr_size
import splam_extract
import splam_clean


def parse_args(args):


    parser = argparse.ArgumentParser(prog='splam', description='\033[1;37msplice junction predictor to improve alignment files (BAM / CRAM)\033[0;0m')
    
    parser.add_argument('-V', '--verbose',
                    action='store_true')  # on/off flag
    parser.add_argument('-v', '--version',
                    action='store_true')  # on/off flag
    parser.add_argument(
        '-c', '--cite', action='store_true')
    

    subparsers = parser.add_subparsers(title='Commands', metavar= "subcommand: 'extract | score | clean'", dest='subcommand')#, required=True)

    # subparsers = parser.add_subparsers(title='Commands', dest='cmd', required=True)


    #############################
    # splam extract subcommands
    #############################
    parser_extract = subparsers.add_parser('extract', help='Extracting all splice junctions from a BAM file')

    #############################
    # splam score subcommands
    #############################
    parser_score= subparsers.add_parser('score', help='Scoring all splice junctions')

    parser_score.add_argument("junction_BED", help="target splice junctions in bed files.")

    parser_score.add_argument('-v', '--verbose',
                    action='store_true')  # on/off flag
    parser_score.add_argument(
        '-o', '--outdir', default="tmp", metavar='score.bed',
        help='write output to file in same format as input; by default, output is written to "score.bed"'
    )
    parser_score.add_argument(
        '-G', '--reference-genome',  metavar='reference_genome.fasta',
        required=True, help='The path to the reference genome.'
    )
    parser_score.add_argument(
        '-m', '--model', metavar='model.pt', 
        required=True, help='the path to the SPLAM! model'
    )
    
    
    #############################
    # splam clean subcommands
    #############################
    parser_clean = subparsers.add_parser('clean', help='Cleaning up spurious splice alignment')


    args_r = parser.parse_args(args)

    return args_r

def main(argv=None):    

    print(
            "====================================================================\n"
            "An accurate spliced alignment pruner and spliced junction predictor.\n"
            "====================================================================\n");
    print("""
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║
  ███████╗██████╔╝██║     ███████║██╔████╔██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝
    """)
    chrs = chr_size.chrs
    args = parse_args(argv)

    print("args: ", args)

    if args.subcommand == "extract":
        # print("extract!!")
        splam_extract.splam_extract(["splam-extract"])
    elif args.subcommand == "score":
        # print("score!!")
        verbose =args.verbose
        outdir = args.outdir
        junction_bed = args.junction_BED
        junction_score_bed = os.path.join(outdir, "junction_score.bed")
        reference_genome = args.reference_genome
        splam_model = args.model

        #################################
        # Step 1: creating donor acceptor bed file.
        #################################
        donor_bed, acceptor_bed = parse.create_donor_acceptor_bed(junction_bed, outdir, chrs)

        #################################
        # Step 2: write donor acceptor fasta file.
        #################################
        donor_fasta = parse.write_donor_acceptor_fasta(donor_bed, reference_genome)
        acceptor_fasta = parse.write_donor_acceptor_fasta(acceptor_bed, reference_genome)

        #################################
        # Step 3: concatenate donor and acceptor into a fasta
        #################################
        junction_fasta = parse.concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta, verbose)
        
        #################################
        # Step 4: splam score junctions
        #################################
        out_score_f = os.path.dirname(junction_fasta) + "/junction_score.bed"
        junction_fasta = prediction.splam_prediction(junction_fasta, junction_score_bed, splam_model)

    elif args.subcommand == "clean":
        # print("clean!!")
        splam_clean.splam_clean(["splam-clean", "-o", "SRR1352129_chr9_sub"])

if __name__ == "__main__":
    main()