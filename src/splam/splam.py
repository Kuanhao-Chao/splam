import pandas as pd
import os 
import argparse
import sys

from splam import prediction, config, parse, chr_size
import splam_extract
import splam_clean

VERSION = "0.1.0"

CITATION = "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM"

def msg(name=None):
    message = '''   splam-extract [arguments] BAM-file(s) \n\n\n
\033[1m\033[91mRequired argument:\033[0m
\t-o / --output\t\tPath to the output directory\n\n
\033[1m\033[94mOptional argument:\033[0m
\t-M / --max-splice <INT>:\tThe maximum intron length for the splice site.
\t-g / --bundle-gap <INT>:\tMinimum locus gap separation value. Reads that are mapped closer than this distance are merged together in the same processing bundle. Default: 100 (bp).\n\n
.
'''
    return message

def parse_args(args):

    parser = argparse.ArgumentParser(prog='splam', description='\033[1;37msplice junction predictor to improve alignment files (BAM / CRAM)\033[0;0m')
    
    # parser.add_argument('-V', '--verbose',
    #                 action='store_true')  # on/off flag
    parser.add_argument('-v', '--version',
                        action='store_true')  # on/off flag
    parser.add_argument('-c', '--citation', 
                        action='store_true')
    

    subparsers = parser.add_subparsers(title='Commands', dest='subcommand')#, required=True)

    # subparsers = parser.add_subparsers(title='Commands', dest='cmd', required=True)


    #############################
    # splam extract subcommands
    #############################
    parser_extract = subparsers.add_parser('extract', help='Extracting all splice junctions from a BAM file')


    #############################
    # splam score subcommands
    #############################
    parser_score= subparsers.add_parser('score', help='Scoring all splice junctions', usage=msg(), add_help=False)

    # parser.add_argument('-h', '--help', action='store_true', dest='help')

    parser_score.add_argument("junction_BED", help="target splice junctions in bed files.")

    parser_score.add_argument('-v', '--verbose',
                    action='store_true')  # on/off flag
    parser_score.add_argument(
        '-o', '--outdir', default="tmp", metavar='score.bed',
        help='write output to file in same fosprmat as input; by default, output is written to "score.bed"'
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


    # args_r = parser.parse_args()
    args_r = parser.parse_known_args()

    return args_r, parser, parser_score

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
    args, parser, parser_score = parse_args(argv)
    
    known_args = args[0]
    unknown_args = args[1]

    if known_args.version:
        print(VERSION)
        exit()
    
    if known_args.citation:
        print(CITATION)
        exit()

    if known_args.subcommand == "extract":
        print("argv: ", sys.argv)
        argv_extract = sys.argv
        argv_extract.pop(0)
        argv_extract[0] = 'splam-extract'
        print("argv_extract: ", argv_extract)
        splam_extract.splam_extract(argv_extract)

    elif known_args.subcommand == "score":
        # parser_score.print_help()
        # if args.help:
        #     print(">> Help!!")
        #     exit()
        print("known_args.help: ", known_args.help)
        if known_args.help:
            parser_score.print_help()

        verbose = known_args.verbose
        outdir = known_args.outdir
        junction_bed = known_args.junction_BED
        junction_score_bed = os.path.join(outdir, "junction_score.bed")
        reference_genome = known_args.reference_genome
        splam_model = known_args.model

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

    elif known_args.subcommand == "clean":
        # print("clean!!")
        splam_clean.splam_clean(["splam-clean", "-o", "SRR1352129_chr9_sub"])
    else:
        parser.print_help()

if __name__ == "__main__":
    main()