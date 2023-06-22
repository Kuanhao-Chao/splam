import pandas as pd
import os 
import argparse
import sys

from splam import prediction, config, parse, chr_size
import splam_extract
import splam_clean

VERSION = "0.1.0"

CITATION = "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM"


def parse_args(args):

    parser = argparse.ArgumentParser(prog='splam', description='\033[1;37msplice junction predictor to improve alignment files (BAM / CRAM)\033[0;0m')
    
    # parser.add_argument('-V', '--verbose',
    #                 action='store_true')  # on/off flag
    parser.add_argument('-v', '--version',
                        action='store_true')  # on/off flag
    parser.add_argument('-c', '--citation', 
                        action='store_true')
    
    subparsers = parser.add_subparsers(title='Commands', dest='subcommand')#, required=True)


    #############################
    # Mode 1: splam extract subcommands
    #############################
    parser_extract = subparsers.add_parser('extract', help='Extracting all splice junctions from a BAM file')
    parser_extract.add_argument("BAM_INPUT", help="target alignment file in BAM format.")
    parser_extract.add_argument('-V', '--verbose',
                    action='store_true',
                    help='running splam in verbose mode')  # on/off flag
    parser_extract.add_argument('-P', '--paired',
                    action='store_true',
                    help='paired-end bundling alignments')  # on/off flag
    parser_extract.add_argument(
        '-n', '--write-junctions-only',
        action='store_true',
        help='only write out splice junction bed file without other temporary files.'
    )
    parser_extract.add_argument(
        '-o', '--outdir', default="tmp_splam_out", metavar='DIR',
        help='the directory where the output file is written to. Default output filename is "junction_score.bed"',
    )
    parser_extract.add_argument(
        '-M', '--max-splice',  metavar='DIST',
        help='maximum splice junction length'
    )
    parser_extract.add_argument(
        '-g', '--bundle-gap',  metavar='GAP',
        help='minimum gap between bundles'
    )


    #############################
    # Mode 2: splam score subcommands
    #############################
    parser_score= subparsers.add_parser('score', help='Scoring all splice junctions')
    parser_score.add_argument("junction_BED", help="target splice junctions in bed files.")
    parser_score.add_argument('-V', '--verbose',
                    action='store_true')  # on/off flag
    parser_score.add_argument(
        '-o', '--outdir', default="tmp_splam_out", metavar='DIR',
        help='the directory where the output file is written to. Default output filename is "junction_score.bed"',
    )
    parser_score.add_argument(
        '-G', '--reference-genome',  metavar='REF.fasta',
        required=True, help='The path to the reference genome.'
    )
    parser_score.add_argument(
        '-m', '--model', metavar='MODEL.pt', 
        required=True, help='the path to the SPLAM! model'
    )
    
    
    #############################
    # Mode 3: splam clean subcommands
    #############################
    parser_clean = subparsers.add_parser('clean', help='Cleaning up spurious splice alignment')
    parser_clean.add_argument(
        '-o', '--outdir', default="tmp_splam_out", metavar='DIR',
        help='the directory where the output file is written to. Default output filename is "junction_score.bed"',
        required=True
    )

    args_r = parser.parse_args()
    # args_r = parser.parse_args()
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

    if args.version:
        print(VERSION)
        exit()
    
    if args.citation:
        print(CITATION)
        exit()

    if args.subcommand == "extract":
        argv_extract = sys.argv
        argv_extract.pop(0)
        argv_extract[0] = 'splam-extract'
        splam_extract.splam_extract(argv_extract)

    elif args.subcommand == "score":
        verbose = args.verbose
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
        out_score_f = os.path.join(outdir, "/junction_score.bed")
        junction_fasta = prediction.splam_prediction(junction_fasta, junction_score_bed, splam_model)

    elif args.subcommand == "clean":
        argv_clean = sys.argv
        argv_clean.pop(0)
        argv_clean[0] = 'splam-clean'
        splam_clean.splam_clean(argv_clean)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()