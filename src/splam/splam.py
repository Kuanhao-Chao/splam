import os
import argparse
import sys

from splam import prediction, config, parse, extract_gff, header, clean_gff
import splam_extract
import splam_clean

VERSION = header.__version__

CITATION = "Kuan-Hao Chao, Mihaela Pertea, and Steven Salzberg, \033[1m\x1B[3mSPLAM: accurate deep-learning-based splice site predictor to clean up spurious spliced alignments\x1B[0m\033[0m, (2023), GitHub repository, https://github.com/Kuanhao-Chao/SPLAM"


def parse_args(args):

    parser = argparse.ArgumentParser(prog='splam', description='\033[1;37msplice junction predictor to improve alignment files (BAM / CRAM)\033[0;0m')
    parser.add_argument('-v', '--version',
                        action='store_true')  # on/off flag
    parser.add_argument('-c', '--citation',
                        action='store_true')

    # Adding subcommand
    subparsers = parser.add_subparsers(title='Commands', dest='subcommand')


    #############################
    # Mode 1: splam extract subcommands
    #############################
    parser_extract = subparsers.add_parser('extract', help='Extracting all splice junctions from an alignment or annotation file')
    parser_extract.add_argument("INPUT", help="target alignment file in BAM format or annotation file in GFF format.")
    parser_extract.add_argument('-V', '--verbose',
                    action='store_true',
                    help='running splam in verbose mode.')  # on/off flag
    parser_extract.add_argument('-P', '--paired',
                    action='store_true',
                    help='bundling alignments in "paired-end" mode.')  # on/off flag
    parser_extract.add_argument(
        '-n', '--write-junctions-only',
        action='store_true',
        help='writing out splice junction bed file only without other temporary files.'
    )
    parser_extract.add_argument(
        '-f', '--file-format', default=None,
        help='the file type for SPLAM to process. It can only be "BAM", "GFF", or "GTF". The default value is "BAM".'
    )
    parser_extract.add_argument(
        '-d', '--database', default=None,
        help='the path to the annotation database built using gffutils. If thie argument is provided, splam loads the database instead of creating a new one.'
    )
    parser_extract.add_argument(
        '-o', '--outdir', default="tmp_out", metavar='DIR',
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
        '-o', '--outdir', default="tmp_out", metavar='DIR',
        help='the directory where the output file is written to. Default output filename is "junction_score.bed"',
    )
    parser_score.add_argument(
        '-b', '--batch-size', default=10, metavar='BATCH',
        help='the number of samples that will be propagated through the network. By default, the batch size is set to 10.'
    )
    parser_score.add_argument(
        '-d', '--device', default="NONE", metavar='pytorch_dev',
        help='the computing device that is used to perform computations on tensors and execute operations in the PyTorch framework. By default, this parameter is detectd automatically.'
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
        '-@', '--threads', default="1", metavar='threads',
        help='Set number of sorting, compression and merging threads. By default, operation is single-threaded.'
    )
    parser_clean.add_argument(
        '-t', '--threshold', default="0.1", metavar='threshold',
        help='The cutoff threshold for identifying spurious splice junctions.'
    )
    parser_clean.add_argument(
        '-n', '--bad-intron-num', default="8", metavar='bad intron num',
        help='The threshold for the number of spurious splice junctions in a transcript determines whether the transcript is considered bad. Default is 8.'
    )
    parser_clean.add_argument('-P', '--paired',
                    action='store_true',
                    help='cleaning up the alignment file in "paired-end" mode.')  # on/off flag
    parser_clean.add_argument(
        '-o', '--outdir', default="tmp_out", metavar='DIR',
        help='the directory where the output file is written to. Default output filename is "junction_score.bed".',
        required=True
    )

    args_r = parser.parse_args()
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

    args, parser, parser_score = parse_args(argv)

    if args.version:
        print(VERSION)
        exit()

    if args.citation:
        print(CITATION)
        exit()

    if args.subcommand == "extract":
        file_format = args.file_format
        input = args.INPUT
        if file_format is None:
            filename, file_extension = os.path.splitext(input)
            if file_extension == ".GTF" or file_extension == ".gtf":
                file_format = "GTF"
            elif file_extension == ".GFF" or file_extension == ".gff":
                file_format = "GFF"
            elif file_extension == ".BAM" or file_extension == ".bam":
                file_format = "BAM"

        if file_format == "GFF" or file_format == "GTF" or file_format == "gff" or file_format == "gtf":
            outdir = args.outdir
            junction_bed = os.path.join(outdir, "junction.bed")
            trans_intron_num_txt = os.path.join(outdir, "intron_num.txt")
            gff_db = args.database
            is_load_gff_db = False

            if gff_db is None:
                gff_db = os.path.join(outdir, "annotation.db")
            else:
                is_load_gff_db = True
            
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            extract_gff.extract_introns(input, gff_db, is_load_gff_db, junction_bed, trans_intron_num_txt)
        
        elif file_format == "BAM" or file_format == "bam":
            argv_extract = sys.argv
            argv_extract.pop(0)
            argv_extract[0] = 'splam-extract'
            splam_extract.splam_extract(argv_extract)

        else:
            print("[ERROR] the input file must be 'BAM', 'GFF', or 'GTF'.")
            sys.exit()


    elif args.subcommand == "score":
        verbose = args.verbose
        outdir = args.outdir
        junction_bed = args.junction_BED
        junction_score_bed = os.path.join(outdir, "junction_score.bed")
        reference_genome = args.reference_genome
        splam_model = args.model
        batch_size = args.batch_size
        device = args.device

        #################################
        # Step 1: creating donor acceptor bed file.
        #################################
        donor_bed, acceptor_bed = parse.create_donor_acceptor_bed(junction_bed, outdir)

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
        junction_fasta = prediction.splam_prediction(junction_fasta, junction_score_bed, splam_model, batch_size, device)

    elif args.subcommand == "clean":
        outdir = args.outdir
        threshold = args.threshold
        bad_intron_num = args.bad_intron_num

        # Running in clean gff file mode
        gff_db = outdir + "/annotation.db"
        if os.path.exists(gff_db):
            clean_gff.clean_gff(outdir, gff_db, threshold, bad_intron_num)
        else:
            argv_clean = sys.argv
            argv_clean.pop(0)
            argv_clean[0] = 'splam-clean'
            splam_clean.splam_clean(argv_clean)

if __name__ == "__main__":
    main()
