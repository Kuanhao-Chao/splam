# Parser for .gff3 file -> extract introns into .bed file
#   For standard gencode, chess3, and MANE files, using UCSC naming seqids
#   Every transcript contains exons, which this code leverages to only look at transcripts and the exons of those
#   Assumes that every exon is mutually exclusive of every other (no overlapping regions)

import os
import gffutils

def run_parser(input_filename):
    
    # define output folder
    output_folder = './Output/'
    input_folder = './Annotations/'
    title = input_filename.split('.')[0]
    input_filename = input_folder + input_filename
    db_filename = output_folder + 'databases/' + title + '.db'
    output_filename = output_folder + title + '_parsed.bed'

    # create database
    create_database(input_filename, db_filename, output_filename)

    # connect to database
    db = gffutils.FeatureDB(db_filename)
    print('Successfully connected to database')
    
    # write db file into bed format
    parse_database(db, output_filename)


def create_database(input_filename, db_filename, output_filename):

    # make output folder
    os.makedirs(os.path.dirname(db_filename), exist_ok=True)

    # get path of the input file
    fin = os.path.abspath(input_filename)

    # generate database if empty (~5 mins / 3.42 mil features)
    if not os.path.exists(db_filename):
        db = gffutils.create_db(fin, db_filename, merge_strategy="create_unique", force=True, \
                                disable_infer_transcripts=True, disable_infer_genes=True, verbose=True)


def handle_duplicate_names(path):
    # pre: path has '/path/to/file(optversion).ext' format
    filename, extension = os.path.splitext(path)
    count = 1
    while os.path.exists(path):
        path = filename + '(' + str(count) + ')' + extension 
        count += 1
    
    return path


def parse_database(db, output_filename):
    print('Parsing file...')
    print(f'Feature types: {list(db.featuretypes())}')

    # handling duplicate output files to avoid overwrite
    output_filename = handle_duplicate_names(output_filename)

    with open(output_filename, 'w') as bed_file:
        cur_chr = ''
        for transcript in db.features_of_type('transcript'):
            # print current chromosome status
            if (cur_chr != transcript.seqid):
                cur_chr = transcript.seqid
                print(f'Reading {cur_chr}')

            # obtain the child exons and sort them by start/end
            child_exons = db.children(transcript, featuretype='exon')
            sorted_exons = sorted(child_exons, key=lambda x: x.start)

            # iterate over consecutive exons to get introns in between
            count = 1
            for ex1, ex2 in zip(sorted_exons[:-1], sorted_exons[1:]):
                intron_chrom = ex1.seqid
                intron_start = ex1.end # BED is 0-indexed
                intron_end = ex2.start - 1 # BED is half-open interval
                name = transcript.id + '__' + str(count)
                
                # write intron to BED file
                bed_file.write(f'{intron_chrom}\t{intron_start}\t{intron_end}\t{name}\t{transcript.score}\t{transcript.strand}\n')

                count += 1

if __name__ == '__main__':
    # define filenames (VARS - can make this an I/O func)
    annotation_files = ['gencode.v43.annotation.gff3', 'refseq.genomic.gff', 'chess3.0.1.gff', 'MANE.GRCh38.v1.1.ensembl_genomic.gff']
    input_filename = 'gencode_all.v43.annotation.gff3'

    run_parser(input_filename)