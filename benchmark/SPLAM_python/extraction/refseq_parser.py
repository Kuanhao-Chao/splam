# Parser for .gff3/.gff file -> extract introns into .bed file
#   For RefSeq genomes which need conversion from RefSeq accession IDs to UCSC
#   Also takes into account that exons may be overlapping in region, in which case will exclude those
#   Also considers that not every exon may have a 'transcript'-labeled parent, so will collect parents
#   of every exon before iterating over exon children of parents, at the cost of runtime

import os
import gffutils
global mappings

def run_parser(input_filename):
    
    # define output folder
    output_folder = './Output/'
    input_folder = './Annotations/'
    title = input_filename.split('.')[0]
    input_filename = input_folder + input_filename
    db_filename = output_folder + 'databases/refseq_filtered.db'
    output_filename = output_folder + title + '_parsed.bed'

    # create database
    create_database(input_filename, db_filename, output_filename)

    # connect to database
    db = gffutils.FeatureDB(db_filename)
    print('Successfully connected to database')
    
    # write db file into bed format
    parse_database(db, output_filename)

def get_mappings(file):
    global mappings
    mappings = {}

    with open('./GRCh38.p14_assembly_report.txt', 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            refseq_accn = columns[6]
            ucsc_name = columns[9]

            # store the key-value pair in the dictionary
            mappings[refseq_accn] = ucsc_name


def chr_transform(d):
    d.seqid = mappings[d.seqid]
    return d


def create_database(input_filename, db_filename, output_filename):

    # make output folder
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    # get path of the input file
    fin = os.path.abspath(input_filename)

    # get the refseq to ucsc mappings
    get_mappings('GRCh38.p14_assembly_report.txt')

    # generate database if empty (~5 mins / 3.42 mil features)
    if not os.path.exists(db_filename):
        db = gffutils.create_db(fin, db_filename, merge_strategy="create_unique", force=True, \
                                disable_infer_transcripts=True, disable_infer_genes=True, transform=chr_transform, verbose=True)


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
        # obtain list of all exons, then obtain a list of all parents of these exons
        print('Generating list of parents...')
        parent_ids = set()
        for exon in db.features_of_type('exon'):
            parents = db.parents(exon, level=1)
            parent_ids.update(parent.attributes['ID'][0] for parent in parents)
        print(f'Found {len(parent_ids)} parents.')

        print('Parsing children...')
        cur_chr = ''
        for parent_id in parent_ids:   
            # print current chromosome status
            if (cur_chr != exon.seqid):
                cur_chr = exon.seqid
                print(f'Reading {cur_chr}')

            # obtain the child exons and sort them by start/end
            parent_feature = db[parent_id]
            child_exons = db.children(parent_feature, level=1, featuretype='exon')
            sorted_exons = sorted(child_exons, key=lambda x: (x.start, x.end))

            # iterate over consecutive exons to get introns in between
            count = 1
            for ex1, ex2 in zip(sorted_exons[:-1], sorted_exons[1:]):
                # filter out overlapping exons
                if int(ex1.end) > int(ex2.start):
                    continue

                intron_chrom = ex1.seqid
                intron_start = ex1.end # BED is 0-indexed
                intron_end = ex2.start - 1 # BED is half-open interval
                name = parent_feature.id + '__' + str(count)
                
                # write intron to BED file
                bed_file.write(f'{intron_chrom}\t{intron_start}\t{intron_end}\t{name}\t{parent_feature.score}\t{parent_feature.strand}\n')

                count += 1

if __name__ == '__main__':
    # define filenames (VARS - can make this an I/O func)
    input_filename = 'refseq.genomic.gff'
    run_parser(input_filename)