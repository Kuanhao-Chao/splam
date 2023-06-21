# Parser for .gff3/.gff file -> extract introns into .bed file

import os
import gffutils
global mappings

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
    output_filename = parse_database(db, output_filename)
    print(f'Complete.\nDatabase file at {db_filename}\nBED file at {output_filename}')


def get_mappings(format_num):
    global mappings
    mappings = {}

    with open('./GRCh38.p14_assembly_report.txt', 'r') as file:  
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            orig_name = columns[format_num]
            ucsc_name = columns[9]

            # store the key-value pair in the dictionary
            mappings[orig_name] = ucsc_name


def chr_transform(d):
    d['orig_seqid'] = d.seqid # preserves original seqid in case needed for data analysis
    d.seqid = mappings.get(d.seqid, d.seqid) # converts seqid to UCSC id (preserves if not found in conversion)
    return d


def create_database(input_filename, db_filename, output_filename):

    # make output folder
    os.makedirs(os.path.dirname(db_filename), exist_ok=True)

    # get path of the input file
    fin = os.path.abspath(input_filename)

    # generate database if empty (~5 mins / 3.42 mil features), 
    if not os.path.exists(db_filename):
        if 'gencode_all.v43.annotation.gff3' in input_filename:
            get_mappings(4) # genbank to ucsc
            db = gffutils.create_db(fin, db_filename, merge_strategy="create_unique", force=True, \
                disable_infer_transcripts=True, disable_infer_genes=True, transform=chr_transform, verbose=True)
        elif 'refseq.genomic.gff' in input_filename:
            get_mappings(6) # refseq to ucsc
            db = gffutils.create_db(fin, db_filename, merge_strategy="create_unique", force=True, \
                disable_infer_transcripts=True, disable_infer_genes=True, transform=chr_transform, verbose=True)
        else: 
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
        # obtain list of all exons, then obtain a list of all parents of these exons
        print('Generating list of parents...')
        parent_features = set()
        cur_parent = 0
        for exon in db.features_of_type('exon'):
            parents = db.parents(exon, level=1)
            parent_features.update(parent for parent in parents)
            if len(parent_features) % 10000 == 0 and len(parent_features) > cur_parent:
                print(f'Found {len(parent_features)} parents.')
                cur_parent = len(parent_features)
        print(f'Found {len(parent_features)} parents.')

        # ordering the list of parents by seqid
        print('Sorting parent list...')
        parent_features = sorted(list(parent_features), key=lambda x: x.seqid)

        # parse through child exons of each parent
        print('Parsing children...')
        cur_chr = ''
        for parent_feature in parent_features:
            # print current chromosome status
            if (cur_chr != parent_feature.seqid):
                cur_chr = parent_feature.seqid
                print(f'Reading {cur_chr}')

            # obtain the child exons and sort them by start/end
            child_exons = db.children(parent_feature, level=1, featuretype='exon')
            sorted_exons = sorted(child_exons, key=lambda x: (x.start, x.end))

            # iterate over consecutive exons to get introns in between
            count = 1
            for ex1, ex2 in zip(sorted_exons[:-1], sorted_exons[1:]):
                # filter out overlapping exons and invalid ids
                if int(ex1.end) > int(ex2.start):
                    print(f'ERROR: Exon {ex1.seqid}:{ex1.start}-{ex1.end} overlaps Exon {ex2.seqid}:{ex2.start}-{ex2.end}. Skipping first exon.')
                    continue #TODO: if this produces a lot of error, switch to while loop and skip second exon instead of first
                elif 'chr' not in ex1.seqid:
                    print(f'ERROR: {ex1.seqid} invalid UCSC-style ID. Skipping pair.')
                    continue

                intron_chrom = ex1.seqid
                intron_start = ex1.end # BED is 0-indexed
                intron_end = ex2.start - 1 # BED is half-open interval
                name = parent_feature.id + '__' + str(count)
                
                # write intron to BED file
                bed_file.write(f'{intron_chrom}\t{intron_start}\t{intron_end}\t{name}\t{parent_feature.score}\t{parent_feature.strand}\n')

                count += 1
        
        return output_filename


if __name__ == '__main__':
    # define filenames
    annotation_files = ['gencode.v43.annotation.gff3', 'gencode_all.v43.annotation.gff3', 'refseq.genomic.gff', 'chess3.0.1.gff', 'MANE.GRCh38.v1.1.ensembl_genomic.gff']
    file_idxs = [1] #CHANGEME
    
    for idx in file_idxs:       
        run_parser(annotation_files[idx])