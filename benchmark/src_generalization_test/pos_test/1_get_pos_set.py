# Parser for .gff NCBI files -> extract introns into .bed file

import os
import gffutils
from progress.bar import Bar

def run(db_name):
    
    # define output folder
    output_folder = './1_output/'
    input_folder = '../data/'
    input_filename = input_folder + db_name + '.gff'
    db_filename = output_folder + 'databases/' + db_name + '.db'
    output_filename = output_folder + db_name + '_introns.bed'

    # create database
    create_database(input_filename, db_filename)

    # connect to database
    db = gffutils.FeatureDB(db_filename)
    print(f'Successfully connected to {db_name} database')
    
    # write db file into bed format
    get_introns(db, output_filename)

    print(f'Complete.\nDatabase file at {db_filename}\nBED file at {output_filename}')


def create_database(input_filename, db_filename):

    # make output folder
    os.makedirs(os.path.dirname(db_filename), exist_ok=True)

    # get path of the input file
    fin = os.path.abspath(input_filename)

    # generate database if empty (~5 mins / 3.42 mil features), 
    if not os.path.exists(db_filename):
        db = gffutils.create_db(fin, db_filename, merge_strategy='create_unique', force=True, \
            disable_infer_transcripts=True, disable_infer_genes=True, verbose=True)


def get_introns(db, output_filename):
    print('Parsing file...')
    print(f'Feature types: {list(db.featuretypes())}')

    with open(output_filename, 'w') as bed_file:        
        # obtain list of all exons, then obtain a list of all parents of these exons
        pbar = Bar('Generating list of parents...', max=len(list(db.features_of_type('exon'))))
        parent_features = set()
        for exon in db.features_of_type('exon'):
            parents = db.parents(exon, level=1)
            parent_features.update(parent for parent in parents)
            pbar.next()
        pbar.finish()
        print(f'Found {len(parent_features)} parents.')

        # ordering the list of parents by seqid
        print('Sorting parent list...')
        parent_features = sorted(list(parent_features), key=lambda x: x.seqid)

        # parse through child exons of each parent
        pbar = Bar('Parsing children...', max=len(parent_features))
        cur_chr = ''
        for parent_feature in parent_features:
            # print current chromosome status
            if (cur_chr != parent_feature.seqid):
                cur_chr = parent_feature.seqid
                print(f'\tReading {cur_chr}')

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

                intron_chrom = ex1.seqid
                intron_start = ex1.end # BED is 0-indexed
                intron_end = ex2.start - 1 # BED is half-open interval
                name = parent_feature.id + '__' + str(count)
                
                # write intron to BED file
                bed_file.write(f'{intron_chrom}\t{intron_start}\t{intron_end}\t{name}\t{parent_feature.score}\t{parent_feature.strand}\n')

                count += 1
            
            pbar.next()
        pbar.finish()


if __name__ == '__main__':

    if os.getcwd() != 'pos_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_generalization_test/pos_test/')

    # define filenames
    annotation_files = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    file_idxs = [0,1,2,3] #CHANGEME
    
    for idx in file_idxs:       
        run(annotation_files[idx])