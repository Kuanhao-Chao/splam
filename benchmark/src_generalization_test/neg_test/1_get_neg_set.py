# gets the gene loci start, end, and strand from all the src_neg

import os
import gffutils

def run(db_name):
    
    # define output folder
    output_folder = './1_output/'
    input_folder = '../data/'
    input_filename = input_folder + db_name + '.gff'
    db_filename = output_folder + 'databases/' + db_name + '.db'
    output_filename = output_folder + db_name + '_genes.bed'

    # create database
    create_database(input_filename, db_filename)

    # connect to database
    db = gffutils.FeatureDB(db_filename)
    print(f'Successfully connected to {db_name} database')
    
    # write db file into bed format
    get_genes(db, output_filename)
    
    print(f'Complete.\nDatabase file at {db_filename}\nBED file at {output_filename}')


def create_database(input_filename, db_filename):

    # make output folder
    os.makedirs(os.path.dirname(db_filename), exist_ok=True)

    # get path of the input file
    fin = os.path.abspath(input_filename)

    # generate database if empty (~5 mins / 3.42 mil features), 
    if not os.path.exists(db_filename):
        db = gffutils.create_db(fin, db_filename, merge_strategy="create_unique", force=True, \
            disable_infer_transcripts=True, disable_infer_genes=True, verbose=True)
        

def get_genes(db, output_filename):
    print('Parsing file...')
    print(f'Feature types: {list(db.featuretypes())}')

    with open(output_filename, 'w') as bed_file:
        i = 0
        for feature in db.features_of_type('gene'):
            if feature['gene_biotype'] != ['protein_coding']:
                continue
            i += 1
            try:
                bed_file.write(db.bed12(feature))
            except ValueError as e: # end of last exon does not match end of feature (or start)
                print(f'Exception on line {i}: {str(e)}')
                line = [feature.seqid, feature.start-1, feature.end, feature.id, feature.score, feature.strand]
                bed_file.write('\t'.join(map(str, line)))
            bed_file.write('\n')
    
    print(f'Total count: {i} protein-coding genes')


if __name__ == '__main__':

    if os.getcwd() != 'neg_test':
        os.chdir('/home/smao10/SPLAM/benchmark/src_generalization_test/neg_test/')

    # define filenames
    annotation_files = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
    file_idxs = [0,1,2,3] #CHANGEME
    
    for idx in file_idxs:       
        run(annotation_files[idx])