import pyfaidx

# use pyfaidx, fix the bed12 indexing...


'''Obtain chromosome size'''
def get_chrom_size(path):
    chrs = {}
    with open(path, 'r') as file:       
        # read the file line by line
        for line in file:  
            if line.startswith('#'):
                continue
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[6]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    
    return chrs


