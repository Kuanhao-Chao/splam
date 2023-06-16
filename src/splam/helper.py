def get_hg38_chrom_size():
    
    chrs = {}
    with open('GRCh38.p14_assembly_report.txt', 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            ucsc_name = columns[9]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[ucsc_name] = chrom_size
    
    return chrs

def get_chrom_size(path):
    chrs = {}
    with open(path, 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            refseq_name = columns[9]
            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    return chrs