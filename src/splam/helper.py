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
    print("chrs: ", chrs)
    return chrs

def get_chrom_size(path, type):
    chrs = {}
    with open(path, 'r') as file:       
        # skip header
        next(file)

        # read the file line by line
        for line in file:  
            # split by tabs
            columns = line.strip().split('\t')
            if type == "refseq":
                refseq_name = columns[6]
            elif type == "chr":
                refseq_name = columns[9]

            chrom_size = int(columns[8])

            # store the key-value pair in the dictionary
            chrs[refseq_name] = chrom_size
    print("chrs: ", chrs)
    return chrs

get_chrom_size("GRCh38.p14_assembly_report.txt", "refseq")
get_chrom_size("GRCh38.p14_assembly_report.txt", "chr")