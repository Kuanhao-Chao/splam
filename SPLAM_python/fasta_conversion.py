# converts refseq accession ids to ucsc chromosome names in .fa file

def run_converter():
    # retrieve the dict for converting between ids
    mappings = get_mappings('./GRCh38.p14_assembly_report.txt')

    # convert ids in the fasta
    convert_fasta('../Dataset/GRCh38_p14_genomic.fa', mappings, '../Dataset/GRCh38_p14_genomic_ucsc.fa')


def get_mappings(conversion_file):

    mappings = {}

    with open(conversion_file, 'r') as file:       
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
    
    return mappings


def convert_fasta(fasta_file, mappings, output_filename):
    # open the .fa file in read mode and the output file in write mode
    with open(fasta_file, 'r') as input_file, open(output_filename, 'w') as output_file:
        ucsc_name = ''
        # read the input file line by line
        for line in input_file:
            if line.startswith('>'):
                # extract the RefSeq accession ID from the header line
                refseq_accn = line.split(' ', 1)[0].lstrip('>')
                print(f'RefSeq ID: {refseq_accn}')
                # convert the RefSeq accession ID to UCSC-style name using the dictionary
                ucsc_name = mappings.get(refseq_accn, 'unknown')
                # replace the RefSeq accession ID with the UCSC-style name in the header line
                modified_line = f'>{ucsc_name} {line.split(" ", 1)[1]}'
                print(f'MODIFIED: {modified_line}')
                output_file.write(modified_line)
            else:
                output_file.write(line)


if __name__ == '__main__':
    run_converter()