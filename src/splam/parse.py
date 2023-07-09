import os
import pybedtools

from splam import config, chr_size

def create_donor_acceptor_bed(junction_bed, junction_dir):
    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    donor_bed = junction_dir+"/juncs/donor.bed"
    acceptor_bed = junction_dir+"/juncs/acceptor.bed"
    os.makedirs(junction_dir+"/juncs/", exist_ok=True)
    fw_donor = open(donor_bed, "w")
    fw_acceptor = open(acceptor_bed, "w")
    
    d_a_bed = junction_dir+"/juncs/d_a.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()

    with open(junction_bed, 'r') as f:
        lines = f.read().splitlines()
        style = ""
        if lines[0].split("\t")[0][:3] == "chr":
            chrs = chr_size.chrs_chr
            style = "'chr*' style"
        else:
            chrs = chr_size.chrs_refseq
            style = "'NCBI RefSeq' style"
        print(f'[Info] Chromosomes in the annotation file is in {style}')

        for line in lines:
            eles = line.split("\t")
            if len(eles) == 1:
                continue
            chr = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]

            if len(eles) == 7:
                trans = eles[6]
                junc_name = f'{junc_name};{score};{trans}'
            else:
                junc_name = f'{junc_name};{score};'

            if (strand == "+"):
                donor = int(eles[1])
                acceptor = int(eles[2])
                splice_junc_len = acceptor - donor
            elif (strand == "-"):
                acceptor = int(eles[1])
                donor = int(eles[2])
                splice_junc_len = donor - acceptor

            flanking_size = config.QUARTER_SEQ_LEN
            if splice_junc_len < config.QUARTER_SEQ_LEN:
                flanking_size = splice_junc_len

            if (strand == "+"):
                donor_s = donor - config.QUARTER_SEQ_LEN
                donor_e = donor + flanking_size
                acceptor_s = acceptor - flanking_size
                acceptor_e = acceptor + config.QUARTER_SEQ_LEN

            elif (strand == "-"):
                donor_s = donor - flanking_size
                donor_e = donor + config.QUARTER_SEQ_LEN
                acceptor_s = acceptor - config.QUARTER_SEQ_LEN
                acceptor_e = acceptor + flanking_size


            if donor_e >= chrs[chr] or acceptor_e >= chrs[chr]:
                continue
            if donor_s < 0 or acceptor_s < 0:
                continue
            new_junc = (chr, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
            if new_junc in JUNCS:
                continue
            else:
                JUNCS.add(new_junc)
                fw_donor.write(f'{chr}\t{str(donor_s)}\t{str(donor_e)}\t{junc_name}\t{score}\t{strand}\n')
                fw_acceptor.write(f'{chr}\t{str(acceptor_s)}\t{str(acceptor_e)}\t{junc_name}\t{score}\t{strand}\n')

                if (strand == "+"):
                    fw_da.write(f'{chr}\t{str(donor)}\t{str(acceptor+1)}\t{junc_name}\t{score}\t{strand}\n')
                elif (strand == "-"):
                    fw_da.write(f'{chr}\t{str(acceptor)}\t{str(donor+1)}\t{junc_name}\t{score}\t{strand}\n')

    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()
    return donor_bed, acceptor_bed

def write_donor_acceptor_fasta(bed_file, reference_genome):

    # Create a BedTool object from the BED file
    bed_f = pybedtools.BedTool(bed_file)

    # Get the sequences from the BED intervals
    sequences = bed_f.sequence(fi=reference_genome, s=True, name=True)

    fasta_f = os.path.splitext(bed_file)[0]+ ".fasta"
    b = sequences.save_seqs(fasta_f)
    return fasta_f

def concatenate_donor_acceptor_fasta(donor_fasta, acceptor_fasta, verbose):
    output_file = os.path.dirname(donor_fasta) + "/junction.fa"
    fw = open(output_file, "w")
    fr_donor = open(donor_fasta, "r")
    fr_acceptor = open(acceptor_fasta, "r")

    # parsing donor and acceptor fa files
    lines_d = fr_donor.read().splitlines()
    lines_a = fr_acceptor.read().splitlines()
    line_num = len(lines_d)

    # initializing stats
    canonical_d_count = 0 # GT
    noncanonical_d_count = 0
    canonical_a_count = 0 # AG
    noncanonical_a_count = 0
    donors = {}
    acceptors = {}
    num_skipped = 0
    for idx in range(0, line_num, 2):
        # PARSE FIRST LINE
        # >chr1:10000-20000(+)
        strand = lines_d[idx][-2]

        eles = lines_d[idx].split(":")
        d_name = eles[0]
        chromosome = eles[2]
        d_splits = eles[3].split("(")
        d_start, d_end = d_splits[0].split("-")

        eles = lines_a[idx].split(":")
        a_name = eles[0]
        chromosome = eles[2]
        a_splits = eles[3].split("(")
        a_start, a_end = a_splits[0].split("-")

        if strand == "+":
            donor_pos = int(d_start) + config.QUARTER_SEQ_LEN
            acceptor_pos = int(a_end) - config.QUARTER_SEQ_LEN
        elif strand == "-":
            donor_pos = int(d_end) - config.QUARTER_SEQ_LEN
            acceptor_pos = int(a_start) + config.QUARTER_SEQ_LEN

        # PARSE SECOND LINE
        idx += 1

        seq_d = lines_d[idx]
        seq_a = lines_a[idx]
        len_d = len(seq_d)
        len_a = len(seq_a)

        if len_d != len_a:
            print(f'Unequal lengths: seq_d {len_d}, seq_a {len_a}')

        if len_d == config.HALF_SEQ_LEN and len_a == config.HALF_SEQ_LEN:
            # combine normally
            x = seq_d + seq_a
        else:
            # pad with repeating N
            x = seq_d + (config.HALF_SEQ_LEN - len_d) * 'N' + (config.HALF_SEQ_LEN - len_a) * 'N' + seq_a

        x = x.upper()
        if (len(x) != int(config.SEQ_LEN)): 
            print("x: ", len(x))
       
        # skip sequence if there are Ns in the sequence
        if x[config.QUARTER_SEQ_LEN] == "N" or x[config.QUARTER_SEQ_LEN+1] == "N" or x[config.QUARTER_SEQ_LEN*3-2] == "N" or x[config.QUARTER_SEQ_LEN*3-1] == "N":
            num_skipped += 1
            continue
        
        # write the final fasta entry as two lines
        fw.write(f'>{chromosome};{str(donor_pos)};{str(acceptor_pos)};{strand};{d_name[1:]}\n')
        fw.write(x + "\n")

        # get stats on the dimers 
        donor_dimer = x[config.QUARTER_SEQ_LEN:config.QUARTER_SEQ_LEN+2]
        acceptor_dimer = x[config.QUARTER_SEQ_LEN*3-2:config.QUARTER_SEQ_LEN*3]

        if donor_dimer not in donors.keys():
            donors[donor_dimer] = 1
        else:
            donors[donor_dimer] += 1

        if acceptor_dimer not in acceptors.keys():
            acceptors[acceptor_dimer] = 1
        else:
            acceptors[acceptor_dimer] += 1

        if (donor_dimer == "GT"):
            canonical_d_count += 1
        else:
            noncanonical_d_count += 1
            
        if (acceptor_dimer == "AG"):
            canonical_a_count += 1
        else:
            noncanonical_a_count += 1

    # output stats
    if verbose:
        print("[Info] Loading splice site statistics ...")
        print("[Info] Canonical donor count: ", canonical_d_count)
        print("[Info] Noncanonical donor count: ", noncanonical_d_count)
        print("[Info] Canonical acceptor count: ", canonical_a_count)
        print("[Info] Noncanonical acceptor count: ", noncanonical_a_count)
        for key, value in donors.items():
            print("\tDonor   : ", key, " (", value, ")")
        for key, value in acceptors.items():
            print("\tAcceptor: ", key, " (", value, ")")

    fw.close()
    fr_acceptor.close()
    fr_donor.close()

    return output_file