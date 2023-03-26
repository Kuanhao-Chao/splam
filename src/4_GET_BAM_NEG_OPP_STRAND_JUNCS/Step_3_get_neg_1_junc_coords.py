import os

def main():
    files = ["./BAM_junctions/neg_hits.bed", "./BAM_junctions/pos_hits.bed"]
    ofile = "./BAM_junctions/junctions_1_cleaned.bed"
    fw = open(ofile, "w")

    for idx in range(len(files)):
        file = files[idx]


        with open(file, "r+") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                gtf_start = int(eles[1])
                gtf_end = int(eles[2])
                gtf_strand = eles[3]
                
                junc_start = int(eles[5])
                junc_end = int(eles[6])
                junc_strand = eles[9]
                eles_new = eles[4:]
                if ( gtf_start<junc_start and gtf_end>junc_end and gtf_strand!=junc_strand ):
                    fw.write("\t".join(eles_new)+"\n")
        
    fw.close()

if __name__ == "__main__":
    main()