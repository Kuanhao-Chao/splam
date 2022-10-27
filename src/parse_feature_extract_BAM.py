import os
import subprocess

def main():
    files = os.listdir("../Dataset/tranx_label_raw/")
    for file in files:
        if file[0:4] == "pos_":
            # print(file)
            f = open("../Dataset/tranx_label_raw/"+file, 'r')
            line = f.readline()
            line = line.split("\t")
            region = line[0]+":"+line[1]+"-"+line[2]
            print(region)

            output_bam = "../Dataset/tranx_bam/" + file[4:-4]+ ".bam"
            subprocess.run(["samtools", "view", "../Dataset/SRR1352415.bam", region, "-o", output_bam ])

# samtools view input.bam "Chr10:18000-45500" > output.bam




            


if __name__ == "__main__":
    main()