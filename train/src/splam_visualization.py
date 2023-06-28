import os
import matplotlib.pyplot as plt

def main():
    LOG_DIR = "./MODEL/SPLAM_v11/LOG/"
    OUT_DIR = LOG_DIR + "IMG/"
    os.makedirs(OUT_DIR, exist_ok=True)

    ##########################################
    # Plotting train & test separately
    ##########################################
    for target in ["TRAIN", "VAL", "TEST"]:
        read_dir = LOG_DIR + target
        write_dir = OUT_DIR + target
        os.makedirs(write_dir, exist_ok=True)
        for (dirpath, dirnames, filenames) in os.walk(read_dir):
            # print(dirpath)
            # print(dirnames)
            # print(filenames)
            for filename in filenames:
                file = dirpath + '/' + filename
                outfile = write_dir + '/' + filename[:-3] + 'png'
                fr = open(file, 'r')
                lines = fr.read().splitlines()
                lines = [float(x) for x in lines]
                x = [*range(0, len(lines))]
                plt.plot(x, lines)
                plt.xlabel("Training batch num")
                plt.ylabel("Score")
                plt.title(filename[:-3])
                plt.savefig(outfile)
                plt.close()
                fr.close()

if __name__ == "__main__":
    main()
