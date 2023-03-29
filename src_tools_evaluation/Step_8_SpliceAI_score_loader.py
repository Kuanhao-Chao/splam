import matplotlib.pyplot as plt
import sys
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main(argv):


    TYPE = argv[1]
    output_file = argv[2]
    # output_files = "outlier_test"
    # paths = ('./models/spliceai{}.h5'.format(x) for x in range(1, 2))
    # models = [load_model(resource_filename('spliceai', x)) for x in paths]
    path = './models/spliceai1.h5'
    print(">> path\t\t: ", path)
    print(">> output_file\t: ", output_file)

    label = '.'
    if output_file == "pos" or output_file == "pos_MANE" or output_file == "pos_ALTS":
        label = '+'

    if output_file == "neg_1" or output_file == "neg_random":
        label = '-'

    # print(">> label\t\t: ", label)
    da_faf = "./dataset/"+output_file+"/spliceai/spliceai."+TYPE+".juncs.seq.fa"


    d_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
    a_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
    n_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"

    name_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.name."+TYPE+"."+output_file+".tsv"

    d_score_fw = open(d_score_tsv_f, "r") 
    a_score_fw = open(a_score_tsv_f, "r") 
    n_score_fw = open(n_score_tsv_f, "r") 

    name_fw = open(name_tsv_f, "r") 


    # scores = np.loadtxt(score_tsv_f)

    lines = d_score_fw.read().splitlines()

    donor_p = np.zeros(400)
    acceptor_p = np.zeros(400)

    
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        # print(len(line))

        donor_ps = np.array(line[0:400]).astype(float)
        acceptor_ps = np.array(line[len(line)-400:len(line)]).astype(float)

        donor_p = donor_p + donor_ps
        acceptor_p = acceptor_p + acceptor_ps

        # print(len(donor_ps))
        # print(len(acceptor_ps))

    lines = a_score_fw.read().splitlines()
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        acceptor_ps = np.array(line[len(line)-400:len(line)]).astype(float)
        acceptor_p = acceptor_p + acceptor_ps




    # print(len(donor_p))
    # print(len(acceptor_p))
    avg_donor_p = donor_p/len(lines)
    avg_acceptor_p = acceptor_p/len(lines)
    scores = np.concatenate((avg_donor_p, avg_acceptor_p))
    plt.plot(list(range(800)), scores)

    plt.show()


    # lines = a_score_fw.read().splitlines()
    # for line in lines:
    #     # scores = np.loadtxt(line)
    #     line = line.split(" ")
    #     print(len(line))

    # lines = n_score_fw.read().splitlines()
    # for line in lines:
    #     # scores = np.loadtxt(line)
    #     line = line.split(" ")
    #     print(len(line))

    d_score_fw.close()
    a_score_fw.close()
    n_score_fw.close()



if __name__ == "__main__":
    main(sys.argv[:])