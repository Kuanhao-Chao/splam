import matplotlib.pyplot as plt
import sys
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score



def print_topl_statistics(y_true, y_pred):
    # Prints the following information: top-kL statistics for k=0.5,1,2,4,
    # auprc, thresholds for k=0.5,1,2,4, number of true splice sites.
    # print ("y_true: ", y_true.shape)
    # print ("y_pred: ", y_pred.shape)

    idx_true = np.nonzero(y_true == 1)[0]
    # print(("idx_true: ", idx_true))
    argsorted_y_pred = np.argsort(y_pred)
    # print(("argsorted_y_pred: ", argsorted_y_pred))
    # print(("argsorted_y_pred.shape: ", argsorted_y_pred.shape))
    sorted_y_pred = np.sort(y_pred)
    # print(("sorted_y_pred: ", sorted_y_pred))
    # print(("sorted_y_pred.shape: ", sorted_y_pred.shape))

    topkl_accuracy = []
    threshold = []

    for top_length in [0.5, 1, 2, 4, 10]:

        idx_pred = argsorted_y_pred[-int(top_length*len(idx_true)):]
        # print(("idx_pred: ", idx_pred))
        
        # print(("np.size(np.intersect1d(idx_true, idx_pred)): ", np.size(np.intersect1d(idx_true, idx_pred))))
        # print(("float(min(len(idx_pred), len(idx_true))): ", float(min(len(idx_pred), len(idx_true)))))
        topkl_accuracy += [np.size(np.intersect1d(idx_true, idx_pred)) \
                  / float(min(len(idx_pred), len(idx_true)))]
        threshold += [sorted_y_pred[-int(top_length*len(idx_true))]] 

    auprc = average_precision_score(y_true, y_pred)

    # print((("%.4f\t\033[91m%.4f\t\033[0m%.4f\t%.4f\t\033%.4f\t\033[94m%.4f\t\033[0m"
    #       + "\t%d") % (
    #       topkl_accuracy[0], topkl_accuracy[1], topkl_accuracy[2],
    #       topkl_accuracy[3], topkl_accuracy[4], auprc,  len(idx_true))))

    print((("\n%.4f\t\033[91m%.4f\t\033[0m%.4f\t%.4f\t\033%.4f\t\033[94m%.4f\t\033[0m"
          + "%.4f\t%.4f\t%.4f\t%.4f\t%d\n\n") % (
          topkl_accuracy[0], topkl_accuracy[1], topkl_accuracy[2],
          topkl_accuracy[3], topkl_accuracy[4], auprc, threshold[0], threshold[1],
          threshold[2], threshold[3], len(idx_true))))
    



def main(argv):

    SUBSET = 400
    os.makedirs("score_plts", exist_ok=True)
    SPLAM_TYPE = argv[1]
    SPLICEAI_TYPE = argv[2]

    output_file = argv[3]
    # output_files = "outlier_test"
    # paths = ('./models/splam{}.h5'.format(x) for x in range(1, 2))
    # models = [load_model(resource_filename('splam', x)) for x in paths]
    print(">> output_file\t: ", output_file)

    label = '.'
    if output_file == "pos" or output_file == "pos_MANE" or output_file == "pos_ALTS":
        label = '+'

    if output_file == "neg_1" or output_file == "neg_random":
        label = '-'

    # print(">> label\t\t: ", label)
    da_faf = "./dataset/"+output_file+"/splam/splam."+SPLAM_TYPE+".juncs.seq.fa"


    d_score_tsv_f = "./splam_result/"+output_file+"/splam_all_seq.score.d."+SPLAM_TYPE+"."+output_file+".tsv"
    a_score_tsv_f = "./splam_result/"+output_file+"/splam_all_seq.score.a."+SPLAM_TYPE+"."+output_file+".tsv"
    n_score_tsv_f = "./splam_result/"+output_file+"/splam_all_seq.score.n."+SPLAM_TYPE+"."+output_file+".tsv"

    name_tsv_f = "./splam_result/"+output_file+"/splam_all_seq.name."+SPLAM_TYPE+"."+output_file+".tsv"

    d_score_fw = open(d_score_tsv_f, "r") 
    a_score_fw = open(a_score_tsv_f, "r") 
    n_score_fw = open(n_score_tsv_f, "r") 

    name_fw = open(name_tsv_f, "r") 


    # scores = np.loadtxt(score_tsv_f)


    donor_p = np.zeros(400)
    acceptor_p = np.zeros(400)

    donor_concat = []
    acceptor_concat = []
    donor_l_concat = []
    acceptor_l_concat = []
    
    lines = d_score_fw.read().splitlines()    
    lines = lines[:SUBSET]
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        donor_ps = np.array(line[0:400]).astype(float)
        donor_p = donor_p + donor_ps
        donor_ls = np.zeros(400)
        donor_ls[200] = 1

        if len(donor_concat) == 0:
            donor_concat = donor_ps
            donor_l_concat = donor_ls
        else:
            donor_concat =  np.concatenate((donor_concat, donor_ps), axis=0)
            donor_l_concat = np.concatenate((donor_l_concat, donor_ls))


        # print(len(donor_ps))
        # print(len(acceptor_ps))

    lines = a_score_fw.read().splitlines()
    lines = lines[:SUBSET]
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        acceptor_ps = np.array(line[len(line)-400:len(line)]).astype(float)
        acceptor_p = acceptor_p + acceptor_ps
        acceptor_ls = np.zeros(400)
        acceptor_ls[200] = 1

        if len(acceptor_concat) == 0:
            acceptor_concat = acceptor_ps
            acceptor_l_concat = acceptor_ls
        else:
            acceptor_concat =  np.concatenate((acceptor_concat, acceptor_ps), axis=0)
            acceptor_l_concat = np.concatenate((acceptor_l_concat, acceptor_ls))


    print(donor_concat.shape)
    print(acceptor_concat.shape)
    print(donor_l_concat.shape)
    print(acceptor_l_concat.shape)
    

    print(len(donor_p))
    print(len(acceptor_p))
    avg_donor_p = donor_p/len(lines)
    avg_acceptor_p = acceptor_p/len(lines)
    # scores = np.concatenate((avg_donor_p, avg_acceptor_p))

    # avg_donor_p = np.pad(avg_donor_p,(0, 400), 'constant')
    # avg_acceptor_p = np.pad(avg_acceptor_p,(400, 0), 'constant')

    scores = avg_donor_p, avg_acceptor_p

    my_dpi = 300
    plt.figure(figsize=(4200/my_dpi, 1200/my_dpi), dpi=my_dpi)
    donor = plt.plot(list(range(400)), avg_donor_p, color="blue", label="SPLAM Donor")
    acceptor = plt.plot(list(range(400, 800)), avg_acceptor_p, color="red", label="SPLAM Acceptor")
















    # print(">> label\t\t: ", label)
    da_faf = "./dataset/"+output_file+"/spliceai/spliceai."+SPLICEAI_TYPE+".juncs.seq.fa"


    d_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.d."+SPLICEAI_TYPE+"."+output_file+".tsv"
    a_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.a."+SPLICEAI_TYPE+"."+output_file+".tsv"
    n_score_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.score.n."+SPLICEAI_TYPE+"."+output_file+".tsv"

    name_tsv_f = "./spliceai_result/"+output_file+"/spliceai_all_seq.name."+SPLICEAI_TYPE+"."+output_file+".tsv"

    d_score_fw = open(d_score_tsv_f, "r") 
    a_score_fw = open(a_score_tsv_f, "r") 
    n_score_fw = open(n_score_tsv_f, "r") 

    name_fw = open(name_tsv_f, "r") 


    # scores = np.loadtxt(score_tsv_f)


    donor_p = np.zeros(400)
    acceptor_p = np.zeros(400)

    donor_concat = []
    acceptor_concat = []
    donor_l_concat = []
    acceptor_l_concat = []
    
    lines = d_score_fw.read().splitlines()    
    lines = lines[:SUBSET]
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        donor_ps = np.array(line[0:400]).astype(float)
        donor_p = donor_p + donor_ps
        donor_ls = np.zeros(400)
        donor_ls[199] = 1

        if len(donor_concat) == 0:
            donor_concat = donor_ps
            donor_l_concat = donor_ls
        else:
            donor_concat =  np.concatenate((donor_concat, donor_ps), axis=0)
            donor_l_concat = np.concatenate((donor_l_concat, donor_ls))


        # print(len(donor_ps))
        # print(len(acceptor_ps))

    lines = a_score_fw.read().splitlines()
    lines = lines[:SUBSET]
    for line in lines:
        # scores = np.loadtxt(line)
        line = line.split(" ")
        acceptor_ps = np.array(line[len(line)-400:len(line)]).astype(float)
        acceptor_p = acceptor_p + acceptor_ps
        acceptor_ls = np.zeros(400)
        acceptor_ls[200] = 1

        if len(acceptor_concat) == 0:
            acceptor_concat = acceptor_ps
            acceptor_l_concat = acceptor_ls
        else:
            acceptor_concat =  np.concatenate((acceptor_concat, acceptor_ps), axis=0)
            acceptor_l_concat = np.concatenate((acceptor_l_concat, acceptor_ls))


    print(donor_concat.shape)
    print(acceptor_concat.shape)
    print(donor_l_concat.shape)
    print(acceptor_l_concat.shape)
    

    print(len(donor_p))
    print(len(acceptor_p))
    avg_donor_p = -donor_p/len(lines)
    avg_acceptor_p = -acceptor_p/len(lines)
    # scores = np.concatenate((avg_donor_p, avg_acceptor_p))

    # avg_donor_p = np.pad(avg_donor_p,(0, 400), 'constant')
    # avg_acceptor_p = np.pad(avg_acceptor_p,(400, 0), 'constant')

    scores = avg_donor_p, avg_acceptor_p

    my_dpi = 300
    # plt.figure(figsize=(4200/my_dpi, 1200/my_dpi), dpi=my_dpi)
    donor = plt.plot(list(range(400)), avg_donor_p, color="purple", label="SpliceAI Donor")
    acceptor = plt.plot(list(range(400, 800)), avg_acceptor_p, color="orange", label="SpliceAI Acceptor")



    plt.xlabel('Indices')
    plt.ylabel('Score')
    # plt.ylim((-1, 1))
    plt.legend()
    plt.tight_layout()
    plt.savefig("score_plts/splam_" + SPLAM_TYPE + "_spliceai_" + SPLICEAI_TYPE + "_" + output_file + ".png", dpi=my_dpi)
    # plt.show()


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

    print ("\n\033[1mDonor:\033[0m")
    print_topl_statistics(donor_l_concat,
    donor_concat)

    print ("\n\033[1mAcceptor:\033[0m")
    print_topl_statistics(acceptor_l_concat,
    acceptor_concat)

    donor_concat = donor_concat[donor_concat<0.1]
    acceptor_concat = acceptor_concat[acceptor_concat<0.1]
    
    # plt.violinplot([donor_concat, acceptor_concat])
    plt.show()

    d_score_fw.close()
    a_score_fw.close()
    n_score_fw.close()



if __name__ == "__main__":
    main(sys.argv[:])