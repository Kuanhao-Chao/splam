import matplotlib.pyplot as plt
import sys
import numpy as np
import os
import pickle
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score


def print_topl_statistics(y_true, y_pred, topl_csv):
    # Prints the following information: top-kL statistics for k=0.5,1,2,4,
    # auprc, thresholds for k=0.5,1,2,4, number of true splice sites.
    idx_true = np.nonzero(y_true == 1)[0]
    argsorted_y_pred = np.argsort(y_pred)
    sorted_y_pred = np.sort(y_pred)
    topkl_accuracy = []
    threshold = []

    topl_fw = open(topl_csv, "w")
    for top_length in [0.5, 1, 2, 4, 10]:
        idx_pred = argsorted_y_pred[-int(top_length*len(idx_true)):]
        topkl_accuracy += [np.size(np.intersect1d(idx_true, idx_pred)) \
                  / float(min(len(idx_pred), len(idx_true)))]
        threshold += [sorted_y_pred[-int(top_length*len(idx_true))]] 

    auprc = average_precision_score(y_true, y_pred)
    print((("\n%.4f\t\033[91m%.4f\t\033[0m%.4f\t%.4f\t\033%.4f\t\033[94m%.4f\t\033[0m"
          + "%.4f\t%.4f\t%.4f\t%.4f\t%d\n\n") % (
          topkl_accuracy[0], topkl_accuracy[1], topkl_accuracy[2],
          topkl_accuracy[3], topkl_accuracy[4], auprc, threshold[0], threshold[1],
          threshold[2], threshold[3], len(idx_true))))
    topl_fw.write(str(topkl_accuracy[0]) + "\t" + str(topkl_accuracy[1]) + "\t" + str(topkl_accuracy[2]) + "\t" + 
                  str(topkl_accuracy[3]) + "\t" + str(topkl_accuracy[4]) + "\t" + str(auprc) + "\t" + 
                  str(threshold[0]) + "\t" + str(threshold[1]) + "\t" + str(threshold[2]) + "\t" + 
                  str(threshold[3]) + "\t" + str(len(idx_true)))
    topl_fw.close()


def main(argv):

    MODEL_VERSION = "SPLAM_v11"
    TYPE = "noshuffle"

    output_files = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
    for output_file in output_files:
        os.makedirs("splam_result/"+MODEL_VERSION+"/"+output_file+"/score_plts", exist_ok=True)
        print(">> output_file\t: ", output_file)

        # Declaring files on scores in the full region
        d_score_tsv_f = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/splam_all_seq.score.d."+TYPE+"."+output_file+".tsv"
        a_score_tsv_f = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/splam_all_seq.score.a."+TYPE+"."+output_file+".tsv"
        n_score_tsv_f = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/splam_all_seq.score.n."+TYPE+"."+output_file+".tsv"

        # Declaring files on scores on donor/acceptor sites
        d_score_ls_f = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/d_scores_"+TYPE+".tsv"
        a_score_ls_f = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/a_scores_"+TYPE+".tsv"

        d_score_fr = open(d_score_tsv_f, "r") 
        a_score_fr = open(a_score_tsv_f, "r") 
        n_score_fr = open(n_score_tsv_f, "r") 

        d_score_fw = open(d_score_ls_f, "w") 
        a_score_fw = open(a_score_ls_f, "w") 

        donor_p = np.zeros(400)
        acceptor_p = np.zeros(400)

        donor_scores = []
        acceptor_scores = []

        donor_concat = []
        acceptor_concat = []
        donor_l_concat = []
        acceptor_l_concat = []
        
        lines = d_score_fr.read().splitlines()    
        print("lines: ", len(lines))
        
        
        
        
        lines = lines[:10000]
        for line in lines:
            # scores = np.loadtxt(line)
            line = line.split(" ")
            donor_p_holder = np.array(line[0:400]).astype(float)
            donor_p = donor_p + donor_p_holder
            donor_ls = np.zeros(400)
            donor_ls[200] = 1

            if len(donor_concat) == 0:
                donor_concat = donor_p_holder
                donor_l_concat = donor_ls
            else:
                donor_concat =  np.concatenate((donor_concat, donor_p_holder), axis=0)
                donor_l_concat = np.concatenate((donor_l_concat, donor_ls))
            # Appending only score for the donor site
            donor_scores.append(donor_p_holder[200])
            d_score_fw.write(str(donor_p_holder[200]) + "\n")

        lines = a_score_fr.read().splitlines()
        lines = lines[:10000]
        for line in lines:
            # scores = np.loadtxt(line)
            line = line.split(" ")
            acceptor_p_holder = np.array(line[len(line)-400:len(line)]).astype(float)
            acceptor_p = acceptor_p + acceptor_p_holder
            acceptor_ls = np.zeros(400)
            acceptor_ls[200] = 1

            if len(acceptor_concat) == 0:
                acceptor_concat = acceptor_p_holder
                acceptor_l_concat = acceptor_ls
            else:
                acceptor_concat =  np.concatenate((acceptor_concat, acceptor_p_holder), axis=0)
                acceptor_l_concat = np.concatenate((acceptor_l_concat, acceptor_ls))
            # Appending only score for the donor site
            acceptor_scores.append(acceptor_p_holder[200])
            a_score_fw.write(str(acceptor_p_holder[200]) + "\n")


        print(donor_concat.shape)
        print(acceptor_concat.shape)
        print(donor_l_concat.shape)
        print(acceptor_l_concat.shape)

        #################################
        # Write out the list of scores of donors / acceptors
        #################################
        # d_score_fw.write("\t".join(donor_scores))
        # a_score_fw.write("\t".join(acceptor_scores))

        with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/d_scores_"+TYPE+".pkl", 'wb') as fw:
            pickle.dump(donor_scores, fw)
        with open("./splam_result/"+MODEL_VERSION+"/"+output_file+"/a_scores_"+TYPE+".pkl", 'wb') as fw:
            pickle.dump(acceptor_scores, fw)

        #########################
        # Calculate average donor / acceptor scores
        #########################
        avg_donor_p = donor_p/len(lines)
        avg_acceptor_p = acceptor_p/len(lines)

        scores = avg_donor_p, avg_acceptor_p
        my_dpi = 300
        plt.figure(figsize=(4200/my_dpi, 1200/my_dpi), dpi=my_dpi)
        donor = plt.plot(list(range(400)), avg_donor_p, color="blue", label="Donor")
        acceptor = plt.plot(list(range(400, 800)), avg_acceptor_p, color="red", label="Acceptor")
        plt.xlabel('Indices')
        plt.ylabel('Score')
        plt.legend()
        plt.tight_layout()
        plt.savefig("splam_result/"+MODEL_VERSION+"/"+output_file+"/score_plts/splam_"+output_file + ".png", dpi=my_dpi)
        # plt.show()

        #########################
        # Writing top-k statistics in file tsv file.
        #########################
        topl_csv = "./splam_result/"+MODEL_VERSION+"/"+output_file+"/topk_statistics.csv"

        print ("\n\033[1mDonor:\033[0m")
        print_topl_statistics(donor_l_concat,
        donor_concat, topl_csv)

        print ("\n\033[1mAcceptor:\033[0m")
        print_topl_statistics(acceptor_l_concat,
        acceptor_concat, topl_csv)

        donor_concat = donor_concat[donor_concat<0.1]
        acceptor_concat = acceptor_concat[acceptor_concat<0.1]
        
        # plt.violinplot([donor_concat, acceptor_concat])
        plt.show()

        d_score_fr.close()
        a_score_fr.close()
        n_score_fr.close()


if __name__ == "__main__":
    main(sys.argv[:])