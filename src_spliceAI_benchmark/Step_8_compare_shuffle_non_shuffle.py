import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():

    splam_v2_j_shuffle_pred_prob = []
    splam_v2_j_shuffle_label_prob = []

    splam_v2_j_noshuffle_pred_prob = []
    splam_v2_j_noshuffle_label_prob = []

    splam_v2_j_nobatch_pred_prob = []
    splam_v2_j_nobatch_label_prob = []


    # for TYPE in ["shuffle", "noshuffle", "nobatch"]:
    with open("./INPUTS/SPLAM_v2/test.shuffle.indices.pkl",'rb') as f:
        shuffle_indices = pickle.load(f)
        print("indices: ", shuffle_indices)
        print("indices: ", len(shuffle_indices))

    with open("./INPUTS/SPLAM_v2/test.noshuffle.indices.pkl",'rb') as f:
        noshuffle_indices = pickle.load(f)
        print("indices: ", noshuffle_indices)
        print("indices: ", len(noshuffle_indices))

    with open("./INPUTS/SPLAM_v2/test.nobatch.indices.pkl",'rb') as f:
        nobatch_indices = pickle.load(f)
        print("indices: ", nobatch_indices)
        print("indices: ", len(nobatch_indices))

    # shuffle_indices = np.array(shuffle_indices)    
    # noshuffle_indices = np.array(noshuffle_indices)    
    # nobatch_indices = np.array(nobatch_indices)    
    # shuffle_indices = shuffle_indices[:5900]
    # noshuffle_indices = noshuffle_indices[:5900]
    # nobatch_indices = nobatch_indices[:5900]

    with open("./INPUT/splam.shuffle.pkl",'rb') as f:
        splam_v2_j_shuffle_pred_prob = pickle.load(f)
        splam_v2_j_shuffle_label_prob = pickle.load(f)

        print("splam_v2_j_shuffle_pred_prob : ", splam_v2_j_shuffle_pred_prob)
        print("splam_v2_j_shuffle_label_prob: ", splam_v2_j_shuffle_label_prob)

        print("splam_v2_j_shuffle_pred_prob : ", len(splam_v2_j_shuffle_pred_prob))
        print("splam_v2_j_shuffle_label_prob: ", len(splam_v2_j_shuffle_label_prob))

    with open("./INPUT/splam.noshuffle.pkl",'rb') as f:
        splam_v2_j_noshuffle_pred_prob = pickle.load(f)
        splam_v2_j_noshuffle_label_prob = pickle.load(f)

        print("splam_v2_j_noshuffle_pred_prob : ", splam_v2_j_noshuffle_pred_prob)
        print("splam_v2_j_noshuffle_label_prob: ", splam_v2_j_noshuffle_label_prob)

        print("splam_v2_j_noshuffle_pred_prob : ", len(splam_v2_j_noshuffle_pred_prob))
        print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

    with open("./INPUT/splam.nobatch.pkl",'rb') as f:
        splam_v2_j_nobatch_pred_prob = pickle.load(f)
        splam_v2_j_nobatch_label_prob = pickle.load(f)

        print("splam_v2_j_nobatch_pred_prob : ", splam_v2_j_nobatch_pred_prob)
        print("splam_v2_j_nobatch_label_prob: ", splam_v2_j_nobatch_label_prob)

        print("splam_v2_j_nobatch_pred_prob : ", len(splam_v2_j_nobatch_pred_prob))
        print("splam_v2_j_nobatch_label_prob: ", len(splam_v2_j_nobatch_label_prob))

    splam_v2_j_shuffle_label_prob = np.array(splam_v2_j_shuffle_label_prob)
    splam_v2_j_shuffle_pred_prob = np.array(splam_v2_j_shuffle_pred_prob)

    splam_v2_j_shuffle_label_prob = np.concatenate([splam_v2_j_shuffle_label_prob, np.zeros(16)])
    splam_v2_j_shuffle_pred_prob = np.concatenate([splam_v2_j_shuffle_pred_prob, np.zeros(16)])

    # print("len(splam_v2_j_shuffle_label_prob): ", len(splam_v2_j_shuffle_label_prob))
    # print("len(splam_v2_j_shuffle_pred_prob): ", len(splam_v2_j_shuffle_pred_prob))


    splam_v2_j_noshuffle_label_prob = np.array(splam_v2_j_noshuffle_label_prob)
    splam_v2_j_noshuffle_pred_prob = np.array(splam_v2_j_noshuffle_pred_prob)


    splam_v2_j_noshuffle_label_prob = np.concatenate([splam_v2_j_noshuffle_label_prob, np.zeros(16)])
    splam_v2_j_noshuffle_pred_prob = np.concatenate([splam_v2_j_noshuffle_pred_prob, np.zeros(16)])

    splam_v2_j_nobatch_label_prob = np.array(splam_v2_j_nobatch_label_prob)
    splam_v2_j_nobatch_pred_prob = np.array(splam_v2_j_nobatch_pred_prob)



    print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

    print("shuffle_indices: ", len(shuffle_indices))
    print("splam_v2_j_noshuffle_label_prob[shuffle_indices]: ", len(splam_v2_j_noshuffle_label_prob[shuffle_indices]))

    print("len(splam_v2_j_shuffle_label_prob): ", len(splam_v2_j_shuffle_label_prob)) 

    # label_diff = np.subtract(splam_v2_j_noshuffle_label_prob[shuffle_indices], splam_v2_j_shuffle_label_prob)

    label_diff = np.subtract(splam_v2_j_noshuffle_label_prob[shuffle_indices], splam_v2_j_shuffle_label_prob)

    # label_diff = splam_v2_j_noshuffle_label_prob[shuffle_indices] - splam_v2_j_shuffle_label_prob

    pred_diff = splam_v2_j_noshuffle_pred_prob[shuffle_indices] - splam_v2_j_shuffle_pred_prob

    plt.bar([*range(len(splam_v2_j_noshuffle_label_prob))], label_diff, width=1)
    plt.savefig("label_diff.png")
    plt.close()

    plt.bar([*range(len(splam_v2_j_noshuffle_label_prob))], pred_diff, width=1)
    plt.savefig("prob_diff.png")
    plt.close()

if __name__ == "__main__":
    main()