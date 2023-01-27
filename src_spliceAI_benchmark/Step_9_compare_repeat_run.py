import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve
import seaborn as sns

def main():

    splam_v2_j_shuffle_pred_prob = []
    splam_v2_j_shuffle_label_prob = []

    splam_v2_j_noshuffle_pred_prob = []
    splam_v2_j_noshuffle_label_prob = []

    splam_v2_j_nobatch_pred_prob = []
    splam_v2_j_nobatch_label_prob = []

    # TARGET = "batch"
    TARGET = "repeat"
    # for TYPE in ["shuffle", "noshuffle", "nobatch"]:
    stacked_array = ""
    prev_prbs = ""
    r_idx = 0
    if TARGET == "repeat":
        batch_size = 100
        for r_idx in range(0, 5):
            print(">> r_idx: ", r_idx)
            print("Reading './INPUT/splam.noshuffle."+str(r_idx)+"."+str(batch_size)+".pkl'")
            with open("./INPUT/splam.noshuffle."+str(r_idx)+"."+str(batch_size)+".pkl",'rb') as f:
                splam_v2_j_noshuffle_pred_prob = pickle.load(f)
                splam_v2_j_noshuffle_label_prob = pickle.load(f)

                splam_v2_j_noshuffle_pred_prob = splam_v2_j_noshuffle_pred_prob[:5000]
                splam_v2_j_noshuffle_label_prob = splam_v2_j_noshuffle_label_prob[:5000]


                print("splam_v2_j_noshuffle_pred_prob : ", splam_v2_j_noshuffle_pred_prob)
                # print("splam_v2_j_noshuffle_label_prob: ", splam_v2_j_noshuffle_label_prob)

                print("splam_v2_j_noshuffle_pred_prob : ", len(splam_v2_j_noshuffle_pred_prob))
                # print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

            splam_v2_j_noshuffle_label_prob = np.array(splam_v2_j_noshuffle_label_prob)
            splam_v2_j_noshuffle_pred_prob = np.array(splam_v2_j_noshuffle_pred_prob)


            splam_v2_j_noshuffle_label_prob = np.concatenate([splam_v2_j_noshuffle_label_prob, np.zeros(16)])
            splam_v2_j_noshuffle_pred_prob = np.concatenate([splam_v2_j_noshuffle_pred_prob, np.zeros(16)])

            if stacked_array == "":
                stacked_array = splam_v2_j_noshuffle_pred_prob
            else:
                stacked_array = np.vstack((stacked_array, splam_v2_j_noshuffle_pred_prob))

            print("Comparing to the previous one: ", splam_v2_j_noshuffle_pred_prob == prev_prbs)
            print("stacked_array", stacked_array)
            print("stacked_array.shape", stacked_array.shape)

            prev_prbs = splam_v2_j_noshuffle_pred_prob
            # print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

    elif TARGET == "batch":
        for batch_size in range(10, 110, 10):
            print(">> batch_size: ", batch_size)
            # print("Reading './INPUTS/SPLAM_v2/test.noshuffle.indices."+str(r_idx)+"."+str(batch_size)+".pkl'")
            # with open("./INPUTS/SPLAM_v2/test.noshuffle.indices."+str(r_idx)+"."+str(batch_size)+".pkl",'rb') as f:
            #     noshuffle_indices = pickle.load(f)
            #     print("indices: ", noshuffle_indices)
            #     print("indices: ", len(noshuffle_indices))

            print("Reading './INPUT/splam.noshuffle."+str(r_idx)+"."+str(batch_size)+".pkl'")
            with open("./INPUT/splam.noshuffle."+str(r_idx)+"."+str(batch_size)+".pkl",'rb') as f:
                splam_v2_j_noshuffle_pred_prob = pickle.load(f)
                splam_v2_j_noshuffle_label_prob = pickle.load(f)

                splam_v2_j_noshuffle_pred_prob = splam_v2_j_noshuffle_pred_prob[:5000]
                splam_v2_j_noshuffle_label_prob = splam_v2_j_noshuffle_label_prob[:5000]


                print("splam_v2_j_noshuffle_pred_prob : ", splam_v2_j_noshuffle_pred_prob)
                # print("splam_v2_j_noshuffle_label_prob: ", splam_v2_j_noshuffle_label_prob)

                print("splam_v2_j_noshuffle_pred_prob : ", len(splam_v2_j_noshuffle_pred_prob))
                # print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

            splam_v2_j_noshuffle_label_prob = np.array(splam_v2_j_noshuffle_label_prob)
            splam_v2_j_noshuffle_pred_prob = np.array(splam_v2_j_noshuffle_pred_prob)


            splam_v2_j_noshuffle_label_prob = np.concatenate([splam_v2_j_noshuffle_label_prob, np.zeros(16)])
            splam_v2_j_noshuffle_pred_prob = np.concatenate([splam_v2_j_noshuffle_pred_prob, np.zeros(16)])

            if stacked_array == "":
                stacked_array = splam_v2_j_noshuffle_pred_prob
            else:
                stacked_array = np.vstack((stacked_array, splam_v2_j_noshuffle_pred_prob))

            print("Comparing to the previous one: ", splam_v2_j_noshuffle_pred_prob == prev_prbs)
            print("stacked_array", stacked_array)
            print("stacked_array.shape", stacked_array.shape)

            prev_prbs = splam_v2_j_noshuffle_pred_prob
            # print("splam_v2_j_noshuffle_label_prob: ", len(splam_v2_j_noshuffle_label_prob))

    average_array = np.average(stacked_array, axis=0)
    print("stacked_array: ", stacked_array)
    print("stacked_array.shape: ", stacked_array.shape)
    print("average_array: ", average_array)
    print("average_array.shape: ", average_array.shape)

    if TARGET == "repeat":
        average_matrix = np.tile(average_array, (5, 1))
    elif TARGET == "batch":
        average_matrix = np.tile(average_array, (10, 1))
    print("average_matrix: ", average_matrix)
    print("average_matrix.shape: ", average_matrix.shape)

    heatmap_array = stacked_array - average_matrix
    print("heatmap_array: ", heatmap_array)

    # average_array = average_array.reshape(5016)
    # plt.imshow(heatmap_array, cmap='hot', interpolation='nearest', aspect="auto")
    s = sns.heatmap(heatmap_array, cmap="PiYG")
    sns.set(font_scale=1.4)
    # , annot=True)
    # , linewidth=0.5)
    s.figure.tight_layout()
    # s.figure.subplots_adjust(left = 0.1, top = 0.9) # change 0.3 to suit your needs.
    plt.xlabel("Samples [0-2999 positives;  3000-4999 negative]")
    if TARGET == "repeat":
        plt.ylabel("Repeat / Rerun (current result - average results)")
        plt.savefig("matrix_repeat.png", dpi=300)
    elif TARGET == "batch":
        plt.ylabel("Batch size (current result - average results)")
        plt.savefig("matrix_batch.png", dpi=300)
    plt.close()


        

    
    # print("len(splam_v2_j_shuffle_label_prob): ", len(splam_v2_j_shuffle_label_prob)) 

    # # label_diff = np.subtract(splam_v2_j_noshuffle_label_prob[shuffle_indices], splam_v2_j_shuffle_label_prob)

    # label_diff = np.subtract(splam_v2_j_noshuffle_label_prob[shuffle_indices], splam_v2_j_shuffle_label_prob)

    # # label_diff = splam_v2_j_noshuffle_label_prob[shuffle_indices] - splam_v2_j_shuffle_label_prob

    # pred_diff = splam_v2_j_noshuffle_pred_prob[shuffle_indices] - splam_v2_j_shuffle_pred_prob

    # plt.bar([*range(len(splam_v2_j_noshuffle_label_prob))], label_diff, width=1)
    # plt.savefig("label_diff.png")
    # plt.close()

    # plt.bar([*range(len(splam_v2_j_noshuffle_label_prob))], pred_diff, width=1)
    # plt.savefig("prob_diff.png")
    # plt.close()

if __name__ == "__main__":
    main()