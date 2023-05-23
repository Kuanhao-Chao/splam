import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import auc, accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
from sklearn import svm


# define a function to calculate precision, recall, F1 score, and queue rate
def calculate_metrics(threshold, true_labels, predict_probabilities):
    predictions = (predict_probabilities >= threshold).astype(int)
    true_positives = np.sum((predictions == 1) & (true_labels == 1))
    false_positives = np.sum((predictions == 1) & (true_labels == 0))
    false_negatives = np.sum((predictions == 0) & (true_labels == 1))
    true_negatives = np.sum((predictions == 0) & (true_labels == 0))
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    f1_score = 2 * (precision * recall) / (precision + recall)
    queue_rate = (true_positives + false_positives) / (true_positives + false_positives + true_negatives + false_negatives)
    return precision, recall, f1_score, queue_rate





def main():
    # for SPLAM_VERSION in ["SPLAM_v11", "SPLAM_v12"]:

    #####################################
    # Declaring parameters for probability & prediction array
    #####################################

    MANE_OR_ALTS = "ALTS"

    for SPLICEAI_VERSION in ["1", "2", "3", "4", "5"]:
        with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.N.merged.BOTH.pkl", "rb") as fr:
            spliceai_N_d_label = pickle.load(fr)
            spliceai_N_d_pred = pickle.load(fr)
            spliceai_N_a_label = pickle.load(fr)
            spliceai_N_a_pred = pickle.load(fr)

            # j_pred_prob = [x.numpy() for x in j_pred_prob]
            # j_pred_prob = [x.numpy() for x in j_pred_prob]

            print("\tspliceai_N_d_label : ", len(spliceai_N_d_label))
            print("\tspliceai_N_d_pred: ", len(spliceai_N_d_pred))
            print("\tspliceai_N_d_label : ", len(spliceai_N_d_label[spliceai_N_d_label == 1]))
            print("\tspliceai_N_d_pred: ", spliceai_N_d_pred)
            print("")
            print("\tspliceai_N_a_label : ", len(spliceai_N_a_label))
            print("\tspliceai_N_a_pred: ", len(spliceai_N_a_pred))
            print("\tspliceai_N_a_pred: ", len(spliceai_N_a_label[spliceai_N_a_label == 1]))
            print("\tspliceai_N_a_pred: ", spliceai_N_a_pred)
            print("")

        with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da.noN.merged.BOTH.pkl", "rb") as fr:
            spliceai_noN_d_label = pickle.load(fr)
            spliceai_noN_d_pred = pickle.load(fr)
            spliceai_noN_a_label = pickle.load(fr)
            spliceai_noN_a_pred = pickle.load(fr)

            print("\tspliceai_noN_d_label : ", len(spliceai_noN_d_label))
            print("\tspliceai_noN_d_pred: ", len(spliceai_noN_d_pred))
            print("\tspliceai_noN_d_pred: ", spliceai_noN_d_pred)
            print("")
            print("\tspliceai_noN_a_label : ", len(spliceai_noN_a_label))
            print("\tspliceai_noN_a_pred: ", len(spliceai_noN_a_pred))
            print("\tspliceai_noN_a_pred: ", spliceai_noN_a_pred)
            print("")
            



        for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:
            #####################################
            # Creating directories for visualization.
            #####################################
            os.makedirs("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/DT_plot/", exist_ok=True)

            predict_probabilities = np.minimum(spliceai_N_d_pred, spliceai_N_a_pred)
            true_labels = spliceai_N_a_label
            plot_DT_plot(true_labels, predict_probabilities, SPLAM_VERSION, "spliceai_N", SPLICEAI_VERSION)

            predict_probabilities = np.minimum(spliceai_noN_d_pred, spliceai_noN_a_pred)
            true_labels = spliceai_noN_a_label
            plot_DT_plot(true_labels, predict_probabilities, SPLAM_VERSION, "spliceai_noN", SPLICEAI_VERSION)

            with open("./splam_result/"+SPLAM_VERSION+"/splam.da.noshuffle.merged.BOTH.pkl",'rb') as f:
                splam_d_label = pickle.load(f)
                splam_d_pred = pickle.load(f)
                splam_a_label = pickle.load(f)
                splam_a_pred = pickle.load(f)
                print("\tsplam_d_label : ", len(splam_d_label))
                print("\tsplam_d_pred: ", len(splam_d_pred))
                print("")
                print("\tsplam_a_label : ", len(splam_a_label))
                print("\tsplam_a_pred: ", len(splam_a_pred))
                print("")

                # generate some sample data for predict probabilities and true labels
                predict_probabilities = np.minimum(splam_d_pred,splam_a_pred)
                true_labels = splam_d_label

                # predict_probabilities = np.minimum(spliceai_noN_d_pred, spliceai_noN_a_pred)
                # true_labels = spliceai_noN_a_label
                plot_DT_plot(true_labels, predict_probabilities, SPLAM_VERSION, "splam", SPLICEAI_VERSION)





def plot_DT_plot(true_labels, predict_probabilities, SPLAM_VERSION, target, SPLICEAI_VERSION):
    # define the range of threshold values to plot
    thresholds = np.arange(0, 0.99, 0.0001)

    # calculate the metrics for each threshold value
    precisions = []
    recalls = []
    f1_scores = []
    # queue_rates = []

    for threshold in thresholds:
        precision, recall, f1_score, queue_rate = calculate_metrics(threshold, true_labels, predict_probabilities)
        precisions.append(precision)
        recalls.append(recall)
        f1_scores.append(f1_score)
        # queue_rates.append(queue_rate)

    print(f1_scores)
    # find the optimal threshold based on F1 score
    optimal_index = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_index]
    optimal_precision = precisions[optimal_index]
    optimal_recall = recalls[optimal_index]
    optimal_f1_score = f1_scores[optimal_index]
    # optimal_queue_rate = queue_rates[optimal_index]
    print(f'Optimal threshold: {optimal_threshold}')
    print(f'Optimal_threshold: {optimal_threshold:.4f}, Precision: {optimal_precision:.2f}, Recall: {optimal_recall:.2f}, F1 score: {optimal_f1_score:.2f}')

    # plot the DT plot
    plt.figure(figsize=(9, 3.5))
    plt.plot(thresholds, precisions, label='Precision')
    plt.plot(thresholds, recalls, label='Recall')
    plt.plot(thresholds, f1_scores, label='F1 Score')
    # plt.plot(thresholds, queue_rates, label='Queue Rate')
    plt.axvline(x=optimal_threshold, linestyle='--', color='r', label='Optimal Threshold (maximum F1 score)')
    plt.xlabel('Threshold', size = 13)
    plt.ylabel('Performance Metrics', size = 13)
    plt.xlim(-0.02, 1.02)
    plt.ylim(0.0, 1.0)
    # plt.legend(loc='upper right', ncol=5, bbox_to_anchor=(1.0, 1.4))
    plt.title(f'Optimal_threshold: {optimal_threshold:.4f}, Precision: {optimal_precision:.2f}, Recall: {optimal_recall:.2f}, F1 score: {optimal_f1_score:.2f}')
    plt.xticks(size = 10)
    plt.yticks(size = 10)
    plt.tight_layout()
    # plt.show()
    plt.savefig("./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/DT_plot/"+target+".png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()