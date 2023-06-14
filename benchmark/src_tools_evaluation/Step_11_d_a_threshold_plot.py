import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay

def plot_scatter_plot(label_d, score_d, label_a, score_a, filename):
    junc_prob = label_d.astype(bool)
    non_junc_prob = (1-label_d).astype(bool)
    # print("\tspliceai_junc_prob: ", junc_prob)
    # print("\tspliceai_junc_prob: ", non_junc_prob)

    fig, ax = plt.subplots()
    ax.set_title(filename)
    ax.set_xlabel("Donor score")
    ax.set_ylabel("Acceptor score")

    junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 0.3)
    non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)    
    ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])

    fig.savefig(filename)
    # fig.close()


def plot_min_linear_plot(label_d, score_d, label_a, score_a, filename):
    junction_score = np.minimum(score_d, score_a)

    junc_prob = label_d.astype(bool)
    non_junc_prob = (1-label_d).astype(bool)
    print("\tspliceai_junc_prob: ", junc_prob)
    print("\tspliceai_junc_prob: ", non_junc_prob)

    fig, ax = plt.subplots()
    ax.set_title(filename)
    # ax.set_xlabel("Donor score")
    # ax.set_ylabel("Acceptor score")

    ar = np.arange(10) # just as an example array
    junc_legend = ax.plot(junction_score[junc_prob], np.zeros_like(junction_score[junc_prob]) + 0.)
                        #   , s = 0.3)
    non_junc_legend = ax.plot(junction_score[non_junc_prob], np.zeros_like(junction_score[non_junc_prob]) + 0.)
                            #   , s = 0.3)
    # junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 0.3)
    # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 0.3)
    ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])
    fig.savefig(filename)
    

# THRESHOLDS = [0.1, 0.01]



# # Helper function used for visualization in the following examples
# def identify_axes(ax_dict, fontsize=48):
#     """
#     Helper to identify the Axes in the examples below.

#     Draws the label in a large font in the center of the Axes.

#     Parameters
#     ----------
#     ax_dict : dict[str, Axes]
#         Mapping between the title / label and the Axes.
#     fontsize : int, optional
#         How big the label should be.
#     """
#     kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
#     for k, ax in ax_dict.items():
#         ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)

def main():

    MANE_OR_ALTS = "ALTS"
    threshold_min = 0.0001
    threshold_max = 0.9999
    THRESHOLDS= np.arange(threshold_min, threshold_max, 0.001)


    splam_threshold = 0.0707
    for TARGET in ["noN", "N"]:
        if TARGET == "noN":
            spliceai_threshold = 0.0003
        elif TARGET == "N":
            spliceai_threshold = 0.0070

        for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:

            for x_axis_rep in ["log", "no_log"]:
                # spliceai_TP__splam_TP_idices
                spliceai_TP__splam_FN = []
                spliceai_FN__splam_TP = []
                spliceai_FN__splam_FN = []

                # spliceai_TN__splam_TN = []
                spliceai_TN__splam_FP = []
                spliceai_FP__splam_TN = []
                spliceai_FP__splam_FP = []
                for threshold in THRESHOLDS:
                    print(">> threshold: ", threshold)
                    #####################################
                    # Creating directories for visualization.
                    #####################################
                    # figure_root = "./IMG/d_a/"

                    #####################################
                    # Declaring parameters for probability & prediction array
                    #####################################
                    spliceai_N_d_pred_prob = []
                    spliceai_N_d_label_prob = []
                    spliceai_N_a_pred_prob = []
                    spliceai_N_a_label_prob = []

                    spliceai_d_pred = []
                    spliceai_d_label = []
                    spliceai_a_pred = []
                    spliceai_a_label = []


                    splam_S_d_pred_prob = []
                    splam_S_d_label_prob = []
                    splam_S_a_pred_prob = []
                    splam_S_a_label_prob = []

                    splam_d_pred = []
                    splam_d_pred = []
                    splam_a_pred = []
                    splam_noS_a_label_prob = []






                    #####################################
                    # Creating directories for visualization.
                    #####################################
                    os.makedirs("IMG/" + SPLAM_VERSION + "/TP_TN_FP_FN/"+TARGET+"/", exist_ok=True)

                    figure_root = "./IMG/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/scatter_plot/"
                    target_figure_root = figure_root+TARGET+"/"
                    # os.makedirs(target_figure_root+"tsv_"+str(threshold)+"/", exist_ok=True)
                    

                    with open("./spliceai_result_"+SPLICEAI_VERSION+"/spliceai.da."+TARGET+".merged.BOTH.pkl", "rb") as fr:
                        spliceai_d_label = pickle.load(fr)
                        spliceai_d_pred = pickle.load(fr)
                        spliceai_a_label = pickle.load(fr)
                        spliceai_a_pred = pickle.load(fr)

                        print("\tspliceai_d_label : ", len(spliceai_d_label))
                        print("\tspliceai_d_pred: ", len(spliceai_d_pred))
                        print("\tspliceai_d_pred: ", spliceai_d_pred)
                        print("")
                        print("\tspliceai_a_label : ", len(spliceai_a_label))
                        print("\tspliceai_a_pred: ", len(spliceai_a_pred))
                        print("\tspliceai_a_pred: ", spliceai_a_pred)
                        print("")


                    #     # spliceai_d_pred = [x.numpy() for x in spliceai_d_pred]
                    #     # spliceai_a_pred = [x.numpy() for x in spliceai_a_pred]
                    # spliceai_d_pred = np.array(spliceai_d_pred)
                    # spliceai_a_pred = np.array(spliceai_a_pred)
                    # spliceai_d_label = np.array(spliceai_d_label)
                    # spliceai_a_label = np.array(spliceai_a_label)
                    # spliceai_d_label = np.array(spliceai_d_label)

                    # print("spliceai_d_pred : ", spliceai_d_pred)
                    # print("spliceai_d_label: ", spliceai_d_label)
                    # print("spliceai_a_pred : ", spliceai_a_pred)
                    # print("spliceai_a_label: ", spliceai_a_label)
                    # print("spliceai_d_label: ", spliceai_d_label)

                    # print("spliceai_d_pred : ", len(spliceai_d_pred))
                    # print("spliceai_d_label: ", len(spliceai_d_label))
                    # print("spliceai_a_pred : ", len(spliceai_a_pred))
                    # print("spliceai_a_label: ", len(spliceai_a_label))

                    # print("spliceai_d_label: ", len(spliceai_d_label))


                    # with open("./splam_result/splam.da.shuffle.pkl",'rb') as f:
                    #     splam_S_d_label_prob = pickle.load(f)
                    #     splam_S_d_pred_prob = pickle.load(f)
                        
                    #     splam_S_a_label_prob = pickle.load(f)
                    #     splam_S_a_pred_prob = pickle.load(f)

                    #     print("splam_S_d_label_prob : ", splam_S_d_label_prob)
                    #     print("splam_S_d_pred_prob: ", splam_S_d_pred_prob)

                    #     print("splam_S_a_label_prob : ", splam_S_a_label_prob)
                    #     print("splam_S_a_pred_prob: ", splam_S_a_pred_prob)

                    

                    with open("./splam_result/"+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/splam.da.noshuffle.merged.BOTH.pkl",'rb') as f:
                        splam_d_label = pickle.load(f)
                        splam_d_pred = pickle.load(f)
                        splam_a_label = pickle.load(f)
                        splam_a_pred = pickle.load(f)
                        print("\tsplam_d_pred : ", len(splam_d_pred))
                        print("\tsplam_d_pred: ", len(splam_d_pred))
                        print("")
                        print("\tsplam_a_label : ", len(splam_a_label))
                        print("\tsplam_a_pred: ", len(splam_a_pred))
                        print("")




                    # # print(spliceai_d_pred >= threshold)
                    # spliceai_N_TP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 1)
                    # spliceai_N_FN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 1)
                    # spliceai_N_FP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 0)
                    # spliceai_N_TN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 0)

                    # print("spliceai_TP_idices : ", len(spliceai_N_junc_name[spliceai_N_TP_idices]))
                    # print("spliceai_FP_idices : ", len(spliceai_N_junc_name[spliceai_N_FP_idices]))
                    # print("spliceai_FN_idices : ", len(spliceai_N_junc_name[spliceai_N_FN_idices]))
                    # print("spliceai_TN_idices : ", len(spliceai_N_junc_name[spliceai_N_TN_idices]))


                    # print(spliceai_d_pred >= threshold)
                    spliceai_TP_idices = (spliceai_d_pred >= threshold) & (spliceai_a_pred >= threshold) & (spliceai_d_label == 1)
                    spliceai_FN_idices = ((spliceai_d_pred < threshold) | (spliceai_a_pred < threshold)) & (spliceai_d_label == 1)
                    spliceai_FP_idices = (spliceai_d_pred >= threshold) & (spliceai_a_pred >= threshold) & (spliceai_d_label == 0)
                    spliceai_TN_idices = ((spliceai_d_pred < threshold) | (spliceai_a_pred < threshold)) & (spliceai_d_label == 0)

                    print("spliceai_TP_idices : ", len(spliceai_d_label[spliceai_TP_idices]))
                    print("spliceai_FP_idices : ", len(spliceai_d_label[spliceai_FP_idices]))
                    print("spliceai_FN_idices : ", len(spliceai_d_label[spliceai_FN_idices]))
                    print("spliceai_TN_idices : ", len(spliceai_d_label[spliceai_TN_idices]))


                    splam_TP_idices = (splam_d_pred >= threshold) & (splam_a_pred >= threshold) & (splam_d_label== 1)
                    splam_FN_idices = ((splam_d_pred < threshold) | (splam_a_pred < threshold)) & (splam_d_label == 1)
                    splam_FP_idices = (splam_d_pred >= threshold) & (splam_a_pred >= threshold) & (splam_d_label == 0)
                    splam_TN_idices = ((splam_d_pred < threshold) | (splam_a_pred < threshold)) & (splam_d_label == 0)

                    print("splam_TP_idices : ", len(splam_d_pred[splam_TP_idices]))
                    print("splam_FP_idices : ", len(splam_d_pred[splam_FP_idices]))
                    print("splam_FN_idices : ", len(splam_d_pred[splam_FN_idices]))
                    print("splam_TN_idices : ", len(splam_d_pred[splam_TN_idices]))

                    spliceai_TP__splam_TP_idices = spliceai_TP_idices & splam_TP_idices
                    spliceai_TP__splam_FN_idices = spliceai_TP_idices & splam_FN_idices
                    spliceai_FN__splam_TP_idices = spliceai_FN_idices & splam_TP_idices
                    spliceai_FN__splam_FN_idices = spliceai_FN_idices & splam_FN_idices

                    spliceai_TN__splam_TN_idices = spliceai_TN_idices & splam_TN_idices
                    spliceai_TN__splam_FP_idices = spliceai_TN_idices & splam_FP_idices
                    spliceai_FP__splam_TN_idices = spliceai_FP_idices & splam_TN_idices
                    spliceai_FP__splam_FP_idices = spliceai_FP_idices & splam_FP_idices

                    print("spliceai_TP__splam_FN_idices : ", len(splam_d_pred[spliceai_TP__splam_FN_idices]))
                    print("spliceai_FN__splam_TP_idices : ", len(splam_d_pred[spliceai_FN__splam_TP_idices]))
                    print("spliceai_FN__splam_FN_idices : ", len(splam_d_pred[spliceai_FN__splam_FN_idices]))

                    print("spliceai_TP__splam_FN_idices : ", len(splam_d_pred[spliceai_TN__splam_FP_idices]))
                    print("spliceai_FN__splam_TP_idices : ", len(splam_d_pred[spliceai_FP__splam_TN_idices]))
                    print("spliceai_FN__splam_FN_idices : ", len(splam_d_pred[spliceai_FP__splam_FP_idices]))

                    spliceai_TP__splam_FN.append(len(splam_d_pred[spliceai_TP__splam_FN_idices]))
                    spliceai_FN__splam_TP.append(len(splam_d_pred[spliceai_FN__splam_TP_idices]))
                    spliceai_FN__splam_FN.append(len(splam_d_pred[spliceai_FN__splam_FN_idices]))

                    spliceai_TN__splam_FP.append(len(splam_d_pred[spliceai_TN__splam_FP_idices]))
                    spliceai_FP__splam_TN.append(len(splam_d_pred[spliceai_FP__splam_TN_idices]))
                    spliceai_FP__splam_FP.append(len(splam_d_pred[spliceai_FP__splam_FP_idices]))

                
                if x_axis_rep == "log":
                    THRESHOLDS_PLT = np.array(THRESHOLDS)
                    THRESHOLDS_PLT = np.log(THRESHOLDS_PLT*1000)
                    print("THRESHOLDS_PLT: ", THRESHOLDS_PLT)
                else:
                    THRESHOLDS_PLT = THRESHOLDS
                plt.figure(figsize=(12, 3))
                plt.xlabel("Threshold ("+str(threshold_min)+" - "+ str(threshold_max) +")")
                plt.ylabel("Number of splice sites")

                plt.plot(THRESHOLDS_PLT, spliceai_FN__splam_TP, label="SPLAM_TP__SpliceAI_FN")
                # plt.plot(THRESHOLDS_PLT, spliceai_FP__splam_TN, label="SPLAM_TN__SpliceAI_FP")

                plt.plot(THRESHOLDS_PLT, spliceai_TP__splam_FN, label="SPLAM_FN__SpliceAI_TP")
                # plt.plot(THRESHOLDS_PLT, spliceai_TN__splam_FP, label="SPLAM_FP__SpliceAI_TN")

                # plt.plot(THRESHOLDS_PLT, spliceai_FN__splam_FN, label="SPLAM_FN__SpliceAI_FN")
                # plt.plot(THRESHOLDS_PLT, spliceai_FP__splam_FP, label="SPLAM_FP__SpliceAI_FP")
                
                plt.axvline(x = splam_threshold, linestyle ='--', color = 'r', label = 'SPLAM_threshold ('+str(splam_threshold)+')')
                plt.axvline(x = spliceai_threshold, linestyle ='--', color = 'b', label = 'SpliceAI_threshold('+str(spliceai_threshold)+')')

                plt.legend(loc="center right", bbox_to_anchor=(1.3, 0.5))
                plt.tight_layout()
                plt.savefig("IMG/" + SPLAM_VERSION + "/TP_TN_FP_FN/"+TARGET+"/threshold_"+x_axis_rep+"_"+str(threshold_min)+"_"+str(threshold_max)+".png", dpi=300)
                # plt.show()
                plt.close()


if __name__ == "__main__":
    main()


