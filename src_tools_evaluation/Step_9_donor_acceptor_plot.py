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


    junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 1)
    non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 1)    
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
                        #   , s = 1)
    non_junc_legend = ax.plot(junction_score[non_junc_prob], np.zeros_like(junction_score[non_junc_prob]) + 0.)
                            #   , s = 1)
    # junc_legend = ax.scatter(score_d[junc_prob], score_a[junc_prob], s = 1)
    # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 1)
    ax.legend([junc_legend, non_junc_legend], ['Junction', 'Non Junction'])
    fig.savefig(filename)
    

THRESHOLDS = [0.1, 0.01]

def main():

    for threshold in THRESHOLDS:
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

        spliceai_d_pred_prob = []
        spliceai_d_label_prob = []
        spliceai_a_pred_prob = []
        spliceai_a_label_prob = []


        splam_S_d_pred_prob = []
        splam_S_d_label_prob = []
        splam_S_a_pred_prob = []
        splam_S_a_label_prob = []

        splam_noS_d_pred_prob = []
        splam_noS_d_label_prob = []
        splam_noS_a_pred_prob = []
        splam_noS_a_label_prob = []

        for TARGET in ["noN", "N"]:
            figure_root = "./IMG/TP_FN_FP_TN/"
            target_figure_root = figure_root+TARGET+"/"
            os.makedirs(target_figure_root+"tsv_"+str(threshold)+"/", exist_ok=True)
            # with open("./spliceai_result/spliceai.N.merged.pkl",'rb') as f:
            #     spliceai_N_d_pred_prob = pickle.load(f)
            #     spliceai_N_d_label_prob = pickle.load(f)
            #     spliceai_N_a_pred_prob = pickle.load(f)
            #     spliceai_N_a_label_prob = pickle.load(f)
            #     spliceai_N_junc_name = pickle.load(f)

            # spliceai_N_d_pred_prob = np.array(spliceai_N_d_pred_prob)
            # spliceai_N_a_pred_prob = np.array(spliceai_N_a_pred_prob)
            # spliceai_N_d_label_prob = np.array(spliceai_N_d_label_prob)
            # spliceai_N_a_label_prob = np.array(spliceai_N_a_label_prob)
            # spliceai_N_junc_name = np.array(spliceai_N_junc_name)

            # print("spliceai_N_d_pred_prob : ", spliceai_N_d_pred_prob)
            # print("spliceai_N_d_label_prob: ", spliceai_N_d_label_prob)
            # print("spliceai_N_a_pred_prob : ", spliceai_N_a_pred_prob)
            # print("spliceai_N_a_label_prob: ", spliceai_N_a_label_prob)
            # print("spliceai_N_junc_name: ", spliceai_N_junc_name)

            # print("spliceai_N_d_pred_prob : ", len(spliceai_N_d_pred_prob))
            # print("spliceai_N_d_label_prob: ", len(spliceai_N_d_label_prob))
            # print("spliceai_N_a_pred_prob : ", len(spliceai_N_a_pred_prob))
            # print("spliceai_N_a_label_prob: ", len(spliceai_N_a_label_prob))

            # print("spliceai_N_junc_name: ", len(spliceai_N_junc_name))



            with open("./spliceai_result/spliceai."+TARGET+".merged.pkl",'rb') as f:
                spliceai_d_pred_prob = pickle.load(f)
                spliceai_d_label_prob = pickle.load(f)
                spliceai_a_pred_prob = pickle.load(f)
                spliceai_a_label_prob = pickle.load(f)
                spliceai_junc_name = pickle.load(f)

                # spliceai_d_pred_prob = [x.numpy() for x in spliceai_d_pred_prob]
                # spliceai_a_pred_prob = [x.numpy() for x in spliceai_a_pred_prob]
            spliceai_d_pred_prob = np.array(spliceai_d_pred_prob)
            spliceai_a_pred_prob = np.array(spliceai_a_pred_prob)
            spliceai_d_label_prob = np.array(spliceai_d_label_prob)
            spliceai_a_label_prob = np.array(spliceai_a_label_prob)
            spliceai_junc_name = np.array(spliceai_junc_name)

            print("spliceai_d_pred_prob : ", spliceai_d_pred_prob)
            print("spliceai_d_label_prob: ", spliceai_d_label_prob)
            print("spliceai_a_pred_prob : ", spliceai_a_pred_prob)
            print("spliceai_a_label_prob: ", spliceai_a_label_prob)
            print("spliceai_junc_name: ", spliceai_junc_name)

            print("spliceai_d_pred_prob : ", len(spliceai_d_pred_prob))
            print("spliceai_d_label_prob: ", len(spliceai_d_label_prob))
            print("spliceai_a_pred_prob : ", len(spliceai_a_pred_prob))
            print("spliceai_a_label_prob: ", len(spliceai_a_label_prob))

            print("spliceai_junc_name: ", len(spliceai_junc_name))


            # with open("./splam_result/splam.da.shuffle.pkl",'rb') as f:
            #     splam_S_d_label_prob = pickle.load(f)
            #     splam_S_d_pred_prob = pickle.load(f)
                
            #     splam_S_a_label_prob = pickle.load(f)
            #     splam_S_a_pred_prob = pickle.load(f)

            #     print("splam_S_d_label_prob : ", splam_S_d_label_prob)
            #     print("splam_S_d_pred_prob: ", splam_S_d_pred_prob)

            #     print("splam_S_a_label_prob : ", splam_S_a_label_prob)
            #     print("splam_S_a_pred_prob: ", splam_S_a_pred_prob)

            
            with open("./splam_result/splam.da.noshuffle.merged.pkl",'rb') as f:
                splam_noS_d_label_prob = pickle.load(f)
                splam_noS_d_pred_prob = pickle.load(f)
                
                splam_noS_a_label_prob = pickle.load(f)
                splam_noS_a_pred_prob = pickle.load(f)
                splam_noS_junc_name = pickle.load(f)

                print("splam_noS_d_label_prob : ", splam_noS_d_label_prob)
                print("splam_noS_d_pred_prob: ", splam_noS_d_pred_prob)

                print("splam_noS_a_label_prob : ", splam_noS_a_label_prob)
                print("splam_noS_a_pred_prob: ", splam_noS_a_pred_prob)
                print("splam_noS_junc_name  : ", splam_noS_junc_name)

            # with open("./spliceai_result/splam.nobatch.0.100.pkl",'rb') as f:
            #     splam_j_nobatch_pred_prob = pickle.load(f)
            #     splam_j_nobatch_label_prob = pickle.load(f)
            #     print("splam_j_nobatch_pred_prob : ", splam_j_nobatch_pred_prob)
            #     print("splam_j_nobatch_label_prob: ", splam_j_nobatch_label_prob)
            #     print("splam_j_nobatch_pred_prob : ", len(splam_j_nobatch_pred_prob))
            #     print("splam_j_nobatch_label_prob: ", len(splam_j_nobatch_label_prob))

            # splam_S_d_label_prob = np.array(splam_S_d_label_prob)
            # splam_S_d_pred_prob = np.array(splam_S_d_pred_prob)
            # splam_S_a_label_prob = np.array(splam_S_a_label_prob)
            # splam_S_a_pred_prob = np.array(splam_S_a_pred_prob)

            splam_noS_d_label_prob = np.array(splam_noS_d_label_prob)
            splam_noS_d_pred_prob = np.array(splam_noS_d_pred_prob)
            splam_noS_a_label_prob = np.array(splam_noS_a_label_prob)
            splam_noS_a_pred_prob = np.array(splam_noS_a_pred_prob)
            splam_noS_junc_name = np.array(splam_noS_junc_name)


            # # print(spliceai_d_pred_prob >= threshold)
            # spliceai_N_TP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 1)
            # spliceai_N_FN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 1)
            # spliceai_N_FP_idices = (spliceai_N_d_pred_prob >= threshold) & (spliceai_N_a_pred_prob >= threshold) & (spliceai_N_d_label_prob == 0)
            # spliceai_N_TN_idices = ((spliceai_N_d_pred_prob < threshold) | (spliceai_N_a_pred_prob < threshold)) & (spliceai_N_d_label_prob == 0)

            # print("spliceai_TP_idices : ", len(spliceai_N_junc_name[spliceai_N_TP_idices]))
            # print("spliceai_FP_idices : ", len(spliceai_N_junc_name[spliceai_N_FP_idices]))
            # print("spliceai_FN_idices : ", len(spliceai_N_junc_name[spliceai_N_FN_idices]))
            # print("spliceai_TN_idices : ", len(spliceai_N_junc_name[spliceai_N_TN_idices]))


            # print(spliceai_d_pred_prob >= threshold)
            spliceai_TP_idices = (spliceai_d_pred_prob >= threshold) & (spliceai_a_pred_prob >= threshold) & (spliceai_d_label_prob == 1)
            spliceai_FN_idices = ((spliceai_d_pred_prob < threshold) | (spliceai_a_pred_prob < threshold)) & (spliceai_d_label_prob == 1)
            spliceai_FP_idices = (spliceai_d_pred_prob >= threshold) & (spliceai_a_pred_prob >= threshold) & (spliceai_d_label_prob == 0)
            spliceai_TN_idices = ((spliceai_d_pred_prob < threshold) | (spliceai_a_pred_prob < threshold)) & (spliceai_d_label_prob == 0)

            print("spliceai_TP_idices : ", len(spliceai_junc_name[spliceai_TP_idices]))
            print("spliceai_FP_idices : ", len(spliceai_junc_name[spliceai_FP_idices]))
            print("spliceai_FN_idices : ", len(spliceai_junc_name[spliceai_FN_idices]))
            print("spliceai_TN_idices : ", len(spliceai_junc_name[spliceai_TN_idices]))


            splam_TP_idices = (splam_noS_d_pred_prob >= threshold) & (splam_noS_a_pred_prob >= threshold) & (splam_noS_d_label_prob == 1)
            splam_FN_idices = ((splam_noS_d_pred_prob < threshold) | (splam_noS_a_pred_prob < threshold)) & (splam_noS_d_label_prob == 1)
            splam_FP_idices = (splam_noS_d_pred_prob >= threshold) & (splam_noS_a_pred_prob >= threshold) & (splam_noS_d_label_prob == 0)
            splam_TN_idices = ((splam_noS_d_pred_prob < threshold) | (splam_noS_a_pred_prob < threshold)) & (splam_noS_d_label_prob == 0)

            print("splam_TP_idices : ", len(splam_noS_junc_name[splam_TP_idices]))
            print("splam_FP_idices : ", len(splam_noS_junc_name[splam_FP_idices]))
            print("splam_FN_idices : ", len(splam_noS_junc_name[splam_FN_idices]))
            print("splam_TN_idices : ", len(splam_noS_junc_name[splam_TN_idices]))

            spliceai_TP__splam_TP_idices = spliceai_TP_idices & splam_TP_idices
            spliceai_TP__splam_FN_idices = spliceai_TP_idices & splam_FN_idices
            spliceai_FN__splam_TP_idices = spliceai_FN_idices & splam_TP_idices
            spliceai_FN__splam_FN_idices = spliceai_FN_idices & splam_FN_idices

            spliceai_TN__splam_TN_idices = spliceai_TN_idices & splam_TN_idices
            spliceai_TN__splam_FP_idices = spliceai_TN_idices & splam_FP_idices
            spliceai_FP__splam_TN_idices = spliceai_FP_idices & splam_TN_idices
            spliceai_FP__splam_FP_idices = spliceai_FP_idices & splam_FP_idices

            # print("spliceai_TP__splam_FN_idices : ", len(splam_noS_junc_name[spliceai_TP__splam_FN_idices]))
            # print("spliceai_FN__splam_TP_idices : ", len(splam_noS_junc_name[spliceai_FN__splam_TP_idices]))
            # print("spliceai_FN__splam_FN_idices : ", len(splam_noS_junc_name[spliceai_FN__splam_FN_idices]))

            # print("spliceai_TP__splam_FN_idices : ", len(splam_noS_junc_name[spliceai_TN__splam_FP_idices]))
            # print("spliceai_FN__splam_TP_idices : ", len(splam_noS_junc_name[spliceai_FP__splam_TN_idices]))
            # print("spliceai_FN__splam_FN_idices : ", len(splam_noS_junc_name[spliceai_FP__splam_FP_idices]))

            # junc_prob = label_d.astype(bool)
            # non_junc_prob = (1-label_d).astype(bool)
            # # print("\tspliceai_junc_prob: ", junc_prob)
            # # print("\tspliceai_junc_prob: ", non_junc_prob)


            for TYPE in ["TP_TP", "TP_FN", "FN_TP", "FN_FN", "TN_TN", "TN_FP", "FP_TN", "FP_FP"]:
                if TYPE == "TP_TP":
                    indices = spliceai_TP__splam_TP_idices
                elif TYPE == "TP_FN":
                    indices = spliceai_TP__splam_FN_idices
                elif TYPE == "FN_TP":
                    indices = spliceai_FN__splam_TP_idices
                elif TYPE == "FN_FN":
                    indices = spliceai_FN__splam_FN_idices
                elif TYPE == "TN_TN":
                    indices = spliceai_TN__splam_TN_idices
                elif TYPE == "TN_FP":
                    indices = spliceai_TN__splam_FP_idices
                elif TYPE == "FP_TN":
                    indices = spliceai_FP__splam_TN_idices
                elif TYPE == "FP_FP":
                    indices = spliceai_FP__splam_FP_idices

                fw = open(target_figure_root+"tsv_"+str(threshold)+"/"+TYPE+"_junc_"+str(threshold)+".tsv", "w")
                for junc in splam_noS_junc_name[indices]:
                    fw.write('\t'.join(junc)+"\n")
                fw.close()




            for TOOL in ["spliceai", "splam"]:
                if TOOL == "spliceai":
                    d_score = spliceai_d_pred_prob
                    a_score = spliceai_a_pred_prob
                elif TOOL == "splam":
                    d_score = splam_noS_d_pred_prob
                    a_score = splam_noS_a_pred_prob
                fig, ax = plt.subplots()
                ax.set_title(TOOL + " score")
                ax.set_xlabel("Donor score")
                ax.set_ylabel("Acceptor score")

                # junc_legend = 

                ax.hlines(y=threshold, xmin=0, xmax=1, linewidth=1, color='r', linestyles="dashed")
                ax.vlines(x=threshold, ymin=0, ymax=1, linewidth=1, color='r', linestyles="dashed")


                splcieai_TP__splam_TP = ax.scatter(d_score[spliceai_TP__splam_TP_idices], a_score[spliceai_TP__splam_TP_idices], s = 0, color="pink")
                splcieai_TP__splam_TP_len = len(d_score[spliceai_TP__splam_TP_idices])

                splcieai_TP__splam_FN = ax.scatter(d_score[spliceai_TP__splam_FN_idices], a_score[spliceai_TP__splam_FN_idices], s = 6, color="blue")
                splcieai_TP__splam_FN_len = len(d_score[spliceai_TP__splam_FN_idices])

                splcieai_FN__splam_TP = ax.scatter(d_score[spliceai_FN__splam_TP_idices], a_score[spliceai_FN__splam_TP_idices], s = 6, color="red")
                spliceai_FN__splam_TP_len = len(d_score[spliceai_FN__splam_TP_idices])

                splcieai_FN__splam_FN = ax.scatter(d_score[spliceai_FN__splam_FN_idices], a_score[spliceai_FN__splam_FN_idices], s = 6, color="green")
                spliceai_FN__splam_FN_len = len(d_score[spliceai_FN__splam_FN_idices])


                splcieai_TN__splam_TN = ax.scatter(d_score[spliceai_TN__splam_TN_idices], a_score[spliceai_TN__splam_TN_idices], s = 0, color="gray")
                splcieai_TN__splam_TN_len = len(d_score[spliceai_TN__splam_TN_idices])

                splcieai_TN__splam_FP = ax.scatter(d_score[spliceai_TN__splam_FP_idices], a_score[spliceai_TN__splam_FP_idices], s = 6, color="orange")
                spliceai_TN__splam_FP_len = len(d_score[spliceai_TN__splam_FP_idices])

                splcieai_FP__splam_TN = ax.scatter(d_score[spliceai_FP__splam_TN_idices], a_score[spliceai_FP__splam_TN_idices], s = 6, color="cyan")
                spliceai_FP__splam_TN_len = len(d_score[spliceai_FP__splam_TN_idices])

                splcieai_FP__splam_FP = ax.scatter(d_score[spliceai_FP__splam_FP_idices], a_score[spliceai_FP__splam_FP_idices], s = 6, color="purple")
                spliceai_FP__splam_FP_len = len(d_score[spliceai_FP__splam_FP_idices])

                # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 1)
                lgd = ax.legend([splcieai_TP__splam_TP, splcieai_TP__splam_FN, splcieai_FN__splam_TP, splcieai_FN__splam_FN, splcieai_TN__splam_TN, splcieai_TN__splam_FP, splcieai_FP__splam_TN, splcieai_FP__splam_FP], ["splcieai_TP__splam_TP ("+str(splcieai_TP__splam_TP_len)+")", "splcieai_TP__splam_FN ("+str(splcieai_TP__splam_FN_len)+")", "splcieai_FN__splam_TP ("+str(spliceai_FN__splam_TP_len)+")", "splcieai_FN__splam_FN ("+str(spliceai_FN__splam_FN_len)+")", "splcieai_TN__splam_TN ("+str(splcieai_TN__splam_TN_len)+")", "splcieai_TN__splam_FP ("+str(spliceai_TN__splam_FP_len)+")", "splcieai_FP__splam_TN ("+str(spliceai_FP__splam_TN_len)+")", "splcieai_FP__splam_FP ("+str(spliceai_FP__splam_FP_len)+")"], loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlim([0, 1])
                plt.ylim([0, 1])
                ax.set_aspect('equal', adjustable='box')
                plt.savefig(target_figure_root + TOOL + "_" + str(threshold) + ".png", dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close()


            # splcieai_TP__splam_FN = ax.scatter(splam_noS_d_pred_prob[spliceai_TP__splam_FN_idices], splam_noS_a_pred_prob[spliceai_TP__splam_FN_idices], s = 5, color="blue")
            # splcieai_FN__splam_TP = ax.scatter(splam_noS_d_pred_prob[spliceai_FN__splam_TP_idices], splam_noS_a_pred_prob[spliceai_FN__splam_TP_idices], s = 5, color="red")
            # splcieai_FN__splam_FN = ax.scatter(splam_noS_d_pred_prob[spliceai_FN__splam_FN_idices], splam_noS_a_pred_prob[spliceai_FN__splam_FN_idices], s = 5, color="green")
            # splcieai_TN__splam_FP = ax.scatter(splam_noS_d_pred_prob[spliceai_TN__splam_FP_idices], splam_noS_a_pred_prob[spliceai_TN__splam_FP_idices], s = 5, color="orange")
            # splcieai_FP__splam_TN = ax.scatter(splam_noS_d_pred_prob[spliceai_FP__splam_TN_idices], splam_noS_a_pred_prob[spliceai_FP__splam_TN_idices], s = 5, color="cyan")
            # splcieai_FP__splam_FP = ax.scatter(splam_noS_d_pred_prob[spliceai_FP__splam_FP_idices], splam_noS_a_pred_prob[spliceai_FP__splam_FP_idices], s = 5, color="purple")
            # # non_junc_legend = ax.scatter(score_d[non_junc_prob], score_a[non_junc_prob], s = 1)
            # ax.legend([splcieai_TP__splam_FN, splcieai_FN__splam_TP, splcieai_FN__splam_FN, splcieai_TN__splam_FP, splcieai_FP__splam_TN, splcieai_FP__splam_FP], ["splcieai_TP__splam_FN", "splcieai_FN__splam_TP", "splcieai_FN__splam_FN", "splcieai_TN__splam_FP", "splcieai_FP__splam_TN", "splcieai_FP__splam_FP"])
            # plt.xlim([0, 1])
            # plt.ylim([0, 1])
            # plt.show()

            # plot_scatter_plot(spliceai_d_label_prob[spliceai_FP_idices], spliceai_d_pred_prob, spliceai_a_label_prob, spliceai_a_pred_prob, figure_root+"spliceai.png")


            # print(spliceai_d_pred_prob[spliceai_d_pred_prob > threshold])

            # spliceai_d_label_prob[spliceai_d_label_prob > threshold]
            # spliceai_a_pred_prob[spliceai_a_pred_prob > threshold]
            # spliceai_a_label_prob[spliceai_a_label_prob > threshold]
            # spliceai_junc_name


            # for x, y in zip(spliceai_junc_name, splam_noS_junc_name):
            #     if (x[0] == y[0] and int(x[1]) == int(y[1]) and int(x[2]) == int(y[2]) and x[3] == y[3]):
            #         print(True)
            #     else:
            #         print(False)

            ################################################
            # spliceai scatter plot
            ################################################
            # plot_scatter_plot(spliceai_d_label_prob, spliceai_d_pred_prob, spliceai_a_label_prob, spliceai_a_pred_prob, figure_root+"spliceai.png")
            # # plot_min_linear_plot(spliceai_d_label_prob, spliceai_d_pred_prob, spliceai_a_label_prob, spliceai_a_pred_prob, figure_root+"spliceai_1d.png")
            # # plot_scatter_plot(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_label_prob, spliceai_N_a_pred_prob, figure_root+"spliceai_N.png")
            # # # plot_min_linear_plot(spliceai_N_d_label_prob, spliceai_N_d_pred_prob, spliceai_N_a_label_prob, spliceai_N_a_pred_prob, figure_root+"spliceai_N_1d.png")

            # # plot_scatter_plot(splam_S_d_label_prob, splam_S_d_pred_prob, splam_S_a_label_prob, splam_S_a_pred_prob, figure_root+"splam_shuffle.png")
            # # # plot_min_linear_plot(splam_S_d_label_prob, splam_S_d_pred_prob, splam_S_a_label_prob, splam_S_a_pred_prob, figure_root+"splam_shuffle_1d.png")

            # plot_scatter_plot(splam_noS_d_label_prob, splam_noS_d_pred_prob, splam_noS_a_label_prob, splam_noS_a_pred_prob, figure_root+"splam_noshuffle.png")
            # # plot_min_linear_plot(splam_noS_d_label_prob, splam_noS_d_pred_prob, splam_noS_a_label_prob, splam_noS_a_pred_prob, figure_root+"splam_noshuffle_1d.png")


if __name__ == "__main__":
    main()


