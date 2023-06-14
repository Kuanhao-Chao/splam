"""
This file evaluate the testing dataset for the SPLAM model.

Parameters:
    MODEL_VERSION   : the version of the model.

    Input directory : "./INPUTS/"+SEQ_LEN+"bp/input_pos/"

    Output directory: "../src_tools_evaluation/splam_result/"+MODEL_VERSION+"/"
"""
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import torch.nn as nn
# from TEST_dataset import *
from splam_dataset_Chromsome import *
from SPLAM import *
# from splam_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from tqdm import tqdm
import warnings
from sklearn.metrics import precision_recall_curve, roc_curve
import pickle
import platform

warnings.filterwarnings("ignore")

MODEL_VERSION = "SPLAM_v11/"
JUNC_THRESHOLD = 0.1

def parse_junction(name):
    # print("name: ", name)
    res = name.split(":")
    strand = name[-2]
    chr_name = res[0]
    if strand == "+":
        start = int(res[1].split("-")[0])+200
        end = int(res[2].split("-")[1].split('(')[0])-200
    elif strand == "-":
        start = int(res[2].split("-")[0])+200
        end = int(res[1].split("-")[1].split('(')[0])-200
    # print("start: ", start)
    # print("end  : ", end)
    return (chr_name, start, end, strand)


def main():
    #############################
    # Global variable definition
    #############################
    EPOCH_NUM = 20
    BATCH_SIZE = 100
    N_WORKERS = 1

    #############################
    # Selecting device
    #############################
    device_str = None
    if torch.cuda.is_available():
        device_str = "cuda"
    else:
        if platform.system() == "Darwin":
            device_str = "mps"
        else:
            device_str = "cpu"
    device = torch.device(device_str)
    print(f"\033[1m[Info]: Use {device} now!\033[0m")

    MODEL = "./MODEL/"+MODEL_VERSION+"splam_14.pt"
    MODEL_OUTPUT_BASE = "../src_tools_evaluation/splam_result/"+MODEL_VERSION+"/"

    print(">> Using model: ", MODEL)

    os.makedirs(MODEL_OUTPUT_BASE, exist_ok=True)
    model = torch.load(MODEL)

    #############################
    # Model Initialization
    #############################
    print(f"[Info]: Finish loading model!",flush = True)
    print("model: ", model)

    #############################
    # Training Data initialization
    #############################

    criterion = nn.BCELoss()

    BATCH_SIZE = 100

    # TARGETS = ["pos", "pos_refseq_protein_isoforms", "pos_refseq_protein_alternative_only", "neg_1", "neg_5"]
    # TARGETS = ["pos", "pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
    TARGETS = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"]
    # TARGETS = ["pos_MANE", "pos_ALTS"]
    # TARGETS = ["neg_1_random"]
    # for shuffle in [True, False]:

    for shuffle in [False]:
        junc_counter = 0
        TYPE = "shuffle" if shuffle else "noshuffle"
        print("########################################")
        print(" Model: ", model)
        print("########################################")

        for target in TARGETS:

            os.makedirs(MODEL_OUTPUT_BASE+target, exist_ok=True)
            d_score_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.score.d."+TYPE+"."+target+".tsv"
            a_score_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.score.a."+TYPE+"."+target+".tsv"
            n_score_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.score.n."+TYPE+"."+target+".tsv"
            name_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.name."+TYPE+"."+target+".tsv"

            d_score_fw = open(d_score_tsv_f, "a")
            a_score_fw = open(a_score_tsv_f, "a")
            n_score_fw = open(n_score_tsv_f, "a")
            name_fw = open(name_tsv_f, "a")

            test_loader = get_eval_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS, shuffle, target)

            # test_iterator = iter(test_loader)
            print(f"[Info]: Finish loading data!", flush = True)
            print("valid_iterator: ", len(test_loader))
            LOG_OUTPUT_TEST_BASE = MODEL_OUTPUT_BASE + "/" + target + "/LOG/"
            os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)


            ############################
            # Log for testing
            ############################
            OUT_SCORE = LOG_OUTPUT_TEST_BASE + TYPE + "_junction_score.bed"
            test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
            test_log_acc = LOG_OUTPUT_TEST_BASE + "test_accuracy.txt"

            test_log_A_auprc = LOG_OUTPUT_TEST_BASE + "test_A_auprc.txt"
            test_log_A_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_A_threshold_precision.txt"
            test_log_A_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_A_threshold_recall.txt"
            test_log_D_auprc = LOG_OUTPUT_TEST_BASE + "test_D_auprc.txt"
            test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
            test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"
            test_log_J_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_J_threshold_precision.txt"
            test_log_J_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_J_threshold_recall.txt"

            fw_test_log_loss = open(test_log_loss, 'w')
            fw_test_log_acc = open(test_log_acc, 'w')

            fw_test_log_A_auprc = open(test_log_A_auprc, 'w')
            fw_test_log_A_threshold_precision = open(test_log_A_threshold_precision, 'w')
            fw_test_log_A_threshold_recall = open(test_log_A_threshold_recall, 'w')
            fw_test_log_D_auprc = open(test_log_D_auprc, 'w')
            fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
            fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')
            fw_test_log_J_threshold_precision = open(test_log_J_threshold_precision, 'w')
            fw_test_log_J_threshold_recall = open(test_log_J_threshold_recall, 'w')

            fw_junc_scores = open(OUT_SCORE, 'w')

            epoch_loss = 0
            epoch_acc = 0
            epoch_donor_acc = 0
            epoch_acceptor_acc = 0

            print("**********************")
            print("** Testing Dataset **")
            print("**********************")
            pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")

            A_G_TP = 1e-6
            A_G_FN = 1e-6
            A_G_FP = 1e-6
            A_G_TN = 1e-6
            D_G_TP = 1e-6
            D_G_FN = 1e-6
            D_G_FP = 1e-6
            D_G_TN = 1e-6

            J_G_TP = 1e-6
            J_G_FN = 1e-6
            J_G_FP = 1e-6
            J_G_TN = 1e-6
            #######################################
            # Important => setting model into evaluation mode
            #######################################

            All_Junction_YL_MIN = []
            All_Junction_YP_MIN = []

            All_Junction_YL_AVG = []
            All_Junction_YP_AVG = []

            SPLAM_Donor_YL = []
            SPLAM_Donor_YP = []
            SPLAM_Acceptor_YL = []
            SPLAM_Acceptor_YP = []
            SPLAM_junc_name = []

            model.eval()
            for batch_idx, data in enumerate(test_loader):
                # print("batch_idx: ", batch_idx)
                # DNAs:  torch.Size([40, 800, 4])
                # labels:  torch.Size([40, 1, 800, 3])
                DNAs, labels, chr = data

                # print("chr: ", chr)

                junc_name = map(parse_junction, chr)
                junc_name = list(junc_name)
                # print(chr.splt(":"))
                DNAs = DNAs.to(torch.float32).to(device)
                labels = labels.to(torch.float32).to(device)

                DNAs = torch.permute(DNAs, (0, 2, 1))
                labels = torch.permute(labels, (0, 2, 1))
                loss, yp = model_fn(DNAs, labels, model, criterion)


                #######################################
                # predicting all bp.
                #######################################
                is_expr = (labels.sum(axis=(1,2)) >= 1)
                # print("is_expr: ", is_expr)

                # Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
                Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
                Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
                Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
                Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

                A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
                A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
                D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
                D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()


                np.savetxt(a_score_fw, A_YP, delimiter=" ")
                np.savetxt(d_score_fw, D_YP, delimiter=" ")
                # np.savetxt(n_score_fw, n_scores.reshape((1,len(n_scores))), delimiter=" ")



                junction_labels_min, junction_scores_min = get_junc_scores(D_YL, A_YL, D_YP, A_YP, "min")
                junction_labels_avg, junction_scores_avg = get_junc_scores(D_YL, A_YL, D_YP, A_YP, "avg")
                donor_labels, donor_scores, acceptor_labels, acceptor_scores = get_donor_acceptor_scores(D_YL, A_YL, D_YP, A_YP)


                SPLAM_Donor_YL = np.concatenate((SPLAM_Donor_YL, donor_labels), axis=None)
                SPLAM_Donor_YP = np.concatenate((SPLAM_Donor_YP, donor_scores), axis=None)
                SPLAM_Acceptor_YL = np.concatenate((SPLAM_Acceptor_YL, acceptor_labels), axis=None)
                SPLAM_Acceptor_YP = np.concatenate((SPLAM_Acceptor_YP, acceptor_scores), axis=None)
                SPLAM_junc_name.extend(junc_name)


                All_Junction_YL_MIN.extend(junction_labels_min)
                All_Junction_YP_MIN.extend(junction_scores_min)
                All_Junction_YL_AVG.extend(junction_labels_avg)
                All_Junction_YP_AVG.extend(junction_scores_avg)

                J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)
                A_accuracy, A_auprc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
                D_accuracy, D_auprc = print_top_1_statistics(Donor_YL, Donor_YP)
                A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, JUNC_THRESHOLD, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
                D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, JUNC_THRESHOLD, D_G_TP, D_G_FN, D_G_FP, D_G_TN)

                batch_loss = loss.item()
                epoch_loss += loss.item()
                epoch_donor_acc += D_accuracy
                epoch_acceptor_acc += A_accuracy

                for idx in range(len(junc_name)):
                    chr_name, start, end, strand = junc_name[idx]
                    if strand == '+':
                        fw_junc_scores.write(chr_name+ '\t'+ str(start) + '\t' + str(end) + '\tJUNC_' + str(junc_counter) + '\t' + str(donor_labels[idx]) + '\t'+ strand + '\t' + str(donor_scores[idx]) + '\t' + str(acceptor_scores[idx]) + '\n')
                    elif strand == '-':
                        fw_junc_scores.write(chr_name+ '\t'+ str(start) + '\t' + str(end) + '\tJUNC_' + str(junc_counter) + '\t' + str(donor_labels[idx]) + '\t'+ strand+ '\t' + str(donor_scores[idx]) + '\t' + str(acceptor_scores[idx]) + '\n')
                    junc_counter += 1


                pbar.update(1)
                pbar.set_postfix(
                    batch_id=batch_idx,
                    idx_test=len(test_loader)*BATCH_SIZE,
                    loss=f"{batch_loss:.6f}",
                    # accuracy=f"{batch_acc:.6f}",
                    # A_accuracy=f"{A_accuracy:.6f}",
                    # D_accuracy=f"{D_accuracy:.6f}",
                    A_auprc = f"{A_auprc:.6f}",
                    D_auprc = f"{D_auprc:.6f}",
                    # A_TP=A_TP,
                    # A_FN=A_FN,
                    # A_FP=A_FP,
                    # A_TN=A_TN,
                    # D_TP=D_TP,
                    # D_FN=D_FN,
                    # D_FP=D_FP,
                    # D_TN=D_TN,
                    A_Precision=f"{A_TP/(A_TP+A_FP+1e-6):.6f}",
                    A_Recall=f"{A_TP/(A_TP+A_FN+1e-6):.6f}",
                    D_Precision=f"{D_TP/(D_TP+D_FP+1e-6):.6f}",
                    D_Recall=f"{D_TP/(D_TP+D_FN+1e-6):.6f}",
                    J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
                    J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
                )

                fw_test_log_loss.write(str(batch_loss)+ "\n")
                fw_test_log_A_auprc.write(str(A_auprc)+ "\n")
                fw_test_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1e-6):.6f}\n")
                fw_test_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1e-6):.6f}\n")
                fw_test_log_D_auprc.write(str(D_auprc)+ "\n")
                fw_test_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1e-6):.6f}\n")
                fw_test_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1e-6):.6f}\n")
                fw_test_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
                fw_test_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")


            pbar.close()
            fw_junc_scores.close()

            print(f'Epoch {batch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Donor top-k Acc: {epoch_donor_acc/len(test_loader):.3f} | Acceptor top-k Acc: {epoch_acceptor_acc/len(test_loader):.3f}')
            print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
            print(f'Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
            print(f'Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
            print("\n\n")


            print("All_Junction_YP_MIN: ", len(All_Junction_YP_MIN))
            print("All_Junction_YL_MIN: ", len(All_Junction_YL_MIN))
            print("All_Junction_YP_AVG: ", len(All_Junction_YP_AVG))
            print("All_Junction_YL_AVG: ", len(All_Junction_YL_AVG))

            print("Donor_YL   : ", len(SPLAM_Donor_YL))
            print("Donor_YP   : ", len(SPLAM_Donor_YP))
            print("Acceptor_YL: ", len(SPLAM_Acceptor_YL))
            print("Acceptor_YP: ", len(SPLAM_Acceptor_YP))

            with open(MODEL_OUTPUT_BASE + "/splam." + TYPE + "." + target + ".pkl", 'wb') as f:
                pickle.dump(All_Junction_YP_MIN, f)
                pickle.dump(All_Junction_YL_MIN, f)
                pickle.dump(All_Junction_YP_AVG, f)
                pickle.dump(All_Junction_YL_AVG, f)

            with open(MODEL_OUTPUT_BASE + "/splam.da." + TYPE + "." + target + ".pkl", 'wb') as f:
                pickle.dump(SPLAM_Donor_YL, f)
                pickle.dump(SPLAM_Donor_YP, f)
                pickle.dump(SPLAM_Acceptor_YL, f)
                pickle.dump(SPLAM_Acceptor_YP, f)
                pickle.dump(SPLAM_junc_name, f)



            d_score_fw.close()
            a_score_fw.close()
            n_score_fw.close()
            name_fw.close()
            ############################
            # Plotting ROC / PR curves
            ############################
            # plot_roc_curve(All_Junction_YL, All_Junction_YP, "Acceptor")
            # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_roc_"+TYPE+".png", dpi=300)
            # plt.close()
            # plot_pr_curve(All_Junction_YL, All_Junction_YP, "Acceptor")
            # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_pr_"+TYPE+".png", dpi=300)
            # plt.close()



            # ############################
            # # Log for testing
            # ############################
            # test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
            # test_log_D_threshold_accuracy = LOG_OUTPUT_TEST_BASE + "test_D_threshold_accuracy.txt"
            # test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
            # test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"

            # fw_test_log_loss = open(test_log_loss, 'w')
            # fw_test_log_D_threshold_accuracy = open(test_log_D_threshold_accuracy, 'w')
            # fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
            # fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')
            # # test_one_epoch(0, test_loader)
            # print("*********************")
            # print("** Testing Dataset **")
            # print("*********************")
            # epoch_loss = 0
            # epoch_acc = 0
            # pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")

            # J_G_TP = 1e-6
            # J_G_FN = 1e-6
            # J_G_FP = 1e-6
            # J_G_TN = 1e-6
            # All_Junction_YL = []
            # All_Junction_YP = []

            # model.eval()
            # # model.train()
            # with torch.no_grad():
            #     for batch_idx, data in enumerate(test_loader):
            #         # DNAs:  torch.Size([40, 800, 4])
            #         # labels:  torch.Size([40, 1, 800, 3])
            #         DNAs, labels, chr = data
            #         DNAs = DNAs.to(torch.float32).to(device)
            #         labels = labels.to(torch.float32).to(device)

            #         DNAs = torch.permute(DNAs, (0, 2, 1))
            #         loss, accuracy, yps = model_fn(DNAs, labels, model, criterion)


            #         #######################################
            #         # predicting splice / non-splice
            #         #######################################
            #         batch_loss = loss.item()
            #         batch_acc = accuracy
            #         epoch_loss += loss.item()
            #         epoch_acc += accuracy

            #         labels = labels.to("cpu").detach().numpy()
            #         yps = yps.to("cpu").detach().numpy()

            #         J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yps, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)

            #         All_Junction_YL.extend(labels)
            #         All_Junction_YP.extend(yps)

            #         pbar.update(1)
            #         pbar.set_postfix(
            #             epoch=batch_idx,
            #             idx_train=len(test_loader)*BATCH_SIZE,
            #             loss=f"{batch_loss:.6f}",
            #             accuracy=f"{batch_acc:.6f}",
            #             J_Precision=f"{J_TP/(J_TP+J_FP+1.e-10):.6f}",
            #             J_Recall=f"{J_TP/(J_TP+J_FN+1.e-10):.6f}"
            #         )
            #         fw_test_log_loss.write(str(batch_loss)+ "\n")
            #         fw_test_log_D_threshold_accuracy.write(str(batch_acc)+ "\n")
            #         fw_test_log_D_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1.e-10):.6f}\n")
            #         fw_test_log_D_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1.e-10):.6f}\n")
            # pbar.close()

            # TYPE = "shuffle" if shuffle else "noshuffle"

            # print("All_Junction_YP: ", All_Junction_YP)
            # print("All_Junction_YL: ", All_Junction_YL)

            # with open(MODEL_OUTPUT_BASE + "/splam." + TYPE + ".pkl", 'wb') as f:
            #     pickle.dump(All_Junction_YP, f)
            #     pickle.dump(All_Junction_YL, f)


        # ############################
        # # Plotting ROC / PR curves
        # ############################
        # plot_roc_curve(All_Junction_YL_MIN, All_Junction_YP_MIN, "Acceptor")
        # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_roc_"+TYPE+"_MIN.png", dpi=300)
        # plt.close()
        # plot_pr_curve(All_Junction_YL_MIN, All_Junction_YP_MIN, "Acceptor")
        # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_pr_"+TYPE+"_MIN.png", dpi=300)
        # plt.close()

        # plot_roc_curve(All_Junction_YL_AVG, All_Junction_YP_AVG, "Acceptor")
        # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_roc_"+TYPE+"_AVG.png", dpi=300)
        # plt.close()
        # plot_pr_curve(All_Junction_YL_AVG, All_Junction_YP_AVG, "Acceptor")
        # plt.savefig(MODEL_OUTPUT_BASE+"output_800bp_pr_"+TYPE+"_AVG.png", dpi=300)
        # plt.close()

        # print(f'Epoch {0+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Acc: {epoch_acc/len(test_loader):.3f}')
        # print(f'Expected #prediction: {len(test_loader)*BATCH_SIZE+0:03}')
        # print(f'Junction Precision  : {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
        # print("")
        # print("\n\n")


def plot_pr_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_roc_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')


if __name__ == "__main__":
    main()
