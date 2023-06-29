import os, sys
import h5py
import matplotlib.pyplot as plt


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from TEST_dataset import *
# from SpliceNN import *
# from SpliceNN_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
from tqdm import tqdm
import warnings

from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def plot_roc_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

def main(argv):
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
    input_sequence = 'CGATCTGACGTGGGTGTCATCGCATTATCGATATTGCAT'
    # Replace this with your custom sequence

    BATCH_SIZE = 100
    context = 10000
    SEQ_LEN = "600"
    paths = ('./models/spliceai{}.h5'.format(x) for x in range(1, 6))
    print("paths: ", paths)
    models = [load_model(resource_filename('spliceai', x)) for x in paths]

    data_dir='./'

    h5f = h5py.File("./INPUTS/dataset.h5", 'r')
    print(h5f.keys())

    X = h5f["X"]
    # [0:200]
    Y = h5f["Y"][0]
    # [0:200]

    # print("Y.shape: ", Y.shape)

    Y_pred = models[0].predict(X, batch_size=BATCH_SIZE)
    print(Y_pred.shape)
    # print(Y[:][150])
    # print(Y[:][450])
    # print(Y_pred[:][150])
    # print(Y_pred[:][450])


    is_expr = (Y.sum(axis=(1,2)) >= 1)
    # print("is_expr: ", is_expr)
    # print("is_expr: ", Y.sum(axis=(0,1)))
    # print("is_expr: ", Y.sum(axis=(0)))

    # Acceptor_YL = labels[is_expr, 1, :].flatten().numpy()
    Acceptor_YL = Y[:, :, 1].flatten()
    Acceptor_YP = Y_pred[:, :, 1].flatten()
    Donor_YL = Y[:, :, 2].flatten()
    Donor_YP = Y_pred[:, :, 2].flatten()

    # print("Acceptor_YL: ", Acceptor_YL.shape)
    # print("Acceptor_YP: ", Acceptor_YP.shape)

    # print("Donor_YL: ", Donor_YL.shape)
    # print("Donor_YP: ", Donor_YP.shape)

    A_YL = np.array(Y[:, :, 1])
    A_YP = np.array(Y_pred[:, :, 1])
    D_YL = np.array(Y[:, :, 2])
    D_YP = np.array(Y_pred[: , :, 2])

    # print("A_YL: ", A_YL.shape)
    # print("A_YP: ", A_YP.shape)

    # print("D_YL: ", D_YL.shape)
    # print("D_YP: ", D_YP.shape)

    plot_roc_curve(Acceptor_YL, Acceptor_YP, "Acceptor")
    plot_roc_curve(Donor_YL, Donor_YP, "Donor")
    plt.savefig("output_SRR1352129.png", dpi=300)

    for i in range(1, 10, 1):
        threshold = i * 0.1
        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, threshold, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        
        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, threshold, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
        D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, threshold, D_G_TP, D_G_FN, D_G_FP, D_G_TN)


        A_G_TP = int(A_G_TP)
        A_G_FN = int(A_G_FN)
        A_G_FP = int(A_G_FP)
        A_G_TN = int(A_G_TN)

        D_G_TP = int(D_G_TP)
        D_G_FN = int(D_G_FN)
        D_G_FP = int(D_G_FP)
        D_G_TN = int(D_G_TN)
        # print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Donor Acc: {epoch_donor_acc/len(train_loader):.3f} | Acceptor Acc: {epoch_acceptor_acc/len(train_loader):.3f}')
        print(">>>>> Threshold: ", threshold)
        print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
        print(f'Donor Accuracy   : {A_accuracy:.5f} | Donor AUC   : {A_auc:.5f} | Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
        print(f'Acceptor Accuracy: {D_accuracy:.5f} | Acceptor AUC: {D_auc:.5f} | Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
        # print ("Learning rate: %.5f" % (get_lr(optimizer)))
        print(">>>>>>>>>>>>>>>\n\n")



    # for key in h5f.keys():
    #     X = h5f[key]
    #     print(Y)



    # x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
    # y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

    # acceptor_prob = y[0, :, 1]
    # donor_prob = y[0, :, 2]

    # print("acceptor_prob: ", acceptor_prob)
    # print("donor_prob   : ", donor_prob)


if __name__ == "__main__":
    main(sys.argv[1:])
