import tensorflow as tf
from tensorflow.keras import backend as K
import os, sys
import h5py
import matplotlib.pyplot as plt
from utils import *
import gc 
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# from TEST_dataset import *
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
import pickle 

# config = tf.compat.v1.ConfigProto(log_device_placement=True)
# config.gpu_options.allow_growth = True
# session = tf.Session(config=config)
# K.set_session(session)


CONSTANT_SIZE = 10
def seq_name(seq):
    return seq[1:]

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
    BATCH_SIZE = int(argv[0])
    COUNTER = 0
    # Replace this with your custom sequence

    # BATCH_SIZE = 100
    context = 10000
    SEQ_LEN = "600"
    # paths = ('./models/spliceai{}.h5'.format(x) for x in range(1, 2))
    # print("paths: ", paths)
    # models = [load_model(resource_filename('spliceai', x)) for x in paths]
    path = './models/spliceai1.h5'
    print("path: ", path)
    model = load_model(resource_filename('spliceai', path))

    data_dir='./'





    os.makedirs("./INPUT/", exist_ok=True)
    # h5f2 = h5py.File("./INPUT/dataset.h5", 'w')

    da_faf = "./OUTPUT/spliceai.juncs.seq.fa"

    Y_ALL = []
    # Y_batch = [[] for t in range(1)]
        
    #################################
    ## Processing 'POSITIVE' samples
    #################################
    pidx = 0
    # max_len = 0

    # 0.1 - 0.9 (9)
    # D_TP = [0]*9
    # D_TN = [0]*9
    # D_FP = [0]*9
    # D_FN = [0]*9

    # A_TP = [0]*9
    # A_TN = [0]*9
    # A_FP = [0]*9
    # A_FN = [0]*9

    # J_TP = [0]*9
    # J_TN = [0]*9
    # J_FP = [0]*9
    # J_FN = [0]*9

    COUNTER = 0

    POS_NUM = 3000
    
    pidx = BATCH_SIZE-400
    if pidx == 0:
        d_pred_prob = []
        d_label_prob = []
        a_pred_prob = []
        a_label_prob = []
    else: 
        with open("./INPUT/spliceai_"+str(pidx)+".pkl",'rb') as f:
            d_pred_prob = pickle.load(f)
            d_label_prob = pickle.load(f)
            a_pred_prob = pickle.load(f)
            a_label_prob = pickle.load(f)

    with open(da_faf, "r") as f:
        print("Processing ", da_faf)
        lines = f.read().splitlines()
        seq_name = ""
        seq = ""
        # for pidx in range(len(lines)):
        while pidx < len(lines):
            print(pidx)
            if pidx == BATCH_SIZE:
                with open("./INPUT/spliceai_"+str(pidx)+".pkl", 'wb') as f: 
                    # d_pred_prob = []
                    # d_label_prob = [1]*POS_NUM + [0]*POS_NUM
                    # a_pred_prob = []
                    # a_label_prob = [1]*POS_NUM + [0]*POS_NUM

                    pickle.dump(d_pred_prob, f)
                    pickle.dump(d_label_prob, f)
                    pickle.dump(a_pred_prob, f)
                    pickle.dump(a_label_prob, f)
                exit()    
            # print(line)
            if pidx % 2 == 0:
                # seq_name = seq_name(line)
                pass
            elif pidx % 2 == 1:
                seq = lines[pidx]
                # print(seq)
                X, Y = create_datapoints(seq, '+')
                # if X.shape[0] > max_len:
                #     max_len = X.shape[0] 
                    # print(COUNTER, "\tX.shape: ", X.shape, "\tY.shape: ", Y[0].shape)
                X = X[None, :]
                Y = np.array(Y)

                X = tf.convert_to_tensor(X, dtype=tf.float32)
                Y = tf.convert_to_tensor(Y, dtype=tf.float32)
                

                Y_pred = model(X)
                K.clear_session()
                # del model
                gc.collect()

                # Y_ALL.append(Y_pred)
                COUNTER += 1
                # print(Y_pred)
                # print(Y_pred.shape)

                donor_p = Y_pred[0][100-1][2]
                acceptor_p = Y_pred[0][len(Y_pred[0])-1-100+1][1]
                
                d_pred_prob.append(donor_p)
                a_pred_prob.append(acceptor_p)
                if pidx//2 < POS_NUM:
                    d_label_prob.append(1)
                    a_label_prob.append(1)
                else:
                    d_label_prob.append(0)
                    a_label_prob.append(0)

                # for idx in range(0,9,1):
                #     threshold = (idx+1)*0.1
                    # print("threshold:  ", threshold)
                    # print("Donor   : ", np.where(Y_pred[:, :, 2] > threshold))
                    # print("Acceptor: ", np.where(Y_pred[:, :, 1] > threshold))

                    # if pidx//2 < POS_NUM:
                    #     # print("INNNN!!POS_NUM")
                    #     if donor_p > threshold:
                    #         D_TP[idx]+=1
                    #     else:
                    #         D_FN[idx]+=1

                    #     if acceptor_p > threshold:
                    #         A_TP[idx]+=1
                    #     else:
                    #         A_FN[idx]+=1

                    #     if donor_p > threshold and acceptor_p > threshold:
                    #         J_TP[idx]+=1
                    #     else:
                    #         J_FN[idx]+=1
                    # else:
                    #     # print("OUTOUT!!POS_NUM")
                    #     if donor_p > threshold:
                    #         D_FP[idx]+=1
                    #     else:
                    #         D_TN[idx]+=1

                    #     if acceptor_p > threshold:
                    #         A_FP[idx]+=1
                    #     else:
                    #         A_TN[idx]+=1

                    #     if donor_p > threshold and acceptor_p > threshold:
                    #         J_FP[idx]+=1
                    #     else:
                    #         J_TN[idx]+=1
            


            # if COUNTER > 10:
            #     break
            pidx += 1


            if pidx %100 == 0:
                print("pidx: ", pidx)
                # print("Delet model & reload")
                # K.clear_session()
                # del model
                # gc.collect()
                # model = load_model(resource_filename('spliceai', path))
                # with open("./INPUT/spliceai_"+str(pidx)+".pkl", 'wb') as f: 
                #     # d_pred_prob = []
                #     # d_label_prob = [1]*POS_NUM + [0]*POS_NUM
                #     # a_pred_prob = []
                #     # a_label_prob = [1]*POS_NUM + [0]*POS_NUM
                #     d_label_prob[1500]=0
                #     a_label_prob[1500]=0
                #     pickle.dump(d_pred_prob, f)
                #     pickle.dump(d_label_prob, f)
                #     pickle.dump(a_pred_prob, f)
                #     pickle.dump(a_label_prob, f)


            # print("pidx: ", pidx)
            # if pidx > 20:
            #     break
    # print("D_TP: ", D_TP)
    # print("D_FN: ", D_FN)
    # print("D_FP: ", D_FP)
    # print("D_TN: ", D_TN)
    # print()
    # print("A_TP: ", A_TP)
    # print("A_FN: ", A_FN)
    # print("A_FP: ", A_FP)
    # print("A_TN: ", A_TN)
    # print()
    # print("J_TP: ", J_TP)
    # print("J_FN: ", J_FN)
    # print("J_FP: ", J_FP)
    # print("J_TN: ", J_TN)
    # print()

    with open("./INPUT/spliceai_"+str(pidx)+".pkl", 'wb') as f: 
        # d_pred_prob = []
        # d_label_prob = [1]*POS_NUM + [0]*POS_NUM
        # a_pred_prob = []
        # a_label_prob = [1]*POS_NUM + [0]*POS_NUM

        pickle.dump(d_pred_prob, f)
        pickle.dump(d_label_prob, f)
        pickle.dump(a_pred_prob, f)
        pickle.dump(a_label_prob, f)



    # # h5f = h5py.File("./INPUTS/dataset.h5", 'r')
    # # print(h5f.keys())

    # # X = h5f["X"]
    # # # [0:200]
    # # Y = h5f["Y"][0]
    # # # [0:200]

    # # print("Y.shape: ", Y.shape)

    # Y_pred = models[0].predict(X, batch_size=BATCH_SIZE)
    # print(Y_pred.shape)
    # # print(Y[:][150])
    # # print(Y[:][450])
    # # print(Y_pred[:][150])
    # # print(Y_pred[:][450])


    # is_expr = (Y.sum(axis=(1,2)) >= 1)
    # # print("is_expr: ", is_expr)
    # # print("is_expr: ", Y.sum(axis=(0,1)))
    # # print("is_expr: ", Y.sum(axis=(0)))

    # # Acceptor_YL = labels[is_expr, 1, :].flatten().numpy()
    # Acceptor_YL = Y[:, :, 1].flatten()
    # Acceptor_YP = Y_pred[:, :, 1].flatten()
    # Donor_YL = Y[:, :, 2].flatten()
    # Donor_YP = Y_pred[:, :, 2].flatten()

    # # print("Acceptor_YL: ", Acceptor_YL.shape)
    # # print("Acceptor_YP: ", Acceptor_YP.shape)

    # # print("Donor_YL: ", Donor_YL.shape)
    # # print("Donor_YP: ", Donor_YP.shape)

    # A_YL = np.array(Y[:, :, 1])
    # A_YP = np.array(Y_pred[:, :, 1])
    # D_YL = np.array(Y[:, :, 2])
    # D_YP = np.array(Y_pred[: , :, 2])

    # # print("A_YL: ", A_YL.shape)
    # # print("A_YP: ", A_YP.shape)

    # # print("D_YL: ", D_YL.shape)
    # # print("D_YP: ", D_YP.shape)

    # plot_roc_curve(Acceptor_YL, Acceptor_YP, "Acceptor")
    # plot_roc_curve(Donor_YL, Donor_YP, "Donor")
    # plt.savefig("output_SRR1352129.png", dpi=300)

    # for i in range(1, 10, 1):
    #     threshold = i * 0.1
    #     J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, threshold, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        
    #     A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
    #     D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
    #     A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, threshold, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
    #     D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, threshold, D_G_TP, D_G_FN, D_G_FP, D_G_TN)


    #     A_G_TP = int(A_G_TP)
    #     A_G_FN = int(A_G_FN)
    #     A_G_FP = int(A_G_FP)
    #     A_G_TN = int(A_G_TN)

    #     D_G_TP = int(D_G_TP)
    #     D_G_FN = int(D_G_FN)
    #     D_G_FP = int(D_G_FP)
    #     D_G_TN = int(D_G_TN)
    #     # print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Donor Acc: {epoch_donor_acc/len(train_loader):.3f} | Acceptor Acc: {epoch_acceptor_acc/len(train_loader):.3f}')
    #     print(">>>>> Threshold: ", threshold)
    #     print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    #     print(f'Donor Accuracy   : {A_accuracy:.5f} | Donor AUC   : {A_auc:.5f} | Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
    #     print(f'Acceptor Accuracy: {D_accuracy:.5f} | Acceptor AUC: {D_auc:.5f} | Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
    #     # print ("Learning rate: %.5f" % (get_lr(optimizer)))
    #     print(">>>>>>>>>>>>>>>\n\n")



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
