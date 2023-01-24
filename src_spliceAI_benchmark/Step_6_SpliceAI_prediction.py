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
    BATCH_SIZE_BASE = 400

    TYPE = argv[1]

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

    os.makedirs("./INPUT/", exist_ok=True)
    output_files = ["pos/", "neg_can/", "neg_noncan/", "neg_1/"]
    cum_lines = 0

    all_lines = []
    for output_file in output_files:
        # os.makedirs("./INPUT/"+output_file, exist_ok=True)
        da_faf = "./OUTPUT/"+output_file+"spliceai."+TYPE+".juncs.seq.fa"
        with open(da_faf, "r") as f:
            print("Processing ", da_faf)
            lines = f.read().splitlines()
            all_lines = all_lines + lines
    print("all_lines: ", len(all_lines))

    #################################
    ## Processing 'POSITIVE' samples
    #################################
    COUNTER = 0
    pidx = 0
    
    pidx = BATCH_SIZE-BATCH_SIZE_BASE
    if pidx == 0:
        d_pred_prob = []
        d_label_prob = []
        a_pred_prob = []
        a_label_prob = []
    else: 
        with open("./INPUT/spliceai."+TYPE+".pkl",'rb') as f:
            d_pred_prob = pickle.load(f)
            d_label_prob = pickle.load(f)
            a_pred_prob = pickle.load(f)
            a_label_prob = pickle.load(f)

    seq = ""
    all_length = len(all_lines)
    while pidx < all_length:
        print(pidx)
        if pidx == BATCH_SIZE:
            with open("./INPUT/spliceai."+TYPE+".pkl", 'wb') as f: 
                pickle.dump(d_pred_prob, f)
                pickle.dump(d_label_prob, f)
                pickle.dump(a_pred_prob, f)
                pickle.dump(a_label_prob, f)
            exit()    
        if pidx % 2 == 0:
            # seq_name = seq_name(line)
            pass
        elif pidx % 2 == 1:
            seq = all_lines[pidx]
            X, Y = create_datapoints(seq, '+')
            X = X[None, :]
            Y = np.array(Y)
            X = tf.convert_to_tensor(X, dtype=tf.float32)
            Y = tf.convert_to_tensor(Y, dtype=tf.float32)

            Y_pred = model(X)
            K.clear_session()
            # del model
            gc.collect()
            COUNTER += 1
            # print(Y_pred)
            print("X.shape     : ", X.shape)
            print("Y_pred.shape: ", Y_pred.shape)

            # print("Donor: ", Y_pred[0][200-2:200+2])
            # print("Acceptor: ", Y_pred[0][len(Y_pred)-200-2:len(Y_pred)-200+2])

            donor_p = Y_pred[0][200-1][2]
            acceptor_p = Y_pred[0][len(Y_pred)-200-1][1]
            print("donor_p   : ", donor_p)
            print("acceptor_p: ", acceptor_p)

            d_pred_prob.append(donor_p)
            a_pred_prob.append(acceptor_p)
            if output_file == "./OUTPUT/pos/":
                d_label_prob.append(1)
                a_label_prob.append(1)
            else:
                d_label_prob.append(0)
                a_label_prob.append(0)
        pidx += 1


        if pidx %100 == 0:
            print("pidx: ", pidx)
        # if pidx > 31:
        #     break


    with open("./INPUT/spliceai."+TYPE+".pkl", 'wb') as f: 
        pickle.dump(d_pred_prob, f)
        pickle.dump(d_label_prob, f)
        pickle.dump(a_pred_prob, f)
        pickle.dump(a_label_prob, f)

if __name__ == "__main__":
    main(sys.argv[1:])
