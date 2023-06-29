import tensorflow as tf
from tensorflow.keras import backend as K
import os, sys
import h5py
import matplotlib.pyplot as plt
from utils import *
import gc 

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt; plt.rcdefaults()
from tqdm import tqdm
import warnings

from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve
import pickle 

CONSTANT_SIZE = 10
def seq_name(seq):
    return seq[1:]

def main(argv):
    BATCH_SIZE = int(argv[0])
    BATCH_SIZE_BASE = 200

    TYPE = argv[1]
    # output_file = argv[2]
    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    splam_dir = "splam/"
    for output_file in output_files:
        path = './models/spliceai1.h5'
        print(">> path\t\t: ", path)
        print(">> output_file\t: ", output_file)
        model = load_model(resource_filename('spliceai', path))

        os.makedirs("./spliceai_result/", exist_ok=True)
        all_lines = []

        label = '.'

        if output_file == "pos" or output_file == "output_file" or output_file == "pos_ALTS":
            label = '+'
        if output_file == "neg_1" or output_file == "neg_random":
            label = '-'

        # print(">> label\t\t: ", label)
        da_faf = "./dataset/"+output_file+"/spliceai/spliceai."+TYPE+".juncs.seq.fa"
        print(">> da_faf\t\t: ", da_faf)
        print(">> pkl file\t\t: ", "./spliceai_result/spliceai."+TYPE+"."+output_file+".pkl")
        
        with open(da_faf, "r") as f:
            print("Processing : ", da_faf)
            lines = f.read().splitlines()
            all_lines = lines
        print("all_lines  : ", len(all_lines))

        #################################
        ## Processing 'POSITIVE' samples
        #################################
        COUNTER = 0
        pidx = 0
        
        pidx = BATCH_SIZE-BATCH_SIZE_BASE
        if pidx == 0:
            d_pred = []
            d_label = []
            a_pred = []
            a_label = []
            junction_name = []
        else: 
            with open("./spliceai_result/spliceai."+TYPE+"."+output_file+".pkl",'rb') as f:
                d_pred = pickle.load(f)
                d_label = pickle.load(f)
                a_pred = pickle.load(f)
                a_label = pickle.load(f)
                junction_name = pickle.load(f)

        seq = ""
        while pidx < len(all_lines):
            print(pidx)
            if pidx == BATCH_SIZE:
                with open("./spliceai_result/spliceai."+TYPE+"."+output_file+".pkl", 'wb') as f: 
                    pickle.dump(d_pred, f)
                    pickle.dump(d_label, f)
                    pickle.dump(a_pred, f)
                    pickle.dump(a_label, f)
                    pickle.dump(junction_name, f)
                exit()    
            if pidx % 2 == 0:
                chr, start, end, strand = all_lines[pidx].split(";")
                chr = chr[1:]
                junction_name.append((chr, start, end, strand))

                pass
            elif pidx % 2 == 1:
                seq = all_lines[pidx]

                X, Y = create_datapoints(seq, label)
                
                X = X[None, :]
                Y = np.array(Y)
                X = tf.convert_to_tensor(X, dtype=tf.float32)
                Y = tf.convert_to_tensor(Y, dtype=tf.float32)

                Y_pred = model.predict(X)
                
                K.clear_session()
                # del model
                gc.collect()
                COUNTER += 1
                print("X.shape     : ", X.shape)
                print("Y_pred.shape: ", Y_pred.shape)

                donor_p = Y_pred[0][200-1][2]
                acceptor_p = Y_pred[0][len(Y_pred)-200-1][1]
                print("(chr, start, end, strand): ", (chr, start, end, strand))
                print("donor_p    : ", donor_p)
                print("acceptor_p : ", acceptor_p)

                d_pred.append(donor_p)
                a_pred.append(acceptor_p)

                if output_file == "pos":
                    d_label.append(1)
                    a_label.append(1)
                else:
                    d_label.append(0)
                    a_label.append(0)
            pidx += 1
            print("====================")

            if pidx %100 == 0:
                print("pidx: ", pidx)
            # if pidx >= 10:
            #     break

        with open("./spliceai_result/spliceai."+TYPE+"."+output_file+".pkl", 'wb') as f: 
            pickle.dump(d_pred, f)
            pickle.dump(d_label, f)
            pickle.dump(a_pred, f)
            pickle.dump(a_label, f)
            pickle.dump(junction_name, f)

if __name__ == "__main__":
    main(sys.argv[1:])
