import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
# from Step_7_util import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve

def main():

    ratio = 12000
    os.makedirs("./IMG/spliceai/", exist_ok=True)
    os.makedirs("./IMG/splam/", exist_ok=True)
    
    d_pred_prob = []
    d_label_prob = []
    a_pred_prob = []
    a_label_prob = []

    splam_d_pred_prob = []
    splam_d_label_prob = []
    splam_a_pred_prob = []
    splam_a_label_prob = []

    for TYPE in ["N", "noN"]:
        with open("./INPUT/spliceai."+TYPE+".pkl",'rb') as f:
            d_pred_prob = pickle.load(f)
            d_label_prob = pickle.load(f)[:ratio//2]
            a_pred_prob = pickle.load(f)
            a_label_prob = pickle.load(f)[:ratio//2]

            d_pred_prob = [x.numpy() for x in d_pred_prob]
            a_pred_prob = [x.numpy() for x in a_pred_prob]

            # d_label_prob = [float(i) for i in d_label_prob]
            # a_label_prob = [float(i) for i in a_label_prob]

            # print("d_pred_prob : ", d_pred_prob)
            # print("d_label_prob: ", d_label_prob)
            # print("a_pred_prob : ", a_pred_prob)
            # print("a_label_prob: ", a_label_prob)

            print("d_pred_prob : ", len(d_pred_prob))
            print("d_label_prob: ", len(d_label_prob))
            print("a_pred_prob : ", len(a_pred_prob))
            print("a_label_prob: ", len(a_label_prob))
        

if __name__ == "__main__":
    main()