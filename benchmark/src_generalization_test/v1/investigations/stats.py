# If we choose a score threshold of 0.8, then the
# (1)sensitivity / (2)accuracy of Splam for (a)donor / (b)acceptor sites / (c)splice junctions  would be
# X%, Y%, and Z% for chimpanzee, mouse, and Arabidopsis 

import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
import pandas as pd
import itertools
from util import *
from sklearn.metrics import auc, accuracy_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, PrecisionRecallDisplay
from sklearn import svm
threshold = 0.15

###### INPUT DATA #########
def read_inputs(db):
    # positive
    noN_pos_df = pd.read_csv(f'../1_pos_test/data/aggregate/avg_data.noN.{db}.csv')
    noN_pos_df['true_label'] = 1
    print(len(noN_pos_df))

    # negative
    noN_neg_df = pd.read_csv(f'../2_neg_test/data/aggregate/avg_data.noN.{db}.csv')
    noN_neg_df['true_label'] = 0
    print(len(noN_neg_df))

    noN_pos_df = noN_pos_df.sample(n=10000, random_state=1091)
    noN_neg_df = noN_neg_df.sample(n=10000, random_state=5802)

    noN_merge_df = pd.concat([noN_pos_df, noN_neg_df], axis=0)

    print("noN_pos_df: ", len(noN_pos_df))
    print("noN_neg_df: ", len(noN_neg_df))
    print("noN_merge_df: ", len(noN_merge_df))
    print(noN_merge_df.head())

    return noN_merge_df

####### CALCULATE METRICS ###########
def calculate_metrics(df, site, model):
    true_labels = df['true_label']
    
    if site == 'donor':
        if model == 'splam':
            scores = df['d_score_splam']
        elif model == 'spliceai':
            scores = df['d_score_spliceai']
    elif site == 'acceptor':
        if model == 'splam':
            scores = df['a_score_splam']
        elif model == 'spliceai':
            scores = df['a_score_spliceai']
    elif site == 'both':
        if model == 'splam':
            t1 = df['d_score_splam']
            t2 = df['a_score_splam']
            scores = pd.concat([t1, t2], axis=0)
            true_labels = pd.concat([true_labels, true_labels], axis=0)
        elif model == 'spliceai':
            t1 = df['d_score_spliceai']
            t2 = df['a_score_spliceai']
            scores = pd.concat([t1, t2], axis=0)
            true_labels = pd.concat([true_labels, true_labels], axis=0)
            
    # Apply threshold to scores
    predicted_labels = scores >= threshold

    # Calculate sensitivity (true positive rate)
    true_positives = np.sum(np.logical_and(predicted_labels, true_labels))
    total_positives = np.sum(true_labels)
    sensitivity = true_positives / total_positives

    # Calculate accuracy
    total_samples = len(true_labels)
    correct_predictions = np.sum(predicted_labels == true_labels)
    accuracy = correct_predictions / total_samples

    return sensitivity, accuracy

###### RUNNER ######
def run():

    dbs = ["NHGRI_mPanTro3", "GRCm39", "TAIR10"]
    sites = ['donor', 'acceptor', 'both']
    models = ['splam', 'spliceai']

    with open(f'./result{threshold}.csv', 'w') as f:
        f.write('Database,Site,Model,Sensitivity,Accuracy\n')
        for db in dbs:
            df = read_inputs(db)
            for site, model in itertools.product(sites, models):
                print('-'*120)
                print(f'Calculating for database: {db}, site: {site}, model: {model}')
                sensitivity, accuracy = calculate_metrics(df, site, model)
                print(f'\tSensitivity: {sensitivity}\n\tAccuracy: {accuracy}\n')
                
                f.write(f'{db},{site},{model},{sensitivity},{accuracy}\n')


if __name__ == '__main__':
    run()