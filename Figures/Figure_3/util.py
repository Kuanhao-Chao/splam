import matplotlib.pyplot as plt
import numpy as np 
import warnings 

def calc_precision_recall(threshold, y_true, y_pred):
    # Instantiate counters
    TP = 0
    FP = 0
    FN = 0
    TN = 0

    # Determine whether each prediction is TP, FP, TN, or FN
    for i in range(len(y_true)): 
        if y_true[i]==1 and y_pred[i]>=threshold:
            TP += 1
        if y_true[i]==1 and y_pred[i]<threshold:
            FN += 1
        if y_true[i]==0 and y_pred[i]>=threshold:
            FP += 1
        if y_true[i]==0 and y_pred[i]<threshold:
            TN += 1

    try:
        precision = TP / (TP + FP)
    except:
        precision = 1
    
    try:
        recall = TP / (TP + FN)
    except:
        recall = 1

    return precision, recall


def stable_cumsum(arr, axis=None, rtol=1e-05, atol=1e-08):
    """Use high precision for cumsum and check that final value matches sum.
    Warns if the final cumulative sum does not match the sum (up to the chosen
    tolerance).
    Parameters
    ----------
    arr : array-like
        To be cumulatively summed as flat.
    axis : int, default=None
        Axis along which the cumulative sum is computed.
        The default (None) is to compute the cumsum over the flattened array.
    rtol : float, default=1e-05
        Relative tolerance, see ``np.allclose``.
    atol : float, default=1e-08
        Absolute tolerance, see ``np.allclose``.
    Returns
    -------
    out : ndarray
        Array with the cumulative sums along the chosen axis.
    """
    out = np.cumsum(arr, axis=axis, dtype=np.float64)
    expected = np.sum(arr, axis=axis, dtype=np.float64)
    if not np.all(
        np.isclose(
            out.take(-1, axis=axis), expected, rtol=rtol, atol=atol, equal_nan=True
        )
    ):
        warnings.warn(
            "cumsum was found to be unstable: "
            "its last element does not correspond to sum",
            RuntimeWarning,
        )
    return out

def _binary_clf_curve(y_true, y_prob, pos_label=None, sample_weight=None):
    desc_score_indices = np.argsort(y_prob, kind="mergesort")[::-1]
    y_prob = y_prob[desc_score_indices]
    y_true = y_true[desc_score_indices]

    distinct_value_indices = np.where(np.diff(y_prob))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_prob.size - 1]

    # for threshold_idx in threshold_idxs:
    #     print("threshold_idx: ", threshold_idx)

    # accumulate the true positives with decreasing threshold
    tps = stable_cumsum(y_true)[threshold_idxs]
    fps = stable_cumsum((1 - y_true))[threshold_idxs]
    return fps, tps, y_prob[threshold_idxs]


def precision_recall_curve_self(y_true, y_prob, label, option):
    y_true = np.array(y_true)
    y_prob = np.array(y_prob)
    if option == "sklearn":

        fps, tps, thresholds = _binary_clf_curve(
            y_true, y_prob, pos_label=label
        )
        ps = tps + fps
        # Initialize the result array with zeros to make sure that precision[ps == 0]
        # does not contain uninitialized values.
        precision = np.zeros_like(tps)
        np.divide(tps, ps, out=precision, where=(ps != 0))

        # When no positive label in y_true, recall is set to 1 for all thresholds
        # tps[-1] == 0 <=> y_true == all negative labels
        if tps[-1] == 0:
            warnings.warn(
                "No positive class found in y_true, "
                "recall is set to one for all thresholds."
            )
            recall = np.ones_like(tps)
        else:
            recall = tps / tps[-1]

        # reverse the outputs so recall is decreasing
        sl = slice(None, None, -1)
        return np.hstack((precision[sl], 1)), np.hstack((recall[sl], 0)), thresholds[sl]
    elif option == "self":
        # 0.1 - 0.9 (9)
        TPs = [0]*1001
        TNs = [0]*1001
        FPs = [0]*1001
        FNs = [0]*1001
        l2_precision_scores = []
        l2_recall_scores = []
        thresholds = []
        for idx in range(0, 1001, 1):
            threshold = (idx)/1000
            thresholds.append(threshold)
            print("threshold: ", threshold)

            labels_1 = np.where(y_true == 1)
            labels_0 = np.where(y_true == 0)
            thre = np.where(y_prob >= threshold)
            thre_0 = np.where(y_prob < threshold)

            TP = len(np.intersect1d(labels_1, thre))
            FN = len(np.setdiff1d(labels_1, thre))
            FP = len(np.setdiff1d(thre, labels_1))
            TN = len(np.intersect1d(labels_0, thre_0))
            # TNs = len(true_y) - TPs - FNs - FPs
            print("\tDonor TP: ", TP)
            print("\tDonor FN: ", FN)
            print("\tDonor FP: ", FP)
            print("\tDonor TN: ", TN)
            TPs[idx] = TP
            TNs[idx] = TN
            FPs[idx] = FP
            FNs[idx] = FN
            try:
                precision = TP / (TP + FP)
            except:
                precision = 1
            
            try:
                recall = TP / (TP + FN)
            except:
                recall = 1
            l2_precision_scores.append(precision)
            l2_recall_scores.append(recall)

        # # Containers for true positive / false positive rates
        # l2_precision_scores = []
        # l2_recall_scores = []

        # # Define probability thresholds to use, between 0 and 1
        # probability_thresholds = np.linspace(0,1,num=1000)

        # # Find true positive / false positive rate for each threshold
        # for p in probability_thresholds:
        #     precision, recall = calc_precision_recall(p, y_true, y_prob_J)
                
        #     l2_precision_scores.append(precision)
        #     l2_recall_scores.append(recall)

        return l2_precision_scores, l2_recall_scores, thresholds        



def _binary_clf_curve_J_level(y_true, y_prob, pos_label=None, sample_weight=None):
    desc_score_indices = np.argsort(y_prob, kind="mergesort")[::-1]

    y_prob = y_prob[desc_score_indices]
    y_true = y_true[desc_score_indices]

    distinct_value_indices = np.where(np.diff(y_prob))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_prob.size - 1]

    # for threshold_idx in threshold_idxs:
    #     print("threshold_idx: ", threshold_idx)

    # accumulate the true positives with decreasing threshold
    tps = stable_cumsum(y_true)[threshold_idxs]
    fps = stable_cumsum((1 - y_true))[threshold_idxs]
    return fps, tps, y_prob[threshold_idxs]


def precision_recall_curve_J_level_self(y_true, y_prob_d, y_prob_a, label, option, choice):
    y_true = np.array(y_true)
    y_prob_d = np.array(y_prob_d)
    y_prob_a = np.array(y_prob_a)

    if choice == "min":
        y_prob_J = np.minimum(y_prob_d, y_prob_a)
    elif choice == "avg":
        y_prob_J = np.mean([y_prob_d, y_prob_a], axis = 0)

    # y_prob_J = np.maximum(y_prob_d, y_prob_a)
    # print("y_prob_J: ", y_prob_J)

    if option == "sklearn":
        fps, tps, thresholds = _binary_clf_curve_J_level(
            y_true, y_prob_J, pos_label=label
        )
        ps = tps + fps
        # Initialize the result array with zeros to make sure that precision[ps == 0]
        # does not contain uninitialized values.
        precision = np.zeros_like(tps)
        np.divide(tps, ps, out=precision, where=(ps != 0))

        # When no positive label in y_true, recall is set to 1 for all thresholds
        # tps[-1] == 0 <=> y_true == all negative labels
        if tps[-1] == 0:
            warnings.warn(
                "No positive class found in y_true, "
                "recall is set to one for all thresholds."
            )
            recall = np.ones_like(tps)
        else:
            recall = tps / tps[-1]

        # reverse the outputs so recall is decreasing
        sl = slice(None, None, -1)
        return np.hstack((precision[sl], 1)), np.hstack((recall[sl], 0)), thresholds[sl]
    elif option == "self":
        # 0.1 - 0.9 (9)
        TPs = [0]*1001
        TNs = [0]*1001
        FPs = [0]*1001
        FNs = [0]*1001
        l2_precision_scores = []
        l2_recall_scores = []
        thresholds = []
        for idx in range(0, 1001, 1):
            threshold = (idx)/1000
            thresholds.append(threshold)
            print("threshold: ", threshold)

            labels_1 = np.where(y_true == 1)
            labels_0 = np.where(y_true == 0)
            thre = np.where(y_prob_J >= threshold)
            thre_0 = np.where(y_prob_J < threshold)

            TP = len(np.intersect1d(labels_1, thre))
            FN = len(np.setdiff1d(labels_1, thre))
            FP = len(np.setdiff1d(thre, labels_1))
            TN = len(np.intersect1d(labels_0, thre_0))
            # TNs = len(true_y) - TPs - FNs - FPs
            print("\tDonor TP: ", TP)
            print("\tDonor FN: ", FN)
            print("\tDonor FP: ", FP)
            print("\tDonor TN: ", TN)
            TPs[idx] = TP
            TNs[idx] = TN
            FPs[idx] = FP
            FNs[idx] = FN
            try:
                precision = TP / (TP + FP)
            except:
                precision = 1
            
            try:
                recall = TP / (TP + FN)
            except:
                recall = 1
            l2_precision_scores.append(precision)
            l2_recall_scores.append(recall)

        # # Containers for true positive / false positive rates
        # l2_precision_scores = []
        # l2_recall_scores = []

        # # Define probability thresholds to use, between 0 and 1
        # probability_thresholds = np.linspace(0,1,num=1000)

        # # Find true positive / false positive rate for each threshold
        # for p in probability_thresholds:
        #     precision, recall = calc_precision_recall(p, y_true, y_prob_J)
                
        #     l2_precision_scores.append(precision)
        #     l2_recall_scores.append(recall)

        return l2_precision_scores, l2_recall_scores, thresholds        
    
def calc_FPR_TPR(threshold, y_true, y_pred):
    # Instantiate counters
    TP = 0
    FP = 0
    FN = 0
    TN = 0

    # Determine whether each prediction is TP, FP, TN, or FN
    for i in range(len(y_true)): 
        if y_true[i]==1 and y_pred[i]>=threshold:
            TP += 1
        if y_true[i]==1 and y_pred[i]<threshold:
            FN += 1
        if y_true[i]==0 and y_pred[i]>=threshold:
            FP += 1
        if y_true[i]==0 and y_pred[i]<threshold:
            TN += 1
    try:
        FPR = FP / (TN + FP)
    except:
        FPR = 1
    
    try:
        TPR = TP / (TP + FN)
    except:
        TPR = 1
    return FPR, TPR


def ROC_curve_self(y_true, y_prob, label, option):
    y_true = np.array(y_true)
    y_prob = np.array(y_prob)
    if option == "sklearn":
        pass
    elif option == "self":

        # 0.1 - 0.9 (9)
        TPs = [0]*1001
        TNs = [0]*1001
        FPs = [0]*1001
        FNs = [0]*1001
        l2_FPR_scores = []
        l2_TPR_scores = []
        thresholds = []
        for idx in range(0, 1001, 1):
            threshold = (idx)/1000
            thresholds.append(threshold)
            print("threshold: ", threshold)

            labels_1 = np.where(y_true == 1)
            labels_0 = np.where(y_true == 0)
            thre = np.where(y_prob >= threshold)
            thre_0 = np.where(y_prob < threshold)

            TP = len(np.intersect1d(labels_1, thre))
            FN = len(np.setdiff1d(labels_1, thre))
            FP = len(np.setdiff1d(thre, labels_1))
            TN = len(np.intersect1d(labels_0, thre_0))
            # TNs = len(true_y) - TPs - FNs - FPs
            print("\tDonor TP: ", TP)
            print("\tDonor FN: ", FN)
            print("\tDonor FP: ", FP)
            print("\tDonor TN: ", TN)
            TPs[idx] = TP
            TNs[idx] = TN
            FPs[idx] = FP
            FNs[idx] = FN
            try:
                FPR = FP / (TN + FP)
            except:
                FPR = 1
            
            try:
                TPR = TP / (TP + FN)
            except:
                TPR = 1
            l2_FPR_scores.append(FPR)
            l2_TPR_scores.append(TPR)



        # # Containers for true positive / false positive rates
        # l2_FPR_scores = []
        # l2_TPR_scores = []

        # # Define probability thresholds to use, between 0 and 1
        # probability_thresholds = np.linspace(0,1,num=1000)

        # # Find true positive / false positive rate for each threshold
        # for p in probability_thresholds:
        #     FPR, TPR = calc_FPR_TPR(p, y_true, y_prob)
                
        #     l2_FPR_scores.append(FPR)
        #     l2_TPR_scores.append(TPR)
        return l2_FPR_scores, l2_TPR_scores, thresholds


def roc_curve_J_level_self(y_true, y_prob_d, y_prob_a, label, option, choice):
    y_true = np.array(y_true)
    y_prob_d = np.array(y_prob_d)
    y_prob_a = np.array(y_prob_a)

    if choice == "min":
        y_prob_J = np.minimum(y_prob_d, y_prob_a)
    elif choice == "avg":
        y_prob_J = np.mean([y_prob_d, y_prob_a], axis = 0)


    if option == "sklearn":
        fps, tps, thresholds = _binary_clf_curve_J_level(
            y_true, y_prob_J, pos_label=label
        )

        # if drop_intermediate and len(fps) > 2:
        #     optimal_idxs = np.where(
        #         np.r_[True, np.logical_or(np.diff(fps, 2), np.diff(tps, 2)), True]
        #     )[0]
        #     fps = fps[optimal_idxs]
        #     tps = tps[optimal_idxs]
        #     thresholds = thresholds[optimal_idxs]

        # Add an extra threshold position
        # to make sure that the curve starts at (0, 0)
        tps = np.r_[0, tps]
        fps = np.r_[0, fps]
        thresholds = np.r_[thresholds[0] + 1, thresholds]

        if fps[-1] <= 0:
            warnings.warn(
                "No negative samples in y_true, false positive value should be meaningless",
                # UndefinedMetricWarning,
            )
            fpr = np.repeat(np.nan, fps.shape)
        else:
            fpr = fps / fps[-1]

        if tps[-1] <= 0:
            warnings.warn(
                "No positive samples in y_true, true positive value should be meaningless",
                # UndefinedMetricWarning,
            )
            tpr = np.repeat(np.nan, tps.shape)
        else:
            tpr = tps / tps[-1]

        return fpr, tpr, thresholds
    elif option == "self":
        # 0.1 - 0.9 (9)
        TPs = [0]*1001
        TNs = [0]*1001
        FPs = [0]*1001
        FNs = [0]*1001
        l2_FPR_scores = []
        l2_TPR_scores = []
        thresholds = []
        for idx in range(0, 1001, 1):
            threshold = (idx)/1000
            thresholds.append(threshold)
            print("threshold: ", threshold)

            labels_1 = np.where(y_true == 1)
            labels_0 = np.where(y_true == 0)
            thre = np.where(y_prob_J >= threshold)
            thre_0 = np.where(y_prob_J < threshold)

            TP = len(np.intersect1d(labels_1, thre))
            FN = len(np.setdiff1d(labels_1, thre))
            FP = len(np.setdiff1d(thre, labels_1))
            TN = len(np.intersect1d(labels_0, thre_0))
            # TNs = len(true_y) - TPs - FNs - FPs
            print("\tDonor TP: ", TP)
            print("\tDonor FN: ", FN)
            print("\tDonor FP: ", FP)
            print("\tDonor TN: ", TN)
            TPs[idx] = TP
            TNs[idx] = TN
            FPs[idx] = FP
            FNs[idx] = FN
            TPs[idx] = TP
            TNs[idx] = TN
            FPs[idx] = FP
            FNs[idx] = FN
            try:
                FPR = FP / (TN + FP)
            except:
                FPR = 1
            
            try:
                TPR = TP / (TP + FN)
            except:
                TPR = 1
            l2_FPR_scores.append(FPR)
            l2_TPR_scores.append(TPR)

        # # Containers for true positive / false positive rates
        # l2_FPR_scores = []
        # l2_TPR_scores = []

        # # Define probability thresholds to use, between 0 and 1
        # probability_thresholds = np.linspace(0,1,num=1000)

        # # Find true positive / false positive rate for each threshold
        # for p in probability_thresholds:
        #     FPR, TPR = calc_FPR_TPR(p, y_true, y_prob_J)
                
        #     l2_FPR_scores.append(FPR)
        #     l2_TPR_scores.append(TPR)
        return l2_FPR_scores, l2_TPR_scores, thresholds


