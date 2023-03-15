###############################################################################
# This file contains code to test the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
from keras.models import load_model
from utils import *
from constants import *

assert int(sys.argv[1]) in [80, 400, 2000, 10000]
CL = int(sys.argv[1])

###############################################################################
# Load model and test data
###############################################################################

BATCH_SIZE = 6
version = [1, 2, 3, 4, 5]

model = [[] for v in range(len(version))]

for v in range(len(version)):
    # Model provided by Illumina
#     model[v] = load_model('/home/kh.chao/Projects/PR3_SpliceAI_testing/models/SpliceNet' + str(CL)
#                           + '_g' + str(version[v]) + '.h5')
    # Model released on github
    model[v] = load_model('/Users/chaokuan-hao/Documents/Projects/PR3_SpliceAI_testing/models/SpliceNet'+sys.argv[1]+'_c' + str(version[v]) + '.h5')

# h5f = h5py.File(data_dir + 'dataset'
#                 + '_' + sys.argv[1] + '_' + sys.argv[2]
#                 + '.h5', 'w')

h5f = h5py.File(data_dir + 'dataset_test_0.h5', 'r')

num_idx = len(list(h5f.keys()))//2

###############################################################################
# Model testing
###############################################################################

start_time = time.time()

output_class_labels = ['Null', 'Acceptor', 'Donor']
# The three neurons per output correspond to no splicing, splice acceptor (AG)
# and splice donor (GT) respectively.

for output_class in [1, 2]:

    Y_true = [[] for t in range(1)]
    Y_pred = [[] for t in range(1)]

    for idx in range(num_idx):

        X = h5f['X' + str(idx)][:]
        Y = h5f['Y' + str(idx)][:]

        Xc, Yc = clip_datapoints(X, Y, CL, 1)
        # print("Xc: ", Xc)

        Yps = [np.zeros(Yc[0].shape) for t in range(1)]

        for v in range(len(version)):

            Yp = model[v].predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            for t in range(1):
                Yps[t] += Yp[t]/len(version)
        # Ensemble averaging (mean of the ensemble predictions is used)

        for t in range(1):

            is_expr = (Yc[t].sum(axis=(1,2)) >= 1)

            Y_true[t].extend(Yc[t][is_expr, :, output_class].flatten())
            Y_pred[t].extend(Yps[t][is_expr, :, output_class].flatten())

    print("\n\033[1m%s:\033[0m" % (output_class_labels[output_class]))

    for t in range(1):

        Y_true[t] = np.asarray(Y_true[t])
        Y_pred[t] = np.asarray(Y_pred[t])

        print(("Y_pred[t]: ", Y_pred[t]))

        print_topl_statistics(Y_true[t], Y_pred[t])


h5f.close()

print("--- %s seconds ---" % (time.time() - start_time))
print("--------------------------------------------------------------")

###############################################################################

