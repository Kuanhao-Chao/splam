###############################################################################
# This file contains the code to train the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
# import tensorflow.keras.backend as kb
import tensorflow as tf
# from spliceai import *
# from SpliceNN import *
from SpliceNN_dataset import *
from spliceai import *
from utils_SpliceNN import *
from multi_gpu import *
from constants import *
from SpliceNN_keras import *

assert int(sys.argv[1]) in [80, 400, 2000, 10000]

###############################################################################
# Model
###############################################################################

L = 32
N_GPUS = 2

if int(sys.argv[1]) == 80:
    W = np.asarray([11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1])
    BATCH_SIZE = 18*N_GPUS
elif int(sys.argv[1]) == 400:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
    BATCH_SIZE = 18*N_GPUS
elif int(sys.argv[1]) == 2000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10])
    BATCH_SIZE = 12*N_GPUS
elif int(sys.argv[1]) == 10000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21, 41, 41, 41, 41])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10, 25, 25, 25, 25])

    BATCH_SIZE = 6*N_GPUS

# Hyper-parameters:
# L: Number of convolution kernels
# W: Convolution window size in each residual unit
# AR: Atrous rate in each residual unit

CL = 2 * np.sum(AR*(W-1))
assert CL <= CL_max and CL == int(sys.argv[1])
print("\033[1mContext nucleotides: %d\033[0m" % (CL))
print("\033[1mSequence length (output): %d\033[0m" % (SL))

# model = SpliceAI(L, W, AR)
# model.summary()
# model_m = make_parallel(model, N_GPUS)
# model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')






embed_dim = 5000  # Embedding size for each token
num_heads = 2  # Number of attention heads
ff_dim = 32  # Hidden layer size in feed forward network inside transformer

inputs = layers.Input(shape=(embed_dim,))
# embedding_layer = TokenAndPositionEmbedding(maxlen, vocab_size, embed_dim)
# x = embedding_layer(inputs)
transformer_block = TransformerBlock(embed_dim, num_heads, ff_dim)
x = transformer_block(inputs)
x = layers.GlobalAveragePooling1D()(x)
x = layers.Dropout(0.1)(x)
x = layers.Dense(20, activation="relu")(x)
x = layers.Dropout(0.1)(x)
outputs = layers.Dense(2, activation="softmax")(x)
model = keras.Model(inputs=inputs, outputs=outputs)
model.compile(optimizer="adam", loss="sparse_categorical_crossentropy", metrics=["accuracy"])


###############################################################################
# Training and validation
###############################################################################

# h5f = h5py.File(data_dir + 'data/dataset_train_all.h5', 'r')
# print(h5f.keys())
# num_idx = len(list(h5f.keys()))//2
# idx_all = np.random.permutation(num_idx)
# idx_train = idx_all[:int(0.9*num_idx)]
# idx_valid = idx_all[int(0.9*num_idx):]

# print("len(idx_train): ", len(idx_train))
# print("idx_train: ", idx_train)
# print("idx_valid: ", idx_valid)


# EPOCH_NUM = 10*len(idx_train)
# start_time = time.time()

# for epoch_num in range(EPOCH_NUM):
    
#     idx = np.random.choice(idx_train)

#     X = h5f['X' + str(idx)][:]
#     Y = h5f['Y' + str(idx)][:]

#     Xc, Yc = clip_datapoints_spliceAI(X, Y, CL, N_GPUS)
#     model_m.fit(Xc, Yc, batch_size=BATCH_SIZE, verbose=1)
    
#     if (epoch_num+1) % len(idx_train) == 0:
#         # Printing metrics (see utils.py for details)
#         model.save('./Models/SpliceAI' + sys.argv[1]
#                    + '_c' + sys.argv[2] + '.h5')

#         print("--------------------------------------------------------------")
#         print("\n\033[1mValidation set metrics:\033[0m")

#         Y_true_1 = [[] for t in range(1)]
#         Y_true_2 = [[] for t in range(1)]
#         Y_pred_1 = [[] for t in range(1)]
#         Y_pred_2 = [[] for t in range(1)]

#         for idx in idx_valid:

#             X = h5f['X' + str(idx)][:]
#             Y = h5f['Y' + str(idx)][:]

#             Xc, Yc = clip_datapoints_spliceAI(X, Y, CL, N_GPUS)
#             Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

#             if not isinstance(Yp, list):
#                 Yp = [Yp]

#             for t in range(1):

#                 is_expr = (Yc[t].sum(axis=(1,2)) >= 1)
#                 print("Yc[t]: ", Yc[t].shape)
#                 print("Yc[t].sum(axis=(1,2)): ", Yc[t].sum(axis=(1,2)))
#                 print("Yc[t].sum(axis=(1,2)): ", Yc[t].sum(axis=(1,2)).shape)
#                 print("Yp[t]: ", Yp[t].shape)
#                 print("Yp[t].sum(axis=(1,2)): ", Yp[t].sum(axis=(1,2)))
#                 print("Yp[t].sum(axis=(1,2)): ", Yp[t].sum(axis=(1,2)).shape)

#                 Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
#                 Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
#                 Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
#                 Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())

#         print ("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                   np.asarray(Y_pred_1[t]))

#         print ("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                   np.asarray(Y_pred_2[t]))

#         print ("\n\033[1mTraining set metrics:\033[0m")

#         Y_true_1 = [[] for t in range(1)]
#         Y_true_2 = [[] for t in range(1)]
#         Y_pred_1 = [[] for t in range(1)]
#         Y_pred_2 = [[] for t in range(1)]

#         for idx in idx_train[:len(idx_valid)]:

#             X = h5f['X' + str(idx)][:]
#             Y = h5f['Y' + str(idx)][:]

#             Xc, Yc = clip_datapoints_spliceAI(X, Y, CL, N_GPUS)
#             Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

#             if not isinstance(Yp, list):
#                 Yp = [Yp]

#             for t in range(1):

#                 is_expr = (Yc[t].sum(axis=(1,2)) >= 1)

#                 Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
#                 Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
#                 Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
#                 Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())

#         print ("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                   np.asarray(Y_pred_1[t]))

#         print ("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                   np.asarray(Y_pred_2[t]))

#         print ("Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
#         print ("--- %s seconds ---" % (time.time() - start_time))
#         start_time = time.time()

#         print ("--------------------------------------------------------------")

#         model.save('./Models/SpliceAI' + sys.argv[1]
#                    + '_c' + sys.argv[2] + '.h5')

#         if (epoch_num+1) >= 6*len(idx_train):
#             kb.set_value(model_m.optimizer.lr,
#                          0.5*kb.get_value(model_m.optimizer.lr))
#             # Learning rate decay

# h5f.close()

###############################################################################