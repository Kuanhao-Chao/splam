###############################################################################
# This file has the functions necessary to create the SpliceAI model.
###############################################################################

from keras.models import Model
from keras.layers import Input
from keras.layers.core import Activation
from keras.layers.convolutional import Conv1D
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import add, Cropping1D
import keras.backend as kb
import numpy as np
import tensorflow as tf


def ResidualUnit(l, w, ar):
    # Residual unit proposed in "Identity mappings in Deep Residual Networks"
    # by He et al.

    def f(input_node):

        bn1 = BatchNormalization()(input_node)
        act1 = Activation('relu')(bn1)
        conv1 = Conv1D(filters=l, kernel_size=(w,), dilation_rate=(ar,), padding='same')(act1)
        bn2 = BatchNormalization()(conv1)
        act2 = Activation('relu')(bn2)
        conv2 = Conv1D(filters=l, kernel_size=(w,), dilation_rate=(ar,), padding='same')(act2)
        output_node = add([conv2, input_node])

        return output_node

    return f


def SpliceAI(L, W, AR):
    # L: Number of convolution kernels
    # W: Convolution window size in each residual unit
    # AR: Atrous rate in each residual unit

    assert len(W) == len(AR)

    CL = 2 * np.sum(AR*(W-1))

    input0 = Input(shape=(None, 4))
    conv = Conv1D(filters=L, kernel_size=(1,))(input0)
    skip = Conv1D(filters=L, kernel_size=(1,))(conv)

    for i in range(len(W)):
        conv = ResidualUnit(L, W[i], AR[i])(conv)
        
        if (((i+1) % 4 == 0) or ((i+1) == len(W))):
            # Skip connections to the output after every 4 residual units
            dense = Conv1D(filters=L, kernel_size=(1,))(conv)
            skip = add([skip, dense])

    skip = Cropping1D(cropping=(int(CL/2), int(CL/2)))(skip)

    output0 = [[] for t in range(1)]

    for t in range(1):
        output0[t] = Conv1D(filters=3, kernel_size=(1,), activation='softmax')(skip)
    
    model = Model(inputs=input0, outputs=output0)

    print("output0: ", output0)
    # print("output0.size(): ", output0.size())

    return model


def categorical_crossentropy_2d(y_true, y_pred):
    # Standard categorical cross entropy for sequence outputs
    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.cast(y_pred, tf.float32)
    # print("y_true: ", y_true)
    # print("y_pred: ", y_pred)

    return - kb.mean(y_true[:, :, 0]*kb.log(y_pred[:, :, 0]+1e-10)
                   + y_true[:, :, 1]*kb.log(y_pred[:, :, 1]+1e-10)
                   + y_true[:, :, 2]*kb.log(y_pred[:, :, 2]+1e-10))
