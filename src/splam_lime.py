import torch
import numpy as np
import lime
from lime import lime_tabular
import torch.nn as nn
from splam_dataset_Chromsome import *
import tensorflow as tf
import matplotlib.pyplot as plt

def parse_junction(name):
    # print("name: ", name)
    res = name.split(":")
    strand = name[-2]
    chr_name = res[0]
    if strand == "+":
        start = int(res[1].split("-")[0])+200
        end = int(res[2].split("-")[1].split('(')[0])-200
    elif strand == "-":
        start = int(res[2].split("-")[0])+200
        end = int(res[1].split("-")[1].split('(')[0])-200
    # print("start: ", start)
    # print("end  : ", end)
    return (chr_name, start, end, strand)

# class model_wrapper:
#     def __init__(self,model):
            
#             self.model = model
 
#     def predict(self,input_data):        
#         self.pred = self.model.predict(input_data).ravel()
#         return np.array([[1-self.pred,self.pred]]).T[:,:,0]
    
    




# Defining the Grad-CAM algorithm
def grad_cam(layer_name, data):
    grad_model = tf.keras.models.Model(
        [model.inputs], [model.get_layer(layer_name).output, model.output]
    )
    last_conv_layer_output, preds = grad_model(data)
    
    with tf.GradientTape() as tape:
        last_conv_layer_output, preds = grad_model(data)
        pred_index = tf.argmax(preds[0])
        class_channel = preds[:, pred_index]
        
    grads = tape.gradient(class_channel, last_conv_layer_output)
    
    pooled_grads = tf.reduce_mean(grads, axis=(0))
    
    last_conv_layer_output = last_conv_layer_output[0]
    
    heatmap = last_conv_layer_output * pooled_grads
    heatmap = tf.reduce_mean(heatmap, axis=(1))
    heatmap = np.expand_dims(heatmap,0)
    return heatmap


model_path = "../src/MODEL/SPLAM_v11/splam_14.pt"
model = torch.load(model_path)

# Wrap model
# wrapped_mod = model_wrapper(model)

# explainer = lime_tabular.RecurrentTabularExplainer(X_dev,training_labels=tf.keras.utils.to_categorical(y_dev), feature_names=["Inter beat int"],
#                                                    discretize_continuous=False, feature_selection='auto', class_names=['Healthy','Sick'])


BATCH_SIZE = 100
MODEL_VERSION = "SPLAM_v10/"
N_WORKERS = 1
shuffle=True
target = "pos"
device="mps"
criterion = nn.BCELoss()
# TARGETS = ["pos", "pos_refseq_protein_alts", "neg_1", "neg_random"]


test_loader = get_eval_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS, shuffle, target)


# Explaining the predictions don on the test set. Here we only look at sequences that are classified as abnormal/sick
label = ["healthy", "sick"]
cnt = 0





# Class activation map from the input layer to the last Conv. layer
layer_name = "last_conv"
cnt = 0
# for i in X_test:
#     data = np.expand_dims(i,0)
#     pred = model.predict(data)[0][0]
#     if  pred > 0.5:
#         heatmap = grad_cam(layer_name,data)
#         print(f"Model prediction = sick ({pred}), True label = {label[int(y_test[cnt])]}")
#         plt.figure(figsize=(30,4))
#         plt.imshow(np.expand_dims(heatmap,axis=2),cmap='Reds', aspect="auto", interpolation='nearest',extent=[0,300,i.min(),i.max()], alpha=0.5)
#         plt.plot(i,'k')
#         plt.colorbar()
#         plt.show()
#     cnt +=1




model.eval()
for batch_idx, data in enumerate(test_loader):
    # print("batch_idx: ", batch_idx)
    # DNAs:  torch.Size([40, 800, 4])
    # labels:  torch.Size([40, 1, 800, 3])
    DNAs, labels, chr = data 

    # print("chr: ", chr)

    junc_name = map(parse_junction, chr)
    junc_name = list(junc_name)
    # print(chr.splt(":"))
    DNAs = DNAs.to(torch.float32).to(device)
    labels = labels.to(torch.float32).to(device)

    DNAs = torch.permute(DNAs, (0, 2, 1))
    labels = torch.permute(labels, (0, 2, 1))
    loss, yp = model_fn(DNAs, labels, model, criterion)

    #######################################
    # predicting all bp.
    #######################################    
    is_expr = (labels.sum(axis=(1,2)) >= 1)
    # print("is_expr: ", is_expr)

    # Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
    Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
    Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
    Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
    Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

    A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
    A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
    D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
    D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()
    print("A_YP: ", A_YP)

    # data = np.expand_dims(i,0)
    # pred = model.predict(data)[0][0]
    # if  pred > 0.5:
    heatmap = grad_cam(layer_name,DNAs)
    print(f"Model prediction = sick ({A_YP}), True label = {1}")
    plt.figure(figsize=(30,4))
    plt.imshow(np.expand_dims(heatmap,axis=2),cmap='Reds', aspect="auto", interpolation='nearest',extent=[0,300,DNAs.min(),DNAs.max()], alpha=0.5)
    plt.plot(i,'k')
    plt.colorbar()
    plt.show()
    cnt +=1



    # pred = model.predict(np.expand_dims(j,axis=0))[0][0]


    # if  pred > 0.5: 
    #     exp = explainer.explain_instance(np.expand_dims(j,axis=0), wrapped_mod.predict, num_features=300)
    #     explanations = exp.as_list()
    #     heatmap = np.zeros([1,300])
    #     for k in explanations:
    #         if k[1] > 0:
    #             heatmap[0,int(k[0].split("-")[-1])] = k[1]
    #     print(f"Model prediction = sick ({pred}), True label = {label[int(y_test[cnt])]}")
    #     plt.figure(figsize=(30,4))
    #     plt.imshow(np.expand_dims(heatmap,axis=2),cmap='Reds', aspect="auto", interpolation='nearest',extent=[0,300,j.min(),j.max()], alpha=0.5)
    #     plt.plot(j,'k', label="V_d")
    #     plt.colorbar()
    #     plt.show()
    # cnt +=1
