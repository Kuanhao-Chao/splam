import matplotlib.pyplot as plt
import numpy as np
import torch.nn as nn
from splam_dataset_Chromsome import *
from splam_utils import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from pytorch_grad_cam.activations_and_gradients import ActivationsAndGradients

import torch.nn.functional as F
from torch.autograd import Variable

def get_activations_and_gradients(model, input):
    # Set model to evaluation mode
    model.eval()

    # Define a hook function to get the activations and gradients
    activations = []
    gradients = []

    def save_activations_hook(module, input, output):
        activations.append(output)

    def save_gradients_hook(module, grad_input, grad_output):
        gradients.append(grad_output[0])

    # Register the hooks to the last convolution layer
    # print(model)
    layer = model.last_cov
    handle_activations = layer.register_forward_hook(save_activations_hook)
    handle_gradients = layer.register_backward_hook(save_gradients_hook)

    # Forward pass
    input_var = Variable(input.unsqueeze(0), requires_grad=True)
    input_var = torch.permute(input_var, (0, 2, 1))
    print("input_var: ", input_var.shape)

    outputs = model(input_var)
    print("output   : ", outputs.shape)

    y_label = np.zeros(outputs.size())
    y_label[:,:,200] = 1
    y_label[:,:,600] = 1

    y_label = torch.tensor(y_label)


    # This is cross entropy loss
    loss = - torch.mean(y_label[:, 0, :]*torch.log(outputs[:, 0, :]+1e-10)
                        + y_label[:, 1, :]*torch.log(outputs[:, 1, :]+1e-10)
                        + y_label[:, 2, :]*torch.log(outputs[:, 2, :]+1e-10))

    # Backward pass to compute gradients
    model.zero_grad()
    loss.backward()

    # Remove the hooks
    handle_activations.remove()
    handle_gradients.remove()

    # Return the activations and gradients as tensors
    activations = activations[0].detach()
    gradients = gradients[0].detach()
    return activations, gradients


def compute_grad_cam(model, input, class_idx=None):
    print("input: ", input.shape)
    # Get the activations and gradients for the input image
    activations, gradients = get_activations_and_gradients(model, input)
    print("activations: ", activations)

    # if all(v for row in gradients for val in row for v in val):
    #     print("Array contains only True values")
    # else:
    #     print("Array contains at least one False value")

    print("gradients  : ", gradients)
    print(gradients == 0)

    # Compute the weights using the global average pooling of gradients
    weights = F.adaptive_avg_pool2d(gradients, 1)
    print("weights: ", weights)

    # Compute the class activation map as the weighted sum of activations
    cam = torch.sum(weights * activations, dim=1, keepdim=True)
    print("cam: ", cam)
    print("cam: ", cam.shape)

    # Apply ReLU and resize to the input size
    # cam = F.relu(cam)
    # cam = F.interpolate(cam, size=input.shape[1:], mode='bilinear', align_corners=False)

    # Normalize the heatmap to [0, 1] range
    cam -= cam.min()
    cam /= cam.max()
    print("cam: ", cam)
    print("cam: ", cam[0,0])

    # # If class_idx is not given, use the predicted class
    # if class_idx is None:
    #     class_idx = torch.argmax(model(input.unsqueeze(0)))

    # Convert the heatmap to a numpy array and return
    return cam[0, 0].numpy()


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

# Example one-hot-encoded DNA sequences
# dna_array = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

model_path = "../src/MODEL/SPLAM_v11/splam_14.pt"
model = torch.load(model_path)


BATCH_SIZE = 200
MODEL_VERSION = "SPLAM_v10/"
N_WORKERS = 1
shuffle=True
target = "pos"
device="cpu"
criterion = nn.BCELoss()
# TARGETS = ["pos", "pos_refseq_protein_alts", "neg_1", "neg_random"]


test_loader = get_eval_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS, shuffle, target)

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


model.to(device)
model.eval()

for batch_idx, data in enumerate(test_loader):
    # print("batch_idx: ", batch_idx)
    # DNAs:  torch.Size([40, 800, 4])
    # labels:  torch.Size([40, 1, 800, 3])
    DNAs, labels, chr = data 

    DNAs = torch.permute(DNAs, (0, 2, 1))
    DNAs.requires_grad = True
    labels = torch.permute(labels, (0, 2, 1))
    # print("chr: ", chr)

    junc_name = map(parse_junction, chr)
    junc_name = list(junc_name)
    # print(chr.splt(":"))
    DNAs = DNAs.to(torch.float32).to(device)
    labels = labels.to(torch.float32).to(device)
    print("DNAs: ", DNAs[0].shape)
    # plt.figure(figsize=(10,2))
    # plt.plot(DNAs[0])
    # # plt.show()
    # plt.savefig("./DNA_plot.png", dpi=300)
    break
    

print("DNAs: ", DNAs.shape)
DNAs = torch.permute(DNAs, (0, 2, 1))
dna_array = DNAs[0].detach().numpy()
print("DNAs: ", DNAs.shape)
# dna_array = np.repeat(dna_array, [1,0], axis=0)




# Create dictionary mapping nucleotide letters to colors
nucleotide_colors = {'A': 'red', 'C': 'blue', 'G': 'green', 'T': 'orange', 'N': 'gray'}

# Calculate the width of each bar to fill the plot
bar_width = 1.0

# Loop over each position in the DNA sequence
for i in range(dna_array.shape[0]):

    # Get the index of the non-zero element in the one-hot-encoded vector
    nucleotide_index = np.argmax(dna_array[i])

    # If the nucleotide is 'N', check if the sum of the one-hot-encoded vector is zero, and set the nucleotide index to 4
    if nucleotide_index == 0 and np.sum(dna_array[i]) == 0:
        nucleotide_index = 4

    # Get the corresponding nucleotide letter
    nucleotide = 'ACGTN'[nucleotide_index]

    # Get the color for this nucleotide
    nucleotide_color = nucleotide_colors[nucleotide]

    # Plot the nucleotide count as a bar chart
    plt.bar(i, 1, color=nucleotide_color, width=bar_width)
    # , edgecolor='pink')



# Compute the Grad-CAM heatmap
heatmap = compute_grad_cam(model, torch.tensor(dna_array, dtype=torch.float32))

# Set the figure size
fig = plt.figure(figsize=(12, 8))

print("heatmap: ", heatmap.shape)
print("heatmap: ", heatmap)

heatmap = heatmap.reshape((1,800))
# Plot the heatmap with linear axis and increased bar height
plt.imshow(heatmap, cmap='coolwarm', aspect='auto', interpolation='nearest', vmin=0, vmax=1)
plt.colorbar()

# Show the plot
plt.show()

# print("heatmap: ", heatmap)
# # Visualize the heatmap overlaid on the input image
# fig, ax = plt.subplots()
# ax.imshow(dna_array.T, cmap=cm.binary)

# heatmap = heatmap / heatmap.max()  # Normalize heatmap to [0, 1]
# heatmap = cm.jet(heatmap)  # Apply color map

# heatmap = (heatmap[:, :, :3] * 255).astype(np.uint8)  # Convert to RGB image
# heatmap = np.flip(heatmap, axis=0)  # Flip image vertically

# ax.imshow(heatmap, alpha=0.5, extent=[0, dna_array.shape[0], 0, 4], aspect='auto')
# plt.show()




# # Set the x-axis and y-axis labels and limits
# plt.xlim(-0.5, dna_array.shape[0] - 0.5)
# plt.ylim(0, 1)
# plt.xlabel('Position')
# plt.ylabel('Nucleotide')

# # Show the plot
# plt.show()