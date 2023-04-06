import torch
import torch.nn.functional as F
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Define the PyTorch model
class ConvNet(torch.nn.Module):
    def __init__(self):
        super(ConvNet, self).__init__()
        self.conv_layer = torch.nn.Conv1d(in_channels=64, out_channels=4, kernel_size=3)
        self.pool_layer = torch.nn.MaxPool1d(kernel_size=2)
        
    def forward(self, x):
        x = F.relu(self.conv_layer(x))
        x = self.pool_layer(x)
        x = F.softmax(x, dim=1)
        return x

model = ConvNet()

# # Load the trained weights (assuming they are saved in 'model.pth')
# model.load_state_dict(torch.load('model.pth'))

# Define a function to compute the Grad-CAM heatmap
def compute_grad_cam(model, input):
    # Compute the output scores and gradients
    scores = model(input)
    one_hot = scores.max(1)[1].unsqueeze(1)
    one_hot = torch.zeros_like(scores).scatter_(1, one_hot, 1)
    scores.backward(gradient=one_hot)
    
    # Compute the feature map gradients
    grads = model.conv_layer.weight.grad.cpu().numpy()
    fmap_grads = np.mean(grads, axis=(2))
    
    # Compute the feature maps
    feature_maps = model.conv_layer(input).cpu().detach().numpy()
    
    # Compute the Grad-CAM heatmap
    cam = np.zeros(feature_maps.shape[2:], dtype=np.float32)
    for i, w in enumerate(model.conv_layer.weight.cpu().detach().numpy()):
        cam += w * fmap_grads[i]
    cam = np.maximum(cam, 0)
    cam = cam / np.max(cam)
    
    return cam

# Define the DNA sequence as a one-hot-encoded array
seq = np.array([
    [1, 0, 0, 0],  # A
    [0, 1, 0, 0],  # C
    [0, 0, 1, 0],  # G
    [0, 0, 0, 1],  # T
    [0, 0, 0, 0],  # N
    [1, 0, 0, 0],  # A
    [0, 1, 0, 0],  # C
    [0, 0, 1, 0],  # G
    [0, 0, 0, 1],  # T
    [0, 0, 0, 0]   # N
], dtype=np.float32)

# Convert the one-hot-encoded array to a PyTorch tensor
input_tensor = torch.tensor(seq, requires_grad=True).unsqueeze(0).transpose(1, 2)

# Compute the Grad-CAM heatmap
cam = compute_grad_cam(model, input_tensor)

# Plot the input sequence and the Grad-CAM heatmap
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
ax1.imshow(seq.transpose(), cmap=cm.get_cmap('Set1', 4))
ax1.set_xticks(np.arange(len(seq)))
ax1.set_xticklabels(['A', 'C', 'G', 'T', 'N']*2)
ax1.set_yticks([0, 1, 2, 3])
ax1.set_yticklabels(['A', 'C', 'G', 'T'])

fig.plot()