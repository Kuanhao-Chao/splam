import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

# Define your model that contains a convolutional layer
class MyModel(nn.Module):
    def __init__(self):
        super(MyModel, self).__init__()
        self.conv1 = nn.Conv1d(1, 16, kernel_size=5, stride=1, padding=2)
        self.relu = nn.ReLU()
        self.fc1 = nn.Linear(16 * 50, 10)

    def forward(self, x):
        x = self.conv1(x)
        x = self.relu(x)
        x = x.view(x.size(0), -1)
        x = self.fc1(x)
        return x

# Load your pre-trained model
model = MyModel()
# model.load_state_dict(torch.load('model.pt'))

# Define the hook function to retrieve the gradients from the last convolutional layer
grads = None
def save_grads(grad):
    global grads
    grads = grad

def get_last_conv_layer(model):
    for m in model.modules():
        if isinstance(m, nn.Conv1d):
            m.register_backward_hook(save_grads)
    return m

last_conv_layer = get_last_conv_layer(model)

# Load your input sequence and create a PyTorch tensor
input_seq = np.random.rand(1, 1, 100) # shape is (batch_size, channels, sequence_length)
x = torch.tensor(input_seq, dtype=torch.float)

# Pass the tensor through the model to generate a prediction
output = model(x)

# Compute the gradients of the output with respect to the convolutional layer
output[:, output.argmax()].backward()

# Compute the weights of the activation map using the gradients
alpha_k = grads.mean(dim=(2, 3), keepdim=True)
cam_weights = (alpha_k * last_conv_layer.weight).sum(dim=1, keepdim=True)

# Multiply the weights by the output of the last convolutional layer to obtain the Grad-CAM
grad_cam = (cam_weights * last_conv_layer(x)).sum(dim=1, keepdim=True)

# Visualize the Grad-CAM
plt.imshow(grad_cam[0, 0].detach().numpy(), cmap='jet')
plt.show()