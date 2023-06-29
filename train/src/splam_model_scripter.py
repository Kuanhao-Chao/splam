import torch
from SPLAM import *
from splam_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as pl
import os

MODEL_VERSION = "SPLAM_v11"
device = torch.device("cpu")
model = torch.load("./MODEL/"+MODEL_VERSION+"/splam_14.pt")

# for parameter in model.parameters():
#     print(parameter.shape)

pytorch_total_params = sum(dict((p.data_ptr(), p.numel()) for p in model.parameters()).values())
print("pytorch_total_params: ", pytorch_total_params)

model.to("cpu")
print("model: ", model)
torch.save(model, "./MODEL/"+MODEL_VERSION+"/splam_14_cpu.pt")

# # Trace the model with random data.
# example_input = torch.rand(1, 4, 800, device="cpu")
# model_traced = torch.jit.trace(model, example_input)

# model_traced.save("./MODEL/"+MODEL_VERSION+"/splam_24_scripted.pt")
