import torch
from SPLAM import *
from splam_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as pl
import os

device = torch.device("cpu")
model = torch.load("./MODEL/SPLAM_v1/splam_24.pt")
model.to("cpu")
print("model: ", model)

# Trace the model with random data.
example_input = torch.rand(1, 4, 800, device="cpu")
model_traced = torch.jit.trace(model, example_input)

model_traced.save("./MODEL/SPLAM_v1/splam_24_scripted.pt")