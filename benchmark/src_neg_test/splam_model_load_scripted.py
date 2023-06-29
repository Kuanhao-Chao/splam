import torch
from SPLAM import *
from splam_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as pl
import os

device = torch.device("cpu")

# model = torch.jit.load("./MODEL/SPLAM_v11/splam_14.pt")
# model = model.to("mps")

model = torch.load("./MODEL/SPLAM_v11/splam_14.pt")
model.to("cpu")

# Trace the model with random data.
example_input = torch.rand(1, 4, 800, device="cpu")
model_traced = torch.jit.trace(model, example_input)
# model_scripted = torch.jit.script(model_traced)

model_traced.save("./MODEL/SPLAM_v11/splam_14_script.pt")
# model_scripted.save("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_v18/splam_19_scripted.pt")
