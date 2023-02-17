import torch
from SpliceNN import *
from SpliceNN_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as pl
import os

device = torch.device("cpu")

model = torch.jit.load("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24_scripted.pt")

model = model.to("mps")

# model = torch.load("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24.pt")
# model.to("cpu")
# print("model: ", model)

# # print(model)
# # torch.save(model, "./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24_cpu.pt")

# # example_input = torch.rand(1, 4, 800, device="cpu")
# # model_scripted = torch.jit.script(model, example_input) # Export to TorchScript
# # model_scripted.save("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24_scripted.pt") # Save

# # Trace the model with random data.
# example_input = torch.rand(1, 4, 800, device="cpu")
# model_traced = torch.jit.trace(model, example_input)
# # model_scripted = torch.jit.script(model_traced)

# model_traced.save("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24_scripted.pt")
# # model_scripted.save("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_v18/SpliceNN_19_scripted.pt")