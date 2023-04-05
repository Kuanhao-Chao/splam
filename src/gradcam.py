from pytorch_grad_cam import GradCAM, HiResCAM, ScoreCAM, GradCAMPlusPlus, AblationCAM, XGradCAM, EigenCAM, FullGrad
from pytorch_grad_cam.utils.model_targets import ClassifierOutputTarget
from pytorch_grad_cam.utils.image import show_cam_on_image
from torchvision.models import resnet50
import cv2
import numpy as np
from pytorch_grad_cam import GuidedBackpropReLUModel
from pytorch_grad_cam.utils.image import show_cam_on_image, \
    deprocess_image, \
    preprocess_image
from pytorch_grad_cam.utils.model_targets import ClassifierOutputTarget, SoftmaxOutputTarget, SPLAMOutputTarget



import torch
from torchvision import models
from pytorch_grad_cam import GradCAM, \
    HiResCAM, \
    ScoreCAM, \
    GradCAMPlusPlus, \
    AblationCAM, \
    XGradCAM, \
    EigenCAM, \
    EigenGradCAM, \
    LayerCAM, \
    FullGrad, \
    GradCAMElementWise

import torch.nn as nn
from splam_dataset_Chromsome import *
from splam_utils import *
import matplotlib.pyplot as plt
from pytorch_grad_cam.activations_and_gradients import ActivationsAndGradients




# model = resnet50(pretrained=True)
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
input_tensor = DNAs



target_layers = [model.last_cov]
# image_path = "../../pytorch-grad-cam/examples/cars_segmentation.png"

# rgb_img = cv2.imread(image_path, 1)[:, :, ::-1]
# rgb_img = np.float32(rgb_img) / 255

# print("rgb_img: ", rgb_img.shape)

# input_tensor = preprocess_image(rgb_img,
#                                 mean=[0.485, 0.456, 0.406],
#                                 std=[0.229, 0.224, 0.225])

print("input_tensor: ", input_tensor.shape)



# Note: input_tensor can be a batch tensor with several images!


method = ['gradcam', 'hirescam', 'gradcam++',
                                 'scorecam', 'xgradcam',
                                 'ablationcam', 'eigencam',
                                 'eigengradcam', 'layercam', 'fullgrad']
method = method[-1]
# Construct the CAM object once, and then re-use it on many images:
cam = GradCAM(model=model, target_layers=target_layers, use_cuda=False)

print("cam: ", cam)
# You can also use it within a with statement, to make sure it is freed,
# In case you need to re-create it inside an outer loop:
# with GradCAM(model=model, target_layers=target_layers, use_cuda=args.use_cuda) as cam:
#   ...

# We have to specify the target we want to generate
# the Class Activation Maps for.
# If targets is None, the highest scoring category
# will be used for every image in the batch.
# Here we use ClassifierOutputTarget, but you can define your own custom targets
# That are, for example, combinations of categories, or specific outputs in a non standard model.

targets = [SPLAMOutputTarget()]
# [SoftmaxOutputTarget()]

# input_tensor = input_tensor.detach().numpy()

# input_tensor = np.repeat(input_tensor[:, :, :, np.newaxis], 800, axis=3)
# input_tensor = input_tensor[:, :3, : , :]
# input_tensor = np.float32(input_tensor) / 255
# input_tensor = torch.tensor(input_tensor)
# print("\tinput_tensor: ", input_tensor.shape)

# You can also pass aug_smooth=True and eigen_smooth=True, to apply smoothing.
grayscale_cam = cam(input_tensor=input_tensor, targets=targets)






# In this example grayscale_cam has only one image in the batch:
grayscale_cam = grayscale_cam[0, :]
print("\tgrayscale_cam: ", grayscale_cam.shape)
print("\tgrayscale_cam: ", grayscale_cam)

input_tensor = input_tensor.detach().numpy()

print("input_tensor: ", input_tensor.shape)
print("grayscale_cam: ", grayscale_cam.shape)

# # input_tensor = np.repeat(input_tensor, [1, 1, 2], axis=0)

# input_tensor = np.transpose(input_tensor, (0, 2, 3, 1))
# print("input_tensor: ", input_tensor.shape)
# # print("input_tensor: ", input_tensor)

# # grayscale_cam = np.repeat(grayscale_cam, [1, 1, 2], axis=0)
# grayscale_cam = np.repeat(grayscale_cam[:, :, np.newaxis], 800, axis=2)
# grayscale_cam = grayscale_cam[:3, : , :]
# grayscale_cam = np.float32(grayscale_cam) / 255
# grayscale_cam = np.transpose(grayscale_cam, (1, 2, 0))
# print("grayscale_cam: ", grayscale_cam.shape)
# grayscale_cam = np.ones(grayscale_cam.shape)





# print("input_tensor : ", input_tensor.shape)
# print("grayscale_cam: ", grayscale_cam.shape)

# visualization = show_cam_on_image(input_tensor[0], grayscale_cam, use_rgb=True)
# print("visualization: ", visualization.shape)
# print("visualization: ", visualization)

# # cam_image is RGB encoded whereas "cv2.imwrite" requires BGR encoding.
# cam_image = cv2.cvtColor(visualization, cv2.COLOR_RGB2BGR)




# # # print("cam_image: ", cam_image.shape)

# # gb_model = GuidedBackpropReLUModel(model=model, use_cuda=False)

# # input_tensor = np.transpose(input_tensor, (0, 3, 1, 2))
# # print("input_tensor: ", input_tensor.shape)

# # gb = gb_model(torch.tensor(input_tensor[0]), target_category=None)

# # cam_mask = cv2.merge([grayscale_cam, grayscale_cam, grayscale_cam])
# # cam_gb = deprocess_image(cam_mask * gb)
# # gb = deprocess_image(gb)

# cv2.imwrite(f'{method}_cam.jpg', cam_image)
# # cv2.imwrite(f'{method}_gb.jpg', gb)
# # cv2.imwrite(f'{method}_cam_gb.jpg', cam_gb)