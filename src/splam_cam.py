import argparse
import cv2
import numpy as np
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



from pytorch_grad_cam import GuidedBackpropReLUModel
from pytorch_grad_cam.utils.image import show_cam_on_image, \
    deprocess_image, \
    preprocess_image
from pytorch_grad_cam.utils.model_targets import ClassifierOutputTarget, ClassifierOutputSoftmaxTarget

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

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--use-cuda', action='store_true', default=False,
                        help='Use NVIDIA GPU acceleration')
    parser.add_argument(
        '--image-path',
        type=str,
        default='./examples/both.png',
        help='Input image path')
    parser.add_argument('--aug_smooth', action='store_true',
                        help='Apply test time augmentation to smooth the CAM')
    parser.add_argument(
        '--eigen_smooth',
        action='store_true',
        help='Reduce noise by taking the first principle componenet'
        'of cam_weights*activations')
    parser.add_argument('--method', type=str, default='gradcam',
                        choices=['gradcam', 'hirescam', 'gradcam++',
                                 'scorecam', 'xgradcam',
                                 'ablationcam', 'eigencam',
                                 'eigengradcam', 'layercam', 'fullgrad'],
                        help='Can be gradcam/gradcam++/scorecam/xgradcam'
                             '/ablationcam/eigencam/eigengradcam/layercam')

    args = parser.parse_args()
    args.use_cuda = args.use_cuda and torch.cuda.is_available()
    if args.use_cuda:
        print('Using GPU for acceleration')
    else:
        print('Using CPU for computation')

    return args


# def get_cam_weights(input_tensor: torch.Tensor,
#                     target_layers: List[torch.nn.Module],
#                     targets: List[torch.nn.Module],
#                     activations: torch.Tensor,
#                     grads: torch.Tensor) -> np.ndarray:
#     raise Exception("Not Implemented")

# function to extract grad




# Define the hook function to retrieve the gradients from the last convolutional layer
grads = None
def save_grads(module: nn.Module, grad_input, grad_output):
    print("inside save_grads")
    global grads
    grads = grad_output[0]


def get_last_conv_layer(model):
    m = model.last_cov
    m.register_full_backward_hook(save_grads)
    # for m in model.modules():
    #     if isinstance(m, nn.Conv1d):
    #         m.register_backward_hook(save_grads)
    return m



if __name__ == '__main__':
    """ python cam.py -image-path <path_to_image>
    Example usage of loading an image, and computing:
        1. CAM
        2. Guided Back Propagation
        3. Combining both
    """

    args = get_args()
    methods = \
        {"gradcam": GradCAM,
         "hirescam": HiResCAM,
         "scorecam": ScoreCAM,
         "gradcam++": GradCAMPlusPlus,
         "ablationcam": AblationCAM,
         "xgradcam": XGradCAM,
         "eigencam": EigenCAM,
         "eigengradcam": EigenGradCAM,
         "layercam": LayerCAM,
         "fullgrad": FullGrad,
         "gradcamelementwise": GradCAMElementWise}

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

    last_conv_layer = get_last_conv_layer(model)
    print("last_conv_layer: ", last_conv_layer)

    print("grads: ", grads)

    # # Step 5: Pass the preprocessed image through the model and get the output
    # yp = model(DNAs)
    # loss = categorical_crossentropy_2d(labels, yp, criterion)
    # # yp[:, yp.argmax()].backward()
    # loss.backward()
    #     # retain_graph=True)
    # print("loss: ", loss)
    # print("yp: ", yp)
    # print("yp: ", yp.shape)

    # print("grads: ", grads)
    # print("grads: ", grads.shape)
    # print("grads: ", len(grads))

    # # print(last_conv_layer.weight)
    # print("last_conv_layer.weight.shape: ", last_conv_layer.weight.shape)


    # print("last_conv_layer.weight.grad.shape: ", last_conv_layer.weight.grad.shape)







    # feature_maps = last_conv_layer.forward(DNAs)
    # print("feature_maps: ", feature_maps)

    # # cam_weights = (grads * last_conv_layer.weight).sum(dim=1, keepdim=True)

    # # print("cam_weights: ", cam_weights)
    # # print("cam_weights: ", cam_weights.shape)


    # # alpha_k = torch.mean(grads, dim=(0), keepdim=True)
    # # print("alpha_k: ", alpha_k)
    # # print("alpha_k: ", alpha_k.shape)

    # grads_output_tmp = grads * model(DNAs)

    # grads_output = torch.mean(grads_output_tmp, dim=(0,1))
    # # , keepdim=True)

    # print("grads_output: ", grads_output)
    # print("grads_output: ", grads_output.shape)

    # print("model(DNAs).shape: ", grads_output.shape)
    # # print("model(DNAs).shape: ", grads_output[0][0].shape)

    # plt.plot([*range(800)], grads_output.detach().numpy())
    # plt.show()


    # last_conv_layer.weight = torch.permute(last_conv_layer.weight, (0, 2, 1))
    # Getting the last convolutional layer of the model.




    # # yp[:, yp.argmax()].backward()
    # # Compute the weights of the activation map using the gradients
    # alpha_k = grads.mean(dim=(2), keepdim=True)
    # alpha_k = alpha_k.reshape([8, 3, 1, 1])
    # print("alpha_k: ", alpha_k)
    # print("alpha_k: ", alpha_k.shape)

    # alpha_k = torch.permute(alpha_k, (0, 2, 1, 3))

    # cam_weights = (alpha_k * last_conv_layer.weight).sum(dim=1, keepdim=True)

    # print("alpha_k: ", alpha_k)
    # print("cam_weights: ", cam_weights)
    # print("cam_weights.shape: ", cam_weights.shape)


    # # Calculate gradients
    # one_hot_output = torch.zeros((1, yp.size()[-1]), dtype=torch.float32, device=device)
    # one_hot_output[0][target_class] = 1
    # output.backward(gradient=one_hot_output)

    # # Get gradients and last convolutional layer output
    # gradients = input_tensor.grad.cpu().numpy()[0]
    # last_conv_layer_output = model._modules['features'][-1](input_tensor).cpu().detach().numpy()[0]





    # # Step 6: Calculate the gradients of the output with respect to the last convolutional layer
    # grads = torch.autograd.grad(outputs=yp, inputs=model.last_cov.weight, grad_outputs=torch.ones(yp.size()), create_graph=True)[0]

    # print("grads: ", grads)
    # print("grads: ", grads.shape)


    # # Step 7: Calculate the weights of each feature map in the last convolutional layer
    # weights = torch.mean(grads, axis=(1, 2)).squeeze()

    # print("weights: ", weights)
    # print("weights: ", weights.shape)

    # # Step 8: Compute the class activation map
    # cam = torch.zeros((weights.size(0), 800))

    # print("cam: ", cam)
    # print("cam: ", cam.shape)

    # for i, w in enumerate(weights):
    #     cam[i] = w * model.last_cov.weight[i]

    # cam = torch.sum(cam, axis=0)
    # cam = np.maximum(cam.detach().numpy(), 0)


    # print("cam: ", cam)
    # print("cam: ", cam.shape)



    # last_cov_grad = model.last_cov.weight.grad
    # print("last_cov_grad: ", last_cov_grad)
    # print("last_cov_grad: ", last_cov_grad.shape)


    # last_cov = model.last_cov

    # for name, para in last_cov.named_parameters():
    #     print("name: ", name)
    #     print("para: ", para)
    #     # if 'weight' in name:
    #     #     last_grad_norm += para.grad.norm().cpu().item() 
    # # print("last grad norm", last_grad_norm) 


    # DNAs_grad = DNAs.grad
    # print("DNAs_grad :",  DNAs_grad)
    # print("DNAs_grad.shape :",  DNAs_grad.shape)
    # heatmap = torch.mean(DNAs_grad, dim=1).squeeze()
    # print("heatmap :",  heatmap)
    # print("heatmap :",  heatmap.shape)

    # plt.plot([*range(800)], heatmap[0])
    # plt.show()



    # reshape_transform = None
    # eigen_smooth = False
    # activations_and_grads = ActivationsAndGradients(
    #     model, target_layers, reshape_transform)

    # print("activations_and_grads: ", activations_and_grads)

    # activations_list = [a.cpu().data.numpy()
    #                     for a in activations_and_grads.activations]
    # grads_list = [g.cpu().data.numpy()
    #                 for g in activations_and_grads.gradients]
    
    # print("activations_list: ", activations_list)
    # print("grads_list: ", grads_list)



    # target_size = input_tensor.size(-1), input_tensor.size(-2)
    # print(target_size)

    # cam_per_target_layer = []
    # # Loop over the saliency image from every layer
    # for i in range(len(target_layers)):
    #     target_layer = target_layers[i]
    #     layer_activations = None
    #     layer_grads = None
    #     if i < len(activations_list):
    #         layer_activations = activations_list[i]
    #     if i < len(grads_list):
    #         layer_grads = grads_list[i]

    #     print("layer_grads: ", layer_grads)

    #     # cam = get_cam_image(input_tensor,
    #     #                             target_layer,
    #     #                             targets,
    #     #                             layer_activations,
    #     #                             layer_grads,
    #     #                             eigen_smooth)
    #     # weights = get_cam_weights(input_tensor,
    #     #                                target_layer,
    #     #                                targets,
    #     #                                activations,
    #     #                                grads)
    #     # weighted_activations = weights[:, :, None, None] * activations
    #     # if eigen_smooth:
    #     #     cam = get_2d_projection(weighted_activations)
    #     # else:
    #     #     cam = weighted_activations.sum(axis=1)
    #     # return cam

    #     cam = np.maximum(cam, 0)
    #     scaled = scale_cam_image(cam, target_size)
    #     cam_per_target_layer.append(scaled[:, None, :])


    # # cam_per_layer = compute_cam_per_layer(input_tensor,
    # #                                             targets,
    # #                                             True)







    # # targets = [ClassifierOutputSoftmaxTarget(3)]
    # # for target, output in zip(targets, outputs):
    # #     print("target: ", target)
    # #     print("output: ", output.shape)
    # #     print("target(output): ", target(output))
    # # loss = sum([target(output)
    # #             for target, output in zip(targets, outputs)])   
    # # # loss = sum(loss)
    # # print("loss: ", loss)
    # # loss.backward()


    # # with cam_algorithm(model=model,
    # #                    target_layers=target_layers,
    # #                    use_cuda=args.use_cuda) as cam:

    # #     # AblationCAM and ScoreCAM have batched implementations.
    # #     # You can override the internal batch size for faster computation.
    # #     cam.batch_size = BATCH_SIZE


    # #     grayscale_cam = cam(input_tensor=input_tensor,
    # #                         targets=targets,
    # #                         aug_smooth=args.aug_smooth,
    # #                         eigen_smooth=args.eigen_smooth)

    # #     # Here grayscale_cam has only one image in the batch
    # #     grayscale_cam = grayscale_cam[0, :]

    # #     cam_image = show_cam_on_image(rgb_img, grayscale_cam, use_rgb=True)

    # #     # cam_image is RGB encoded whereas "cv2.imwrite" requires BGR encoding.
    # #     cam_image = cv2.cvtColor(cam_image, cv2.COLOR_RGB2BGR)

    # # gb_model = GuidedBackpropReLUModel(model=model, use_cuda=args.use_cuda)
    # # gb = gb_model(input_tensor, target_category=None)

    # # cam_mask = cv2.merge([grayscale_cam, grayscale_cam, grayscale_cam])
    # # cam_gb = deprocess_image(cam_mask * gb)
    # # gb = deprocess_image(gb)

    # # cv2.imwrite(f'{args.method}_cam.jpg', cam_image)
    # # cv2.imwrite(f'{args.method}_gb.jpg', gb)
    # # cv2.imwrite(f'{args.method}_cam_gb.jpg', cam_gb)
