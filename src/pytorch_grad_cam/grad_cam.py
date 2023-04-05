import numpy as np
from pytorch_grad_cam.base_cam import BaseCAM


class GradCAM(BaseCAM):
    def __init__(self, model, target_layers, use_cuda=False,
                 reshape_transform=None):
        super(
            GradCAM,
            self).__init__(
            model,
            target_layers,
            use_cuda,
            reshape_transform)
        print("Constructor!!")

    def get_cam_weights(self,
                        input_tensor,
                        target_layer,
                        target_category,
                        activations,
                        grads):
        print("grads (get_cam_weights): ", grads)
        print("grads (get_cam_weights): ", grads.shape)
        return np.mean(grads, axis=(1, 2))
        # print("grads: ", grads.shape)
        # return np.mean(grads, axis=(0))
