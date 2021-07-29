import numpy as np
import torch
import torch.nn as nn

from KalmanNet_data import N_T

MSE_test_linear_arr = np.empty([N_T])

def NNTest(SysModel, test_input, test_target):

    # MSE LOSS Function
    loss_fn = nn.MSELoss(reduction='mean')

    Model = torch.load('best-model.pt')

    Model.eval()
    torch.no_grad()

    for j in range(0, N_T):
        # Unrolling Forward Pass

        y_mdl_tst = test_input[j, :, :]

        Model.InitSequence(SysModel.m1x_0, SysModel.m2x_0, SysModel.T)

        x_Net_mdl_tst = Model(y_mdl_tst)

        MSE_test_linear_arr[j] = loss_fn(x_Net_mdl_tst, test_target[j, :, :]).item()

    # Average
    MSE_test_linear_avg = np.mean(MSE_test_linear_arr)
    MSE_test_dB_avg = 10 * np.log10(MSE_test_linear_avg)

    # Print MSE Cross Validation
    print("KalmanNet - MSE Test:", MSE_test_dB_avg, "[dB]")

    return [MSE_test_linear_arr, MSE_test_linear_avg, MSE_test_dB_avg]