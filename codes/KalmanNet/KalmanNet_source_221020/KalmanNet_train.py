import torch
import torch.nn as nn
import numpy as np
import random

from KalmanNet_data import T, N_E, N_CV

# Number of Training Epochs
N_Epochs = 400

# Number of Samples in Batch
N_B = N_CV

# Learning Rate
learning_rate = 8 * 1e-4

# L2 Weight Regularization - Weight Decay
wd = 1e-5


def NNTrain(SysModel, Model, cv_input, cv_target, train_input, train_target):

    # MSE LOSS Function
    loss_fn = nn.MSELoss(reduction='mean')

    # Use the optim package to define an Optimizer that will update the weights of
    # the model for us. Here we will use Adam; the optim package contains many other
    # optimization algoriths. The first argument to the Adam constructor tells the
    # optimizer which Tensors it should update.
    optimizer = torch.optim.Adam(Model.parameters(), lr=learning_rate, weight_decay=wd)

    MSE_cv_linear_batch = np.empty([N_CV])
    MSE_cv_linear_epoch = np.empty([N_Epochs])
    MSE_cv_dB_epoch = np.empty([N_Epochs])

    MSE_train_linear_batch = np.empty([N_B])
    MSE_train_linear_epoch = np.empty([N_Epochs])
    MSE_train_dB_epoch = np.empty([N_Epochs])


    ##############
    ### Epochs ###
    ##############

    MSE_cv_dB_opt = 1000
    MSE_cv_idx_opt = 0

    for ti in range(0, N_Epochs):

        #################################
        ### Validation Sequence Batch ###
        #################################

        # Cross Validation Mode
        Model.eval()

        for j in range(0, N_CV):
            y_cv = cv_input[j, :, :]
            Model.InitSequence(SysModel.m1x_0, SysModel.m2x_0, T)
            x_Net_cv = Model(y_cv)

            # Compute Training Loss
            MSE_cv_linear_batch[j] = loss_fn(x_Net_cv, cv_target[j, :, :]).item()

        # Average
        MSE_cv_linear_epoch[ti] = np.mean(MSE_cv_linear_batch)
        MSE_cv_dB_epoch[ti] = 10 * np.log10(MSE_cv_linear_epoch[ti])

        if(MSE_cv_dB_epoch[ti] < MSE_cv_dB_opt):

            MSE_cv_dB_opt = MSE_cv_dB_epoch[ti]
            MSE_cv_idx_opt = ti
            torch.save(Model, 'best-model.pt')

        ###############################
        ### Training Sequence Batch ###
        ###############################

        # Training Mode
        Model.train()

        # Init Hidden State
        Model.init_hidden()

        Batch_Optimizing_LOSS_sum = 0

        for j in range(0, N_B):
            n_e = random.randint(0, N_E - 1)

            y_training = train_input[n_e, :, :]
            Model.InitSequence(SysModel.m1x_0, SysModel.m2x_0, T)
            x_Net_training = Model(y_training)

            # Compute Training Loss
            LOSS = loss_fn(x_Net_training, train_target[n_e, :, :])
            MSE_train_linear_batch[j] = LOSS.item()

            Batch_Optimizing_LOSS_sum = Batch_Optimizing_LOSS_sum + LOSS

        # Average
        MSE_train_linear_epoch[ti] = np.mean(MSE_train_linear_batch)
        MSE_train_dB_epoch[ti] = 10 * np.log10(MSE_train_linear_epoch[ti])

        ##################
        ### Optimizing ###
        ##################

        # Before the backward pass, use the optimizer object to zero all of the
        # gradients for the variables it will update (which are the learnable
        # weights of the model). This is because by default, gradients are
        # accumulated in buffers( i.e, not overwritten) whenever .backward()
        # is called. Checkout docs of torch.autograd.backward for more details.
        optimizer.zero_grad()

        # Backward pass: compute gradient of the loss with respect to model
        # parameters
        Batch_Optimizing_LOSS_mean = Batch_Optimizing_LOSS_sum / N_B
        Batch_Optimizing_LOSS_mean.backward()

        # Calling the step function on an Optimizer makes an update to its
        # parameters
        optimizer.step()

        ########################
        ### Training Summary ###
        ########################
        print(ti, "MSE Training :", MSE_train_dB_epoch[ti], "[dB]", "MSE Validation :", MSE_cv_dB_epoch[ti], "[dB]")

        if (ti > 1):
            d_train = MSE_train_dB_epoch[ti] - MSE_train_dB_epoch[ti - 1]
            d_cv    = MSE_cv_dB_epoch[ti] - MSE_cv_dB_epoch[ti - 1]
            print("diff MSE Training :", d_train, "[dB]", "diff MSE Validation :", d_cv, "[dB]")

        print("Optimal idx:", MSE_cv_idx_opt, "Optimal :", MSE_cv_dB_opt, "[dB]")

    return [MSE_cv_linear_epoch, MSE_cv_dB_epoch, MSE_train_linear_epoch, MSE_train_dB_epoch]