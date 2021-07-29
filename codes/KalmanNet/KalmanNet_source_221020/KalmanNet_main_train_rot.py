from KalmanNet_data import F_design, H_design, T, m1x_0_design, m2x_0_design, DataGen

from KalmanFilter_test import KFTest
from KalmanNet_build import NNBuild
from KalmanNet_train import NNTrain
from KalmanNet_test import NNTest
from KalmanNet_sysmdl import SystemModel
from KalmanNet_plt import KNPlot_test, NNPlot_train

from datetime import datetime

import numpy as np
import torch

####################
### Design Model ###
####################
q_design = 0.1
r_design = 0.1
SysModel_design = SystemModel(F_design, q_design, H_design, r_design, T)
SysModel_design.InitSequence(m1x_0_design, m2x_0_design)

##################
### Data Model ###
##################
q_data = q_design
r_data = r_design
F_data = F_design
H_data = H_design

SysModel_data = SystemModel(F_data, q_data, H_data, r_data, T)
SysModel_data.InitSequence(m1x_0_design, m2x_0_design)


print("Kalman Net - Train on Rotation")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

# Generate Data
#SysModel_data.GenerateBatch(N_T, 0)
#test_input = SysModel_data.Input
#test_target = SysModel_data.Target

print("Start Gen Data")
[train_input, train_target, cv_input, cv_target, test_input, test_target] = DataGen(SysModel_data)

# Test Kalman Filter with data parameters
[MSE_KF_data_linear_arr, MSE_KF_data_linear_avg, MSE_KF_data_dB_avg] = KFTest(SysModel_data, test_input, test_target)

# Test Kalman Filter with design parameters
[MSE_KF_design_linear_arr, MSE_KF_design_linear_avg, MSE_KF_design_dB_avg] = KFTest(SysModel_design, test_input, test_target)

# Build Neural Network
Model = NNBuild(SysModel_design)

# Train Neural Network
[MSE_cv_linear_epoch, MSE_cv_dB_epoch, MSE_train_linear_epoch, MSE_train_dB_epoch] = \
    NNTrain(SysModel_design, Model, cv_input, cv_target, train_input, train_target)

# Test Neural Network
[MSE_KNet_linear_arr, MSE_KNet_linear_avg, MSE_KNet_dB_avg] = NNTest(SysModel_design, test_input, test_target)

# Plot Training
NNPlot_train(MSE_KF_data_linear_arr, MSE_KF_data_dB_avg,
             MSE_KNet_linear_arr, MSE_KNet_dB_avg,
             MSE_cv_dB_epoch, MSE_train_dB_epoch)


KNPlot_test(MSE_KF_design_linear_arr, MSE_KF_data_linear_arr, MSE_KNet_linear_arr)