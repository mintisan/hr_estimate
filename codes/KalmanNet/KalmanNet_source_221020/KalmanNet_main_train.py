from KalmanNet_sysmdl import SystemModel

from KalmanNet_data import DataGen
from KalmanNet_data import F_design, H_design, T, m1x_0_design, m2x_0_design
#from KalmanNet_data import F_design, H_design, T, F_data, H_data, m1x_0_design, m2x_0_design

from KalmanFilter_test import KFTest
from KalmanNet_build import NNBuild
from KalmanNet_train import NNTrain
from KalmanNet_test import NNTest
from KalmanNet_plt import NNPlot_train

from datetime import datetime

print("Kalman Net Start")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

####################
### Design Model ###
####################
q_design = 1
r_design = 1
SysModel_design = SystemModel(F_design, q_design, H_design, r_design, T)
SysModel_design.InitSequence(m1x_0_design, m2x_0_design)

#####################
### Generate Data ###
#####################
print("Start Gen Data")
[train_input, train_target, cv_input, cv_target, test_input, test_target] = DataGen(SysModel_design)

print("KF Design")
# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_design, test_input, test_target)

# Build Neural Network
Model = NNBuild(SysModel_design)

# Train Neural Network
[MSE_cv_linear_epoch, MSE_cv_dB_epoch, MSE_train_linear_epoch, MSE_train_dB_epoch] = NNTrain(SysModel_design, Model, cv_input, cv_target, train_input, train_target)

# Test Neural Network
[MSE_test_linear_arr, MSE_test_linear_avg, MSE_test_dB_avg] = NNTest(SysModel_design, test_input, test_target)

# Plot
NNPlot_train(MSE_KF_linear_arr, MSE_KF_dB_avg,
             MSE_test_linear_arr, MSE_test_dB_avg,
             MSE_cv_dB_epoch, MSE_train_dB_epoch)
