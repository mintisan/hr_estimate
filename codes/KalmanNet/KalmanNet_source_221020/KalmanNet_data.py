import torch
import h5py
import math

from KalmanNet_sysmdl import SystemModel

#######################
### Size of DataSet ###
#######################

# Number of Training Examples
N_E = 20000

# Number of Cross Validation Examples
N_CV = 100

# Number of Test Examples
N_T = 1

##################
### Design #10 ###
##################
F10 = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

H10 = torch.tensor([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

#############
### 2 x 2 ###
#############
m = 2
n = 2
F_design = F10[0:m, 0:m]
#H_design = torch.eye(2)
H_design = H10[0:n, 10-m:10]
m1x_0_design = torch.tensor([[5.0], [-5.0]])
m2x_0_design = 0 * 0 * torch.eye(m)

############
### Data ###
############
theta_deg = 1 * 10
theta_rad = (2 / 360) * (math.pi) * theta_deg
F_rot = torch.tensor([[math.cos(theta_rad), -math.sin(theta_rad)], [math.sin(theta_rad), math.cos(theta_rad)]])
F_data = F_rot * F_design

#############
### 5 x 5 ###
#############
#m = 5
#n = 5
#F_design = F10[0:m, 0:m]
#H_design = H10[0:n, 10-m:10]
#m1x_0_design = torch.zeros(m, 1)
#m1x_0_design = torch.tensor([[1.0], [-1.0], [2.0], [-2.0], [0.0]])
#m2x_0_design = 0 * 0 * torch.eye(m)

# Length of Time Series Sequence
T = 10



def DataGen(SysModel_data):

    ##################################
    ### Generate Training Sequence ###
    ##################################
    SysModel_data.GenerateBatch(N_E, 0.5)
    training_input = SysModel_data.Input
    training_target = SysModel_data.Target

    ####################################
    ### Generate Validation Sequence ###
    ####################################
    SysModel_data.GenerateBatch(N_CV, 0.5)
    cv_input = SysModel_data.Input
    cv_target = SysModel_data.Target

    ##############################
    ### Generate Test Sequence ###
    ##############################
    SysModel_data.GenerateBatch(N_T, 0.5)
    test_input = SysModel_data.Input
    test_target = SysModel_data.Target

    # return
    return [training_input, training_target, cv_input, cv_target, test_input, test_target]