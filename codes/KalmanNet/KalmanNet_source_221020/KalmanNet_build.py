from KalmanNet_nn import KalmanNetNN

def NNBuild(SysModel):

    Model = KalmanNetNN()

    Model.InitSystemDynamics(SysModel.F, SysModel.H)

    # Number of neurons in the 1st hidden layer
    H1_KNet = (SysModel.m + SysModel.n) * (10) * 8

    # Number of neurons in the 2nd hidden layer
    H2_KNet = (SysModel.m * SysModel.n) * 1 * (4);

    Model.InitKGainNet(H1_KNet, H2_KNet)

    return Model