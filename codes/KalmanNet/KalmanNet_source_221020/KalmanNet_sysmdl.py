import torch
import numpy as np

class SystemModel:

    def __init__(self, F, q, H, r, T):

        ####################
        ### Motion Model ###
        ####################
        self.F = F
        self.m = self.F.size()[0]

        self.q = q
        self.Q = q * q * torch.eye(self.m)

        #########################
        ### Observation Model ###
        #########################
        self.H = H
        self.n = self.H.size()[0]

        self.r = r
        self.R = r * r * torch.eye(self.n)

        ################
        ### Sequence ###
        ################
        # Assign T
        self.T = T

        # Pre allocate an array for current state
        self.x = torch.empty(size=[self.m, self.T])

        # Pre allocate an array for current observation
        self.y = torch.empty(size=[self.n, self.T])

    #####################
    ### Init Sequence ###
    #####################
    def InitSequence(self, m1x_0, m2x_0):

        self.m1x_0 = m1x_0
        self.m2x_0 = m2x_0


    #########################
    ### Update Covariance ###
    #########################
    def UpdateCovariance_Gain(self, q, r):

        self.q = q
        self.Q = q * q * torch.eye(self.m)

        self.r = r
        self.R = r * r * torch.eye(self.n)

    def UpdateCovariance_Matrix(self, Q, R):

        self.Q = Q

        self.R = R


    #########################
    ### Generate Sequence ###
    #########################
    def GenerateSequence(self, Q_gen, R_gen):

        # Set x0 to be x previous
        self.x_prev = self.m1x_0

        # Generate Sequence Iteratively
        for t in range(0, self.T):
            ########################
            #### State Evolution ###
            ########################
            xt = self.F.matmul(self.x_prev)

            # Process Noise
            mean = torch.zeros(self.m)
            eq = np.random.multivariate_normal(mean, Q_gen, 1)
            eq = torch.transpose(torch.tensor(eq), 0, 1)
            eq = eq.type(torch.float)

            # Additive Process Noise
            xt = xt.add(eq)

            ################
            ### Emission ###
            ################
            yt = self.H.matmul(xt)

            # Observation Noise
            mean = torch.zeros(self.n)
            er = np.random.multivariate_normal(mean, R_gen, 1)
            er = torch.transpose(torch.tensor(er), 0, 1)

            # Additive Observation Noise
            yt = yt.add(er)

            ########################
            ### Squeeze to Array ###
            ########################

            # Save Current State to Trajectory Array
            self.x[:, t] = torch.squeeze(xt)

            # Save Current Observation to Trajectory Array
            self.y[:, t] = torch.squeeze(yt)

            ################################
            ### Save Current to Previous ###
            ################################
            self.x_prev = xt

    ######################
    ### Generate Batch ###
    ######################
    def GenerateBatch(self, size, gain):

        # Allocate Empty Array for Input
        self.Input = torch.empty(size, self.n, self.T)

        # Allocate Empty Array for Target
        self.Target = torch.empty(size, self.m, self.T)

        ### Generate Examples
        for i in range(0, size):
            # Generate Sequence

            #[Q_gen, R_gen] = self.sampling(self.q, self.r, gain)

            self.GenerateSequence(self.Q, self.R)

            # Training sequence input
            self.Input[i, :, :] = self.y

            # Training sequence output
            self.Target[i, :, :] = self.x


    def sampling(self, q, r, gain):

        if (gain != 0):
            gain_q = 0.1
            #aq = gain * q * np.random.randn(self.m, self.m)
            aq = gain_q * q * torch.eye(self.m)
            #aq = gain_q * q * torch.tensor([[1.0, 1.0], [1.0, 1.0]])
        else:
            aq = 0

        Aq = q * torch.eye(self.m) + aq
        Q_gen = np.transpose(Aq) * Aq

        if (gain != 0):
            gain_r = 0.5
            #ar = gain * r * np.random.randn(self.n, self.n)
            ar = gain_r * r * torch.eye(self.n)
            #ar = gain_r * r * torch.tensor([[1.0, 1.0], [1.0, 1.0]])

        else:
            ar = 0

        Ar = r * torch.eye(self.n) + ar
        R_gen = np.transpose(Ar) * Ar

        return [Q_gen, R_gen]