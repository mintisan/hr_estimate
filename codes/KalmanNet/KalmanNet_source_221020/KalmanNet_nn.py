"""# **Class: KalmanNet**"""

import torch
import torch.nn as nn
import torch.nn.functional as func

class KalmanNetNN(torch.nn.Module):

    ###################
    ### Constructor ###
    ###################
    def __init__(self):
        super().__init__()

    ######################################
    ### Initialize Kalman Gain Network ###
    ######################################
    def InitKGainNet(self, H1, H2):
        # Input Dimensions
        D_in = self.m + self.n  # x(t-1), y(t)

        # Output Dimensions
        D_out = self.m * self.n;  # Kalman Gain

        ###################
        ### Input Layer ###
        ###################
        # Linear Layer
        self.KG_l1 = torch.nn.Linear(D_in, H1, bias=True)

        # ReLU (Rectified Linear Unit) Activation Function
        self.KG_relu1 = torch.nn.ReLU()

        ###########
        ### GRU ###
        ###########
        # Input Dimension
        self.input_dim = H1
        # Hidden Dimension
        self.hidden_dim = ((self.n * self.n) + (self.m * self.m)) * 10 * 1
        # Number of Layers
        self.n_layers = 1
        # Batch Size
        self.batch_size = 1
        # Input Sequence Length
        self.seq_len_input = 1
        # Hidden Sequence Length
        self.seq_len_hidden = self.n_layers

        # batch_first = False
        # dropout = 0.1 ;

        # Initialize a Tensor for GRU Input
        # self.GRU_in = torch.empty(self.seq_len_input, self.batch_size, self.input_dim)

        # Initialize a Tensor for Hidden State
        self.hn = torch.randn(self.seq_len_hidden, self.batch_size, self.hidden_dim)

        # Iniatialize GRU Layer
        self.rnn_GRU = nn.GRU(self.input_dim, self.hidden_dim, self.n_layers)

        ####################
        ### Hidden Layer ###
        ####################
        self.KG_l2 = torch.nn.Linear(self.hidden_dim, H2, bias=True)

        # ReLU (Rectified Linear Unit) Activation Function
        self.KG_relu2 = torch.nn.ReLU()

        ####################
        ### Output Layer ###
        ####################
        self.KG_l3 = torch.nn.Linear(H2, D_out, bias=True)

    ##################################
    ### Initialize System Dynamics ###
    ##################################
    def InitSystemDynamics(self, F, H):
        # Set State Evolution Matrix
        self.F = F
        self.F_T = torch.transpose(F, 0, 1);
        self.m = self.F.size()[0];

        # Set Observation Matrix
        self.H = H
        self.H_T = torch.transpose(H, 0, 1);
        self.n = self.H.size()[0];

    ###########################
    ### Initialize Sequence ###
    ###########################
    def InitSequence(self, M1_0, M2_0, T):
        self.m1x_posterior = M1_0

        self.T = T
        self.x_out = torch.empty(self.m, T)

        self.state_process_posterior_0 = M1_0

    ######################
    ### Compute Priors ###
    ######################
    def step_prior(self):
        # Predict the 1-st moment of x
        self.m1x_prior = torch.matmul(self.F, self.m1x_posterior);

        # Predict the 1-st moment of y
        self.m1y = torch.matmul(self.H, self.m1x_prior);

        self.state_process_prior_0 = torch.matmul(self.F, self.state_process_posterior_0)
        self.obs_process_0 = torch.matmul(self.H, self.state_process_prior_0)

    ##############################
    ### Kalman Gain Estimation ###
    ##############################
    def step_KGain_est(self, y):
        # Reshape and Normalize m1x Posterior
        m1x_post_0 = self.m1x_posterior - self.state_process_posterior_0
        #m1x_reshape = torch.squeeze(self.m1x_posterior)
        m1x_reshape = torch.squeeze(m1x_post_0)
        m1x_norm = func.normalize(m1x_reshape, p=2, dim=0, eps=1e-12, out=None);

        # Normalize y
        my_0 = y - torch.squeeze(self.obs_process_0)
        y_norm = func.normalize(my_0, p=2, dim=0, eps=1e-12, out=None);
        #y_norm = func.normalize(y, p=2, dim=0, eps=1e-12, out=None);

        # KGain Net Input
        KGainNet_in = torch.cat([m1x_norm, y_norm], dim=0)

        # Kalman Gain Network Step
        KG = self.KGain_step(KGainNet_in)

        # Reshape Kalman Gain to a Matrix
        self.KGain = torch.reshape(KG, (self.m, self.n))

    #######################
    ### Kalman Net Step ###
    #######################
    def KNet_step(self, y):
        # Compute Priors
        self.step_prior()

        # Compute Kalman Gain
        self.step_KGain_est(y)

        # Innovation
        y_obs = torch.unsqueeze(y, 1)
        dy = y_obs - self.m1y;

        # Compute the 1-st posterior moment
        INOV = torch.matmul(self.KGain, dy)
        self.m1x_posterior = self.m1x_prior + INOV;

        self.state_process_posterior_0 = self.state_process_prior_0

        # return
        return torch.squeeze(self.m1x_posterior)

    ########################
    ### Kalman Gain Step ###
    ########################
    def KGain_step(self, KGainNet_in):

        ###################
        ### Input Layer ###
        ###################
        L1_out = self.KG_l1(KGainNet_in);
        La1_out = self.KG_relu1(L1_out);

        ###########
        ### GRU ###
        ###########
        GRU_in = torch.empty(self.seq_len_input, self.batch_size, self.input_dim)
        GRU_in[0, 0, :] = La1_out
        GRU_out, self.hn = self.rnn_GRU(GRU_in, self.hn)
        GRU_out_reshape = torch.reshape(GRU_out, (1, self.hidden_dim))

        ####################
        ### Hidden Layer ###
        ####################
        L2_out = self.KG_l2(GRU_out_reshape);
        La2_out = self.KG_relu2(L2_out);

        ####################
        ### Output Layer ###
        ####################
        L3_out = self.KG_l3(La2_out);
        return L3_out;

    ###############
    ### Forward ###
    ###############
    def forward(self, y):
        for t in range(0, self.T):
            self.x_out[:, t] = self.KNet_step(y[:, t])

        return self.x_out

    #########################
    ### Init Hidden State ###
    #########################
    def init_hidden(self):
        weight = next(self.parameters()).data
        hidden = weight.new(self.n_layers, self.batch_size, self.hidden_dim).zero_()
        self.hn = hidden.data