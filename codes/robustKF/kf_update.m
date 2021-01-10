%************************************************************
% File:         kf_update.m
% Date:         January 8, 2008
% Author:       Jo-Anne Ting
% Description:
%  Implements the update step of the standard Kalman filter.
%  Calculates the posterior mean and covariance of the state.
%
% Inputs:
%
%  x_prev       : value of the previous state
%  z            : current observed data sample 
%  P            : propagated state covariance 
%  C            : observation matrix (assumed to be constant)
%  R            : observation noise (assumed to be constant)
% 
% Outputs:
%
%  x            : value of the current state
%  S            : posterior covariance of the prediction
%  P            : posterior covariance of the current state
%
%************************************************************
function  [x, S, P] = kf_update(x_prev, z, P, C, R)

dim = size(P,1);

S = C * P * C' + R;
K = P * C * inv(S);
x = x_prev + K * (z - C * x_prev);
P = (eye(dim) - K * C) * P;
