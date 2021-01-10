%************************************************************
% File:         wrKF_learn.m
% Date:         January 8, 2008
% Author:       Jo-Anne Ting
% Description:
%  Implements the weighted outlier-robust Kalman filter (as
%  described in Ting et al., 2007) where the system
%  dynamics (A, C, Q, R) are learnt as well
%
% Inputs:
%
%  x_prev       : value of the previous state
%  z            : current observed data sample
%  P            : propagated state covariance
%  A            : state transition matrix
%  C            : observation matrix
%  Q            : state noise variance
%  R            : observation noise variance
%  ss           : sufficient statistics
%
% Outputs:
%
%  x            : value of the current state
%  weight       : posterior weight of the data sample
%  S            : posterior covariance of the prediction
%  P            : posterior covariance of the current state
%  A            : updated state transition matrix
%  C            : updated observation matrix
%  Q            : updated state noise variance
%  R            : updated observation noise variance
%  ss           : updated sufficient statistics
%
%************************************************************
function [x, weight, S, P, A, C, Q, R, ss] = ...
         wrKF_learn(x_prev, z, P, A, C, Q, R, ss) 

% A small constant to avoid computationally problems
SMALL = 1e-6;

% Initial priors for weight of the observed data sample
%-------------------------------------------------------
alpha = 1e-6;  
beta  = 1e-6; 

% Calculate posterior mean and covariance of state
%---------------------------------------------------
oldP = P;
r = z - C*x_prev;
omega = r'*r + trace(C'*C*oldP); 
weight = (alpha + 0.5)/(beta + (0.5/(R+SMALL))*omega); 

S = C*Q*C' + R/(weight+SMALL);
P = inv(weight*(C'/(R+SMALL))*C + 1/Q + SMALL);
x = P/Q*A*x_prev  + weight*P*C'/(R+SMALL)*z;

% Update sufficient statistics for A, C, Q and R
%--------------------------------------------------
ss.sum_wzxT = ss.sum_wzxT + weight*z*x';
ss.sum_wxxT = ss.sum_wxxT + weight*(x*x' + P);
ss.sum_xxold = ss.sum_xxold + x*x_prev;
ss.sum_xxoldT = ss.sum_xxoldT + (x_prev*x_prev' + oldP);

ss.sum_N = ss.sum_N + 1;
ss.sum_wzz = ss.sum_wzz + weight*z'*z;
ss.sum_wzx = ss.sum_wzx + weight*(z.*x);
ss.sum_ExTx = ss.sum_ExTx + (x'*x + P);
ss.sum_Exxold = ss.sum_Exxold + x.*x_prev;

% Calculate new system dynamics
%--------------------------------
A = ss.sum_xxold / ss.sum_xxoldT;
C = ss.sum_wzxT / ss.sum_wxxT;
R = ss.sum_wzz - 2*sum(ss.sum_wzx'*C) + diag(C*ss.sum_wxxT*C);
R = R / ss.sum_N;
Q = ss.sum_ExTx - 2*sum(ss.sum_Exxold'*A) + diag(A*ss.sum_xxoldT*A);
Q = Q / ss.sum_N;
