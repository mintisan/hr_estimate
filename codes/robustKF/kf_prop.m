%************************************************************
% File:         kf_prop.m
% Date:         January 8, 2008
% Author:       Jo-Anne Ting
% Description:
%  Implements the propagation step of the standard Kalman filter
%
% Inputs:
%
%  x_prev  : value of the previous state
%  P       : previous state's covariance
%  A       : state transition matrix (assumed to be constant)
%  Q       : state noise (assumed to be constant)
%
% Outputs:
%
%  x       : propagation of the current state
%  P       : propagation of the current state's covariance
%
%************************************************************
function [x, P] = kf_prop(x_prev, P, A, Q)

% Propagation of state
x = A * x_prev;
  
% Propagation of covariance 
P = A * P * A' + Q;
