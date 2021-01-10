%************************************************************
% File:         kf_update_learn.m
% Date:         January 8, 2008
% Author:       Jo-Anne Ting
% Description:
%  Implements the update step of the standard Kalman filter,
%  except it finds if the Mahalanobis distance of the
%  current data sample exceeds a certain threshold value
%  or not. If this threshold value is exceeded, then 
%  the data sample is treated as an outlier and no update
%  step is performed.
%
%  Additionally, the system dynamics (A, C, Q, R) are learnt
%  using a maximum-likelihood framework. Refer to Ting et
%  al., 2007 for more details.
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
%  gamma        : threshold value (that is found experimentally)
%  ss           : sufficient statistics
%
% Outputs:
%
%  x            : value of the current state
%  S            : posterior covariance of the prediction
%  P            : posterior covariance of the current state
%  A            : updated state transition matrix
%  C            : updated observation matrix
%  Q            : updated state noise variance
%  R            : updated observation noise variance
%  ss           : updated sufficient statistics
%
%************************************************************
function [x, S, P, A, C, Q, R, ss] = ...
         kf_update_learn(x_prev, z, P, A, C, Q, R, gamma, ss)

% A small constant to avoid computationally problems
SMALL = 1e-6;

oldP = P;
r = z - C * x_prev;
S = C * P * C' + R;
if (r' * inv(S + SMALL) * r <= gamma)
  K = P * C' * inv(S + SMALL); 
  x = x_prev + K * r;
  P = P  - P * C' * inv(S + SMALL) * C * P;
else
  x = x_prev;
  P = P; 
end

% Update sufficient statistics only if sample was not
% an outlier

if (r' * inv(S + SMALL) * r <= gamma)
% Update sufficient statistics for A, C, Q and R
%--------------------------------------------------
ss.sum_zxT = ss.sum_zxT + z*x;
ss.sum_xxT = ss.sum_xxT + (x*x' + P);
ss.sum_xxold = ss.sum_xxold + x*x_prev;
ss.sum_xxoldT = ss.sum_xxoldT + (x_prev^2 + oldP);

ss.sum_N = ss.sum_N + 1;
ss.sum_zz = ss.sum_zz + z'*z;
ss.sum_zx = ss.sum_zx + z.*x;
ss.sum_ExTx = ss.sum_ExTx + (x'*x + P);
ss.sum_Exxold = ss.sum_Exxold + x.*x_prev;

% Calculate new system dynamics
%--------------------------------
A = ss.sum_xxold / ss.sum_xxoldT;
C = ss.sum_zxT / ss.sum_xxT;
R = ss.sum_zz - 2*sum(ss.sum_zx'*C) + diag(C*ss.sum_xxT*C);
R = R / ss.sum_N;
Q = ss.sum_ExTx - 2*sum(ss.sum_Exxold'*A) + diag(A*ss.sum_xxoldT*A);
Q = Q / ss.sum_N;
end

