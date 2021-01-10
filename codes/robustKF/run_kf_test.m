%************************************************************
% File: 	run_kf_test.m
% Date: 	January 8, 2008
% Author: 	Jo-Anne Ting
% Description:
%   Compares performances of:
%
%  i)   the standard Kalman filter, 
%  ii)  thresholded Kalman filter (where the threshold is 
%       optimally determined given prior knowledge of the 
%       entire data set), 
%  iii) our proposed weighted outlier-robust Kalman filter 
%       (as described in Ting et al., 2007) 
%
%   on a real-time outlier detection problem, using
%   a data set consisting of LittleDog's orientation
%   information (in quarternion form). 
%  
%   This test assumes that the data is streaming, with a 
%   new data sample arriving at each time step. We 
%   run this test on one of the quarternion coefficients
%   (i.e., the training data has 1-dimensional inputs). 
%
%************************************************************
clear all;
close all;


% Load dog quarternion data
load Qoff.data

N      = size(Qoff,1);    % number of data samples
% Choose one of the 4 quarternion coefficients for this test
numDim = 1;               % number of input dimensions
Yn     = Qoff(:,4);


%------------------------------------------------------------
% 1. Set values of indicator variables for test settings
%------------------------------------------------------------

% Set indicator variable to 1 if want to run algorithm
% (use the plain Kalman filter as a standard comparison)

THRESHOLDED_KF_ON     = 0;   % thresholded Kalman filter
HACKED_WEIGHTED_KF_ON = 1;   % "hacked" weighted Kalman filter
WEIGHTED_ROBUST_KF_ON = 0;   % weighted robust Kalman filter


%------------------------------------------------------------
% 2. Perform initialization needed for all algorithms
%------------------------------------------------------------

% Set the covariance matrix to a small initial value
P = 0.01*ones(numDim);       % for the standard KF
if THRESHOLDED_KF_ON == 1 
  P_thres = P;
end
if HACKED_WEIGHTED_KF_ON == 1
  P_hKF = P;
end
if WEIGHTED_ROBUST_KF_ON == 1
  P_wrKF = P;
end

% Initialize the state estimate, x_hat, to some random value
x_hat(1,:) = rand;           % for the standard KF
if THRESHOLDED_KF_ON == 1 
  x_hat_thres(1,:)   = x_hat(1,1);
end
if HACKED_WEIGHTED_KF_ON == 1
  x_hat_hKF(1,:)   = x_hat(1,1);
end
if WEIGHTED_ROBUST_KF_ON == 1
  x_hat_wrKF(1,:)  = x_hat(1,1);
end

% Initialize sufficient statistics to be collected 
if THRESHOLDED_KF_ON == 1
  ss_thres.sum_zxT = 0;
  ss_thres.sum_xxT = 0;
  ss_thres.sum_xxold = 0;
  ss_thres.sum_xxoldT = 0;
  ss_thres.sum_N = 0;
  ss_thres.sum_zz = 0;
  ss_thres.sum_zx = 0;
  ss_thres.sum_ExTx = 0;
  ss_thres.sum_Exxold = 0;
end
if HACKED_WEIGHTED_KF_ON == 1
  ss_hKF.sum_wzxT = 0;       
  ss_hKF.sum_wxxT = 0;
  ss_hKF.sum_xxold = 0;
  ss_hKF.sum_xxoldT = 0;
  ss_hKF.sum_N = 0;
  ss_hKF.sum_wzz = 0;
  ss_hKF.sum_wzx = 0;
  ss_hKF.sum_ExTx = 0;
  ss_hKF.sum_Exxold = 0;
end
if WEIGHTED_ROBUST_KF_ON == 1
  ss_wrKF.sum_wzxT = 0;                
  ss_wrKF.sum_wxxT = 0;
  ss_wrKF.sum_xxold = 0;
  ss_wrKF.sum_xxoldT = 0;
  ss_wrKF.sum_N = 0;
  ss_wrKF.sum_wzz = 0;
  ss_wrKF.sum_wzx = 0;
  ss_wrKF.sum_ExTx = 0;
  ss_wrKF.sum_Exxold = 0;
end

% Initialize the system matrices

% For the standard Kalman filter
A0 = 1;     % state transition matrix
C0 = 1;     % observation matrix    
Q0 = 1e-4;  % state noise variance
R0 = 1e-4;  % observation noise variance

if THRESHOLDED_KF_ON == 1
  A_thres = 1;
  C_thres = 1;
  R_thres = 1e-4;
  Q_thres = 1e-4;
end
if HACKED_WEIGHTED_KF_ON == 1
  A_hKF = 1;     
  C_hKF = 1;
  Q_hKF = 1e-4;  
  R_hKF = 1e-4;
end
if WEIGHTED_ROBUST_KF_ON == 1
  A_wrKF = 1;      
  C_wrKF = 1;
  Q_wrKF = 1e-4;   
  R_wrKF = 1e-4; 
end


%------------------------------------------------------------
% 3. Run the algorithms 
%------------------------------------------------------------

% Assume one data sample arrives at each time step
for i=1:N

   % Run the standard Kalman filter 
   % (system dynamics used are constant, ie., A0, C0, Q0, R0)
   %---------------------------------------------------------
   [x_hat(i+1,:), P] = kf_prop(x_hat(i,:), P, A0, Q0);
   [x_hat(i+1,:), S, P] = ...
       kf_update(x_hat(i+1,:), Yn(i,:), P, C0, R0);

   % Run the "hacked" weighted Kalman filter 
   % (where systems dynamics, A, C, Q, R, are learnt) 
   %---------------------------------------------------------
   if HACKED_WEIGHTED_KF_ON == 1 
   [x_hat_hKF(i+1,:), hKF_weights(i), S_hKF, P_hKF, A_hKF, ...
    C_hKF, Q_hKF, R_hKF, ss_hKF] = ...
    hackedwKF_learn(x_hat_hKF(i,:), Yn(i,:), P_hKF, A_hKF, ...
                     C_hKF, Q_hKF, R_hKF, ss_hKF);
   end

   % Run the weighted robust Kalman filter
   %---------------------------------------------------------
   if WEIGHTED_ROBUST_KF_ON == 1
   [x_hat_wrKF(i+1,:), wrKF_weights(i), S_wrKF, P_wrKF, A_wrKF, ...
    C_wrKF, Q_wrKF, R_wrKF, ss_wrKF] =...
    wrKF_learn(x_hat_wrKF(i,:), Yn(i,:), P_wrKF, A_wrKF, ...
               C_wrKF, Q_wrKF, R_wrKF, ss_wrKF);
   end

   % Run the thresholded Kalman filter 
   % (threshold value is for the Mahalanobis distance and its
   % optimal value is experimentally found, given the entire
   % data set)
   %---------------------------------------------------------
   if THRESHOLDED_KF_ON == 1 
   % Use the learnt system dynamics of wrKF as initial values
   A_thres = A_wrKF; 
   C_thres = C_wrKF; 
   R_thres = R_wrKF; 
   Q_thres = Q_wrKF;
   
   % Optimal threshold value 
   % (found experimentally for each data set)
   gamma = 0.1;  

   [x_hat_thres(i+1,:), P_thres] = ...
        kf_prop(x_hat_thres(i,:), P_thres, A_thres, Q_thres);
   [x_hat_thres(i+1,:), S_thres, P_thres, A_thres, ...
    C_thres, Q_thres, R_thres, ss_thres] = ...
    kf_update_learn(x_hat_thres(i+1,:), Yn(i,:), P_thres, ...
                    A_thres, C_thres, Q_thres, R_thres, ...
                    gamma, ss_thres);
   end

   % Calculate predicted output
   %---------------------------------------------------------
   z_predicted(i+1,:) = C0 * x_hat(i+1,:);
   if THRESHOLDED_KF_ON == 1
   z_predicted_thres(i+1,:) = C_thres * x_hat_thres(i+1,:);
   end
   if HACKED_WEIGHTED_KF_ON == 1
   z_predicted_hKF(i+1,:) = C_hKF * x_hat_hKF(i+1,:);
   end
   if WEIGHTED_ROBUST_KF_ON == 1
   z_predicted_wrKF(i+1,:) = C_wrKF * x_hat_wrKF(i+1,:);
   end

   % Store the posterior covariance of state & output
   %---------------------------------------------------------
   x_cov(i+1,:) = P;
   z_cov(i+1,:) = S;
   if HACKED_WEIGHTED_KF_ON == 1
   x_cov_hKF(i+1,:) = P_hKF;
   z_cov_hKF(i+1,:) = S_hKF;
   end
   if WEIGHTED_ROBUST_KF_ON == 1
   x_cov_wrKF(i+1,:) = P_wrKF;
   z_cov_wrKF(i+1,:) = S_wrKF;
   end

end  % for i=1:N


%------------------------------------------------------------
% 3. Plot/Visualize results 
%------------------------------------------------------------

X = [1:N]';

% Plot observed data
%--------------------
figure(); 
hold on
plot(X, Yn, 'og');
legend('Observed output')
xlabel('Time step');
ylabel('Output data');
hold off;

% Plot predicted vs. actual observations for standard KF
% and weighted outlier-robust KF
%-----------------------------------------------------------
figure(); 
hold on
plot(X, Yn, 'g.');
plot(X, z_predicted(2:end),'c');
if WEIGHTED_ROBUST_KF_ON == 1
  plot(X, z_predicted_wrKF(2:end),'r');
  legend('Observed output', 'Standard KF', ...
         'Weighted Robust KF')
end
xlabel('Time step');
ylabel('Output data');
hold off;

% Plot predicted vs. actual observations for the thresholded 
% KF, "hacked" weighted KF and weighted outlier-robust KF
%-----------------------------------------------------------
figure(); 
hold on
plot(X, Yn, 'g.');
if THRESHOLDED_KF_ON == 1
  plot(X, z_predicted_thres(2:end),'k');
end
if HACKED_WEIGHTED_KF_ON == 1
  plot(X, z_predicted_hKF(2:end),'c');
end
if WEIGHTED_ROBUST_KF_ON == 1
  plot(X, z_predicted_wrKF(2:end),'r');
end
xlabel('Time step');
ylabel('Output data');
legend('Observed Output', 'Thresholded KF',...
       'Alternative KF', 'Weighted Robust KF')
hold off;
