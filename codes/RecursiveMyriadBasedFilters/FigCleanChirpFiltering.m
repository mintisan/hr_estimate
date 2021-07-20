clear all;
close all;

%% Routine for reproduce the bandpass filtering operations on a clean chirp signal
%   Fig. 5. Bandpass Filtering Operations
%   Reference: 
%
%   [1] Ramirez, J., & Paredes, J. (2016). Recursive Weighted Myriad Based
%   Filters and their Optimizations. IEEE Transactions on Signal
%   Processing, 64(15), 4027-4039.
%
%   Author:
%   Juan Marcos Ramirez, M.S.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   September, 2016
%% Add Paths

addpath('RecursiveMyriadFilters/');
addpath('OthersNonlinearFilters/');

%% Training Stage
% Adaptive algorithms for obtaining the filter coefficients
% Note: FIR filter delay = #coeff/2

% FIR Filter Parameters
FILTER_ORDER        = 96;
firFilterDesign     = fir1(FILTER_ORDER - 1,[.075 .15]);

% Number of coefficients for recursive filters
N_FEEDFORWARD_WEIGHTS   = 64;
N_FEEDBACK_WEIGHTS      = 32;

% Generating the training and desired signals
tempTrainingSignal  = sign(randn(1,5047));
trainingSignal      = tempTrainingSignal(1:5000);
tempDesiredSignal   = filter(firFilterDesign,1,tempTrainingSignal);
desiredSignal       = tempDesiredSignal(48:end);

% Parameters for adaptive algorithms of the recursive myriad based filters
Ko   = 1;       % Initial Linearity Parameter
mu   = 0.001;  % Initial Step Size


% Adaptive Weighted Myriad Filter (Kalluri and Arce):
M = FILTER_ORDER;                % Number of filter coefficients (non-causal)
wmyWeights = (1/(M))*ones(1,M);  % Initial weights
[wmyWeights, wmyK] = adaptiveWMy(trainingSignal, ...
    desiredSignal, wmyWeights, Ko, mu, 'LMS');

M1 = N_FEEDFORWARD_WEIGHTS;    
M2 = N_FEEDBACK_WEIGHTS;

% Adaptive Weighted Median Filter
gwm = (1/(M1+M2))*ones(1,M1);
hwm = (1/(M1+M2))*ones(1,M2);
[gwm,hwm,ewm] = adaptiveRWM(trainingSignal,desiredSignal,...
                gwm,hwm,0.2);

% Adaptive Recursive Weighted Myriad (RWMy) Filter
rwmyH = ones(1,M2)/(M1 + M2);
rwmyG = ones(1,M1)/(M1 + M2);
rwmyK1 = 1;  % Non-recursive linearity parameter
rwmyK2 = 1;  % Recursive linearity parameter
[rwmyG, rwmyH] = adaptiveRWMy(trainingSignal,desiredSignal,rwmyG,rwmyH,...
    rwmyK1, rwmyK2, mu);
            
% Adaptive Recursive Hybrid Myriad (RHMy) Filter
rhmyH = ones(1,M2)/(M1 + M2);
rhmyG = ones(1,M1)/(M1 + M2);
rhmyK1 = 1;    % Non-recursive linearity parameter
[rhmyG, rhmyH] = adaptiveRHMy(trainingSignal,desiredSignal, rhmyG, rhmyH,...
    rhmyK1,mu);

% Adaptive Scaled Recursive Weighted Myriad (SRWMy) Filter
srwmyH = ones(1,M2)/(M1 + M2);
srwmyG = ones(1,M1)/(M1 + M2);
srwmyK1 = 1;  % Non-recursive linearity parameter
srwmyK2 = 1;  % Recursive linearity parameter
[srwmyG, srwmyH, srwmyK1, srwmyK2] = adaptiveSRWMy(trainingSignal, ...
    desiredSignal, srwmyG, srwmyH, srwmyK1, srwmyK2, mu);

% Adaptive Scaled Recursive Hybrid Myriad (SRHMy) Filter
srhmyH = ones(1,M2)/(M1 + M2);
srhmyG = ones(1,M1)/(M1 + M2);
srhmyK1 = 1;    % Non-recursive linearity parameter
[srhmyG, srhmyH, srhmyK1] = adaptiveSRHMy(trainingSignal,desiredSignal, srhmyG, srhmyH,...
    srhmyK1,mu);
 
%% Test stage

% Building the clean chirp signal
fs          = 2000;
n           = 0:1/fs:1.0 - 1/fs;
cleanChirp  = chirp(n,0,1,400);

% 96-tap FIR filtering on the clean chirp signal
tempFirOutput = filter(firFilterDesign, 1, cleanChirp);
firOutput =[tempFirOutput(48:end) tempFirOutput(1:47)];

% IIR filtering on the chirp signal
f = [0 0.075  0.075 0.125 0.150 1];
m = [0 0 1 1 0 0];
warning('off','all');
[b,a] = yulewalk(47,f,m);
tempIirOutput = filter(b,a,cleanChirp);
iirOutput=[tempIirOutput(46:end) tempIirOutput(1:45)];

% Weighted Myriad Filter Output
tic;
wmyOutput = wmyFilter(cleanChirp, wmyWeights, wmyK);
t(1) = toc;

% Dual Weighted Iterative Truncate Mean Filter (DWITM filter)
tic;
delayedChirp        = [zeros(1,FILTER_ORDER - 1) cleanChirp];
tempDwitmOutput     = zeros(1, length(cleanChirp));
for i = 1:length(cleanChirp)
    dwitmSamples        = delayedChirp(i:i + FILTER_ORDER - 1);
    dwitmWeights        = firFilterDesign;
    tempDwitmOutput(i)  = dwitm(dwitmSamples, dwitmWeights);
end
dwitmOutput      = [tempDwitmOutput(round(FILTER_ORDER/ 2):end)...
                    tempDwitmOutput(1:round(FILTER_ORDER / 2) - 1)];
t(2) = toc;

% Recursive Weighted Median Filter
tic;
rwmOutput = rwmFilter(cleanChirp,gwm,hwm);
t(3) = toc;

% Recursive Weighted Myriad Filter
tic;
rwmyOutput = rwmyFilter(cleanChirp, rwmyG, rwmyH, rwmyK1, rwmyK2);
t(4) = toc;

% Scaled Recursive Weighted Myriad Filter
tic;
srwmyOutput = srwmyFilter(cleanChirp, srwmyG, srwmyH, srwmyK1, srwmyK2);
t(5) = toc;

% Recursive Hybrid Myriad Filter
tic;
rhmyOutput = rhmyFilter(cleanChirp, rhmyG, rhmyH, rhmyK1);
t(6)=toc;


% Scaled Recursive Hybrid Myriad Filter
tic;
srhmyOutput = srhmyFilter(cleanChirp, srhmyG, srhmyH, srhmyK1);
t(7)=toc;

%% Display Results
% Mean absolute error (MAE) between the desired signal and the filtered
% signals

MAE(1) = mean(abs(firOutput - iirOutput));
MAE(2) = mean(abs(firOutput' - wmyOutput));
MAE(3) = mean(abs(firOutput - dwitmOutput));
MAE(4) = mean(abs(firOutput - rwmOutput));
MAE(5) = mean(abs(firOutput - rwmyOutput));
MAE(6) = mean(abs(firOutput' - srwmyOutput));
MAE(7) = mean(abs(firOutput - rhmyOutput));
MAE(8) = mean(abs(firOutput' - srhmyOutput));

% This data assignment is for saving signals in a specific format
cleanChirpSignals(:,1)     = n;
cleanChirpSignals(:,2)     = cleanChirp;
cleanChirpSignals(:,3)     = firOutput';
cleanChirpSignals(:,4)     = iirOutput';
cleanChirpSignals(:,5)     = wmyOutput;
cleanChirpSignals(:,6)     = dwitmOutput';
cleanChirpSignals(:,7)     = rwmOutput';
cleanChirpSignals(:,8)     = rwmyOutput';
cleanChirpSignals(:,9)     = srwmyOutput;
cleanChirpSignals(:,10)    = rhmyOutput';
cleanChirpSignals(:,11)    = srhmyOutput;

% Save signals in a .dat file
% save FigFilteringCleanChirp.dat cleanChirpSignals -ascii

figure;
subplot(5,2,1);
plot(n, cleanChirp), axis([0 0.55 -1.2 1.2]);
grid on;
title('Input Signal');

subplot(5,2,2);
plot(n, firOutput), axis([0 0.55 -1.2 1.2]);
grid on;
title('FIR output (desired signal)');

subplot(5,2,3);
plot(n, iirOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['IIR output. MAE = ' num2str(MAE(1),4)]);

subplot(5,2,4);
plot(n, wmyOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['Weighted Myriad. MAE = ' num2str(MAE(2),4)]);

subplot(5,2,5);
plot(n, dwitmOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['DWITM. MAE = ' num2str(MAE(3),4)]);

subplot(5,2,6);
plot(n, rwmOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['Recursive Weighted Median. MAE = ' num2str(MAE(4),4)]);

subplot(5,2,7);
plot(n, rwmyOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['Recursive Weighted Myriad. MAE = ' num2str(MAE(5),4)]);

subplot(5,2,8);
plot(n, srwmyOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['Scaled Recursive Weighted Myriad. MAE = ' num2str(MAE(6),4)]);
 
subplot(5,2,9);
plot(n, rhmyOutput), axis([0 0.55 -1.2 1.2]);
grid on;
title(['Recursive Hybrid Myriad. MAE = ' num2str(MAE(7),4)]);
 
subplot(5,2,10);
plot(n, srhmyOutput), hold on;
axis([0 0.55 -1.2 1.2]);
grid on;
title(['Scaled Recursive Hybrid Myriad. MAE = ' num2str(MAE(8),4)]);