% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine for obtaining the behavior of the adaptive algorithms as they
% progress. Figure 4.
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

disp('---------------------------------------');
disp('This routine could take several minutes');
disp('---------------------------------------');

% Routine for obtaining the behavior of the adaptive algorithms
addpath('RecursiveMyriadFilters/');

% Adaptive recursive myriad based filter parameters
Ki      = 1;
u       = 0.001;
M1      = 64;
M2      = 32;

% Training signal parameters
N       = 5000;

% Simulation parameters
trials  = 10;

% FIR filter parameters
FILTER_ORDER    = 96;
firFilterDesign = fir1(FILTER_ORDER - 1,[.075 .15]);

maeRWMyEnsembleAverage  = zeros(N - M1 + 1,1);
maeRHMyEnsembleAverage  = zeros(N - M1 + 1,1);
maeRMMyEnsembleAverage  = zeros(N - M1 + 1,1);

maeSRWMyEnsembleAverage  = zeros(N - M1 + 1,1);
maeSRHMyEnsembleAverage  = zeros(N - M1 + 1,1);

gSRHMyEnsembleAverage = zeros(N - M1 + 1,M1);
hSRHMyEnsembleAverage = zeros(N - M1 + 1,M2);
KSRHMyEnsembleAverage = zeros(1,N - M1 + 1);

gSRWMyEnsembleAverage = zeros(N - M1 + 1,M1);
hSRWMyEnsembleAverage = zeros(N - M1 + 1,M2);
K1SRWMyEnsembleAverage = zeros(1,N - M1 + 1);
K2SRWMyEnsembleAverage = zeros(1,N - M1 + 1);

disp(['Number of realizations: ' num2str(trials)])
for ii = 1:trials
    tic;
    disp(['Iteration:   ' num2str(ii)]);
    % Training signal settings
    tempTrainingSignal  = sign(randn(1,N + 48));
    trainingSignal      = tempTrainingSignal(1:N);
    tempDesiredSignal   = filter(firFilterDesign,1,tempTrainingSignal);
    desiredSignal       = tempDesiredSignal(48:end);
    
    h = ones(1,M2)/(M1 + M2);
    g = ones(1,M1)/(M1 + M2);    
    
    [~, ~, eRWMy] = adaptiveRWMy(trainingSignal,desiredSignal, g, h, Ki, Ki, u);
    [~, ~, eRHMy] = adaptiveRHMy(trainingSignal,desiredSignal, g, h, Ki, u);
    
    maeRWMyEnsembleAverage = maeRWMyEnsembleAverage + abs(eRWMy)/trials;
    maeRHMyEnsembleAverage = maeRHMyEnsembleAverage + abs(eRHMy)/trials;
    
    [~, ~, ~, ~, eSRWMy, gSRWMy, hSRWMy, K1SRWMy, K2SRWMy] = adaptiveSRWMy(trainingSignal,desiredSignal, g, h, Ki, Ki, u);
    [~, ~, ~, eSRHMy, gSRHMy, hSRHMy, KSRHMy] = adaptiveSRHMy(trainingSignal,desiredSignal, g, h, Ki, u);
    
    maeSRWMyEnsembleAverage = maeSRWMyEnsembleAverage + abs(eSRWMy)/trials;
    maeSRHMyEnsembleAverage = maeSRHMyEnsembleAverage + abs(eSRHMy)/trials;
    
    gSRHMyEnsembleAverage = gSRHMyEnsembleAverage + gSRHMy/trials;
    hSRHMyEnsembleAverage = hSRHMyEnsembleAverage + hSRHMy/trials;
    KSRHMyEnsembleAverage = KSRHMyEnsembleAverage + KSRHMy/trials;
    
    gSRWMyEnsembleAverage = gSRWMyEnsembleAverage + gSRWMy/trials;
    hSRWMyEnsembleAverage = hSRWMyEnsembleAverage + hSRWMy/trials;
    K1SRWMyEnsembleAverage = K1SRWMyEnsembleAverage + K1SRWMy/trials;
    K2SRWMyEnsembleAverage = K2SRWMyEnsembleAverage + K2SRWMy/trials;
    toc;
end

maeRWMyTrial    = abs(eRWMy);
maeRHMyTrial    = abs(eRHMy);
maeSRWMyTrial    = abs(eSRWMy);
maeSRHMyTrial    = abs(eSRHMy);

%% Display Results

n = 1:1:N-M1+1;

n_rand1 = randi([1 M1],1);
n_rand2 = randi([1 M2],1);

subplot(321)
plot(n, gSRWMy(:,n_rand1)), hold on;
plot(n, gSRWMyEnsembleAverage(:,n_rand1), 'k');
plot(n, hSRWMy(:,n_rand2), 'r');
plot(n, hSRWMyEnsembleAverage(:,n_rand2), 'k--');
title('Scaled Recursive Weighted Myriad')
ylabel('Parameter value'), xlabel('Iteration (n)');
legend('g (trial)', 'g (average)', 'h (trial)', 'h (average)');
axis('tight');

subplot(322)
plot(n, K1SRWMy), hold on;
plot(n, K1SRWMyEnsembleAverage, 'k');
plot(n, K2SRWMy, 'r');
plot(n, K2SRWMyEnsembleAverage, 'k--');
title('Scaled Recursive Weighted Myriad')
ylabel('Linearity parameters'), xlabel('Iteration (n)');
legend('K1 (trial)', 'K1 (average)', 'K2 (trial)', 'K2 (average)');
axis('tight');

subplot(323)
plot(n, maeRWMyTrial), hold on;
plot(n, maeRWMyEnsembleAverage, 'k');
legend('Single Trial', 'Ensemble Average');
title('Recursive Weighted Myriad')
ylabel('MAE'), xlabel('Iteration (n)');
axis('tight');

subplot(324)
plot(n, maeRHMyTrial), hold on;
plot(n, maeRHMyEnsembleAverage, 'k');
legend('Single Trial', 'Ensemble Average');
title('Recursive Hybrid Myriad')
ylabel('MAE'), xlabel('Iteration (n)');
axis('tight');

subplot(325)
plot(n, maeSRWMyEnsembleAverage), hold on;
plot(n, maeRWMyEnsembleAverage, 'k');
legend('Scaled', 'Normalized');
title('Recursive Weighted Myriad')
ylabel('MAE'), xlabel('Iteration (n)');
axis('tight');

subplot(326)
plot(n, maeSRHMyEnsembleAverage), hold on;
plot(n, maeRHMyEnsembleAverage, 'k');
legend('Scaled', 'Normalized');
title('Recursive Hybrid Myriad')
ylabel('MAE'), xlabel('Iteration (n)');
axis('tight');