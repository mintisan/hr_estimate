% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Montecarlo Simulation for computing the performace of the proposed
% recursive filters under various performance criteria, for different
% impulsiveness levels of the additive noise. Table I. 
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
%% Add Paths

addpath('RecursiveMyriadFilters/');
addpath('OthersNonlinearFilters/')

%% Starting the Montecarlo simulation

alpha = [0.75 1.00 1.50 2];
trials = 10;
u = 0.001;

% Noiseless chirp signal
fs = 2000;
n = 0:1/fs:0.5 - 1/fs;
s = chirp(n,0,1,400);

% FIR filtering on the clean chirp signal 96-tap
win1 = 96;
wfir1 = fir1(win1-1,[.075 .125]);
d1 = filter(wfir1,1,s);
ax1 = d1(1:47);
d1=[d1(47+1:end) ax1];

% IIR filter parameters
f = [0 0.075  0.075 0.125 0.125 1];
m = [0 0 1 1 0 0];
warning('off');
[a,b] = yulewalk(47,f,m);

% Parameters under test
dispersion = 0.1;

% Linearity parameters
Kwf =[0.35 0.35 0.35 0.35 0.35];
Khf =[0.50 0.50 0.50 0.50 0.50];

for ii = 1 : 5
    if (ii <= 4)
        disp(['Characteristic Exponent (alpha):' num2str(alpha(ii),2)]);
    else
        disp('Laplacian Noise')
    end
    
    disp('Starting the training stage');
    su = sign(randn(1,10047));
    so =  su(1:10000);
    if (ii <= 4)
        s1 = so + astable(1,10000,alpha(ii),0,dispersion,0);
    else
        variance = 2*dispersion;
        s1 = so + laprnd(1,10000,0,sqrt(variance));
    end
    

    % Building the desired signal for adaptive filter algorithms
    win = 96;
    wfir = fir1(win-1,[.075 .125]);
    d = filter(wfir,1,so);
    ax = d(1:47);                           %Revisar fvtool...
    d=[d(47+1:end) ax];

    % Adaptive recursive weighted median filter (RWM filter)
    tic;
    M = 64; N = 32;
    g = (1/(M+N))*ones(1,M);
    h = (1/(M+N))*ones(1,N);
    [gwm,hwm,ewm] = adaptiveRWM(s1,d,g,h,0.01);
    trainingTime(1) = toc;

    % Adaptive (non recursive) weighted myriad filter
    tic;
    Mnr = 96;
    w = (1/(Mnr))*ones(1,Mnr);
    [w1,e1] = adaptiveWMy(s1,d,w,1,u*10,1);
    trainingTime(2) = toc;

    % Adaptive recursive weighted myriad filter (RWMy filter)
    tic;
    [g2,h2] = adaptiveRWMy(s1,d,g,h,Kwf(ii),Kwf(ii),u);
    trainingTime(3) = toc;

    % Adaptive scaled recursive weighted myriad filter (SRWMy filter)
    tic;
    [go,ho,K1,K2] = adaptiveSRWMy(s1,d,g,h,Kwf(ii),Kwf(ii),u);
    trainingTime(4) = toc;

    % Adaptive recursive hybrid myriad filter (RHMy filter)
    tic;
    [g4,h4] = adaptiveRHMy(s1,d,g,h,Khf(ii),u);
    trainingTime(5) = toc;
    
    % Adaptive scaled recursive hybrid myriad filter (SRHMy filter)
    tic;
    [g6,h6,K16] = adaptiveSRHMy(s1,d,g,h,Khf(ii),u);
    trainingTime(6) = toc;

    disp('--------------------------');
    disp(['Filter   Training time (s)'])
    disp('--------------------------');
    disp(['WMy      ' num2str(trainingTime(2),2)]);
    disp(['RWM      ' num2str(trainingTime(1),2)]);
    disp('--------------------------');
    disp(['RWMy     ' num2str(trainingTime(3),2)]);
    disp(['SRWMy    ' num2str(trainingTime(4),2)]);
    disp(['RHMy     ' num2str(trainingTime(5),2)]);
    disp(['SRHMy    ' num2str(trainingTime(6),2)]);
    disp('--------------------------');
    
    mae = zeros(9,trials);
    mse = zeros(9,trials);
    snr = zeros(9,trials);   
    tic;
    disp(['Number of Iterations for the Test Stage: ' num2str(trials)]);
    for jj = 1:trials
        disp(['Iteration:   ' num2str(jj)]);
        if (ii <= 4)
            sn = s + astable(1,length(s),alpha(ii),0,dispersion,0);
        else
            variance = 2*dispersion;
            sn = s + laprnd(1,length(s), 0, sqrt(variance));
        end
        
        % FIR filtering on the noisy chirp signal 96-tap
        d2 = filter(wfir1,1,sn);
        ax2 = d2(1:47); %Revisar fvtool...
        d2=[d2(47+1:end) ax2];
        
        % IIR output
        siir = filter(a,b,sn);
        ax = siir(1:45); %Revisar fvtool...
        siir=[siir(45+1:end) ax];
        
        % Weighted Myriad Filter on the chirp signal
        y1 = wmyFilter(sn,w1,1);
        
        % Recursive Weighted Myriad Filter on the chirp signal
        y2 = rwmyFilter(sn,g2,h2,Kwf(ii),Kwf(ii));
        
        % Scaled Recursive Weighted Myriad Filter on the chirp signal
        y3 = srwmyFilter(sn,go,ho,K1,K2);
        
        % Recursive Hybrid Myriad Filter on the chirp signal
        y4 = rhmyFilter(sn,g4,h4,Khf(ii));
        
        % Scale Recursive Hybrid Myriad Filter on the chirp signal
        y5 = srhmyFilter(sn,g6,h6,K16);
        
        % Dual Weighted Iterative Truncate Mean Filter (DWITM filter)
        si      = [zeros(1,win1-1) sn];
        for i = 1:length(sn)
            samples_dwitm   = si(i:i+win1-1);
            weights_dwitm   = wfir1;
            y6(i)           = dwitm(samples_dwitm, weights_dwitm);
        end
        aux_y6  = y6(1:47); %Revisar fvtool...
        y6      = [y6(47+1:end) aux_y6];
        
        % Recursive Weighted Median Filter on the chirp signal
        y7 = rwmFilter(sn,gwm,hwm);
        
        %--------------------------------------------
        % Computing metrics
        %--------------------------------------------
        
        % Mean Absolute Error (MAE)
        mae(1,jj) = mean(abs(d1-d2));
        mae(2,jj) = mean(abs(d1-siir));
        mae(3,jj) = mean(abs(d1'-y1));
        mae(4,jj) = mean(abs(d1-y2));
        mae(5,jj) = mean(abs(d1'-y3));
        mae(6,jj) = mean(abs(d1-y4));
        mae(7,jj) = mean(abs(d1'-y5));
        mae(8,jj) = mean(abs(d1-y6));
        mae(9,jj) = mean(abs(d1-y7));
        
        % Mean Square Error (MSE)
        mse(1,jj) = mean((d1-d2).^2);
        mse(2,jj) = mean((d1-siir).^2);
        mse(3,jj) = mean((d1'-y1).^2);
        mse(4,jj) = mean((d1-y2).^2);
        mse(5,jj) = mean((d1'-y3).^2);
        mse(6,jj) = mean((d1-y4).^2);
        mse(7,jj) = mean((d1'-y5).^2);
        mse(8,jj) = mean((d1-y6).^2);
        mse(9,jj) = mean((d1-y7).^2);
        
        % Signal to noise ratio (SNR)
        snr(1,jj) = 10 * log10(mean(d1.^2)/mean((d1-d2).^2));
        snr(2,jj) = 10 * log10(mean(d1.^2)/mean((d1-siir).^2));
        snr(3,jj) = 10 * log10(mean(d1.^2)/mean((d1'-y1).^2));
        snr(4,jj) = 10 * log10(mean(d1.^2)/mean((d1-y2).^2));
        snr(5,jj) = 10 * log10(mean(d1.^2)/mean((d1'-y3).^2));
        snr(6,jj) = 10 * log10(mean(d1.^2)/mean((d1-y4).^2));
        snr(7,jj) = 10 * log10(mean(d1.^2)/mean((d1'-y5).^2));
        snr(8,jj) = 10 * log10(mean(d1.^2)/mean((d1-y6).^2));
        snr(9,jj) = 10 * log10(mean(d1.^2)/mean((d1-y7).^2));
    end
    disp('--------------------------');
    toc;
    disp('--------------------------');
    
    mae_mean(:,ii) = mean(mae,2);
    mse_mean(:,ii) = mean(mse,2);
    snr_mean(:,ii) = mean(snr,2);
end

disp('Mean absolute error');
disp('--------------------------------------------------------------');
disp(['Filter   0.75       1.00      1.50       2.00       Laplacian'])
disp('--------------------------------------------------------------');
disp(['FIR      ' num2str(mae_mean(1,:),4)]);
disp(['IIR      ' num2str(mae_mean(2,:),4)]);
disp(['WMy      ' num2str(mae_mean(3,:),4)]);
disp(['DWITM    ' num2str(mae_mean(8,:),4)]);
disp(['RWM      ' num2str(mae_mean(9,:),4)]);
disp('--------------------------------------------------------------');
disp(['RWMy     ' num2str(mae_mean(4,:),4)]);
disp(['SRWMy    ' num2str(mae_mean(5,:),4)]);
disp(['RHMy     ' num2str(mae_mean(6,:),4)]);
disp(['SRHMy    ' num2str(mae_mean(7,:),4)]);
disp('--------------------------------------------------------------');