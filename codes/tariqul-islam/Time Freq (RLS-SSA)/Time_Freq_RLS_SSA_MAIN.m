
clear; clc; close all;
 
addpath(genpath('TestData'));

% Test Dataset IDs
ID = { 'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
       'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
       'TEST_S07_T02', 'TEST_S08_T01'};  
   PD = { 'True_S01_T01', 'True_S02_T01', 'True_S02_T02', 'True_S03_T02', ...
       'True_S04_T02', 'True_S05_T02', 'True_S06_T01', 'True_S06_T02',...
       'True_S07_T02', 'True_S08_T01'};
         
 m_err = zeros(1,1);
 err_std = [];
 ERRK = [];
 BPM_TRUE = [];
 BPM_EST = [];
for idnb = 1
     
    if idnb>13
     load(ID{idnb-13});                          % load test dataset
     load(PD{idnb-13});
    else 
     load([num2str(idnb) '.mat']);
     load([num2str(idnb) 'B.mat']);
     sig = sig(2:end,:);
    end
    
    srate = 125;                             % 125 Hz    
    window   = 8 * srate;                    % window length is 8 seconds
    step     = 2 * srate;                    % step size is 2 seconds
        
    windowNb = (length(sig)-window)/step + 1;  % total number of windows(estimates)
    
    %**********************************************************************
    % Please write your codes as follows (i.e.,inputing data window by window)
    BPM = [];
    N_prev = 0;
    d = 0;
    tic
    for i =   1  : windowNb

        curSegment = (i-1)*step+1 : (i-1)*step+window;
        
        % Your algorithm's code
        %tic
            [BPM(i), N_prev, d]= time_freq_rls_ssa(sig(:,curSegment),BPM,srate,i,idnb, N_prev, BPM0(i), d);
        %toc
    end
    toc

    %**********************************************************************
    
    % Codes to save results
    BPM = BPM';
    err = BPM-BPM0;
    
    m_err(idnb) = mean(abs(err))
    std_err(idnb) = std(err);
    err_std = [err_std; abs(err)];
    
    
    BPM_TRUE = [BPM_TRUE; BPM0];
    BPM_EST = [BPM_EST; BPM];
    ERRK = [ERRK; err];
end

m_err
mean(m_err)
std_err = std(err_std)

