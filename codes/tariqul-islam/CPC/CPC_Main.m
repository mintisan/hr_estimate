clear; clc; close all;
 
% Implements the paper:
% Mohammad Tariqul Islam, Sk. Tanvir Ahmed, Ishmam Zabir, Celia Shahnaz, and Shaikh Anowarul Fattah, 
% "Cascade and Parallel Combination (CPC) of Adaptive Filters for Estimating Heart Rate During Intensive 
% Physical Exercise from Photoplethysmographic Signal", Healthcare Technology Letters, 2017
 
addpath(genpath('TestData'));

% Test Dataset IDs
ID = { 'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
       'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
       'TEST_S07_T02', 'TEST_S08_T01'};        
ID2 = { 'TRUE_S01_T01', 'TRUE_S02_T01', 'TRUE_S02_T02', 'TRUE_S03_T02', ...
       'TRUE_S04_T02', 'TRUE_S05_T02', 'TRUE_S06_T01', 'TRUE_S06_T02',...
       'TRUE_S07_T02', 'TRUE_S08_T01'}; 
 m_err = zeros(1,1);
 m_err_dom = zeros(1,1);
 N_LMS = 1;
 N_RLS = 1;
 see = 0;
 ERRK = [];
 BPM_TRUE = [];
 BPM_EST = [];

for idnb = 1:1
    if idnb>13
     load(ID{idnb-13});                          % load test dataset
     load(ID2{idnb-13});
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
    BPM_DOM = [];
    N_prev  = 0;
    tic
    for i =   1 : windowNb

        curSegment = (i-1)*step+1 : (i-1)*step+window;
        %tic;
        % Your algorithm's code
        [BPM(i), BPM_DOM(i), N_prev]= CPC(sig(:,curSegment),BPM,srate,i,idnb, N_prev, N_LMS, N_RLS,see,BPM0(i),BPM_DOM);
        %toc;
    end
    toc

    %**********************************************************************
    
    % Codes to save results
    BPM = BPM';
    BPM = post_processing(BPM);
    %str = ID{idnb};
%     str2 = ['Saved Data\Result_' num2str(idnb)];
%     save(str2,'BPM');
%     load(['true' str(5:end)]);
    err = BPM-BPM0;
    BPM_TRUE = [BPM_TRUE; BPM0];
    BPM_EST = [BPM_EST; BPM];
    ERRK = [ERRK; err];
    m_err(idnb) = mean(abs(err))
    m_err_dom(idnb) = mean(abs(BPM0-BPM_DOM'))
    
    figure,plot(1:windowNb, [BPM';BPM_DOM;BPM0']);
    legend('Estimated', 'Dominant', 'Original');
    drawnow
end

mean(m_err)
mean(m_err_dom)