%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Codes For: 
%           "Heart Rate Tracking using Wrist-Type
%    Photoplethysmographic (PPG) Signals during Physical
%         Exercise with Simultaneous Accelerometry"
%
% References:
%
% [1] M. Boloursaz, E. Asadi, M. Eskandari, S. Kiani, and F. Marvasti, “Heart Rate Tracking using Wrist-Type Photoplethysmographic (PPG)
%     Signals during Physical Exercise with Simultaneous Accelerometry,” Submitted to IEEE Signal Processing Letters, March 2015
% [2] Farokh Marvasti et al., "A Unified Approach to Sparse Signal Processing," EURASIP Journal on Advances in Signal Processing,2012
%
% Written By: Ehsan Asadi Kangarshahi
%             Mohsen Eskandari
%             Shahrzad Kiani
%             Mahdi Boloursaz Mashhadi
% Affiliation: 
%          Advanced Communications Research Institute (ACRI)
% Electrical Engineering Department, Sharif University of Technology
%                              Tehran, Iran
%
% For any problems, contact me at asadikangarshahi_ehsan@ee.sharif.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
addpath(genpath('..\competition_data'));
sample = zeros(1,200);
estimate = zeros(1,200);
ID = { 'DATA_01_TYPE01', 'DATA_02_TYPE02', 'DATA_03_TYPE02', 'DATA_04_TYPE02', ...
       'DATA_05_TYPE02', 'DATA_06_TYPE02', 'DATA_07_TYPE02', 'DATA_08_TYPE02',...
       'DATA_09_TYPE02', 'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02'};  
ID2  = { 'DATA_01_TYPE01_BPMtrace', 'DATA_02_TYPE02_BPMtrace', 'DATA_03_TYPE02_BPMtrace', 'DATA_04_TYPE02_BPMtrace', ...
       'DATA_05_TYPE02_BPMtrace', 'DATA_06_TYPE02_BPMtrace', 'DATA_07_TYPE02_BPMtrace', 'DATA_08_TYPE02_BPMtrace',...
       'DATA_09_TYPE02_BPMtrace', 'DATA_10_TYPE02_BPMtrace','DATA_11_TYPE02_BPMtrace','DATA_12_TYPE02_BPMtrace'};
for idnb = 1 : 12
    close all;
    load(ID{idnb});
    load(ID2{idnb});%the ground-truth bpm is not used in the algorithm
    srate = 125;                            
    window   = 8 * srate;                    
    step     = 2 * srate;                    
    windowNb = floor((length(sig)-window)/step) + 1;
    BPM = zeros(1,200);
%%
%initialize the HR
    firstwin = sig(:,1:1000);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [t1,t2] = MA_cancellation(firstwin,100,20,1,1000,500);
    t1 = t1';
    t2 = t2';
    t1 = abs(fft(t1,30000)).^2;
    t2 = abs(fft(t2,30000)).^2;
    emt1 = [];
    emt2 = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     plot(t1(1:2000));
    [g11,q11] = findpeaks(t1(200:2000));
    [g22,q22] = findpeaks(t2(200:2000));
    q11 = q11 + 199;
    q22 = q22 + 199;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a1,b1] = sort(g11,'descend');
    [a2,b2] = sort(g22,'descend');
    q1 = q11(b1);
    q2 = q22(b2);
    g1 = a1;
    g2 = a2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (g1(2)/g1(1))<0.3
        if (q1(1)>700)
            q1(1) = q1(1)/2;
        end
        y = q1(1);
      
    elseif ((g1(3)/g1(1))<0.3) && ((abs(q1(1)/2-q1(2))<9) || (abs(q1(2)/2-q1(1))<9))
        h = ceil(q1(1)/q1(2)-0.5);
        y = q1(1) / (h+(h==0));
    
    elseif (g2(2)/g2(1))<0.3
        if (q2(1)>700)
            q2(1) = q2(1)/2;
        end
        y = q2(1);
       
    elseif ((g2(3)/g2(1))<0.3) && ((abs(q2(1)/2-q2(2))<9) || (abs(q2(2)/2-q2(1))<9))
        h = ceil(q2(1)/q2(2)-0.5);
        y = q2(1) / (h+(h==0));
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    else
        g1 = g11;
        g2 = g22;
        q1 = q11;
        q2 = q22;
        for i=1:length(g1)
            temp = [];
            for j=1:length(g1)
                if abs(q1(i)-q1(j))<200
                    temp =[temp,j];
                end
            end
            tg1 = g1(temp);
            emt1 = [emt1,100*g1(i)/max(tg1)];
        end

        for i=1:length(g2)
            temp = [];
            for j=1:length(g2)
                if abs(q2(i)-q2(j))<200
                    temp =[temp,j];
                end
            end
            tg2 = g2(temp);
            emt2 = [emt2,100*g2(i)/max(tg2)];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        k = [];
        we = [q1,q2];
        qr = [emt1,emt2];
        

        for i=1:length(we)
            if qr(i) >= 100
                k = [k,we(i)];
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        moj = 0;
        while (moj==0)
            f = [];
            mini1 = min(k);
            mini = mini1;
            if mini > 700
                p = ceil(mini/420-0.6);
                mini = mini / (p+(p==0));
            end

            for i=1:length(k)
                temp = ceil(k(i)/mini-0.9);
                if temp<3
                    f = [f,k(i)/(temp+(temp==0))]; 
                end
            end

            maxi = 0;
            for i=1:length(f)
                ted = 0;
                for j=1:length(f)
                    if abs(f(i)-f(j))<10
                        ted = ted + 1;
                    end
                end
                if ted > maxi 
                    y = f(i);
                    maxi = ted;
                end
            end
            if (abs(mini-y)<10) && (y > 230) && (maxi~=1)
                moj =1;
            else
                temp2 = [];
                for i=1:length(k)
                    if k(i)~=mini1
                        temp2 = [temp2,k(i)];
                    end
                end
                k = temp2;
            end
            if isempty(k)
                moj = 1;
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
    BPM(1) = y/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%%
Nprv = BPM(1)*4;
figure (1);
stem(BPM0);
hold on;
for n =   2  :  windowNb
    endpoint = (n-1)*step+window; %the algorithm takes the input PPG signal from the beginning to the end of the current window
                                  %i.e it doesn't use recordings of the subsequent windows and hence is casual
    [x2,x3] = MA_cancellation(sig,100,30,n,1000,Nprv);
    [z2,z3] = MA_cancellation(sig,100,30,n,2000,Nprv);
[BPM(n),sample,estimate] = spectral_analysis(sig(:,1:endpoint),Nprv,sample,estimate,n,x2,x3,z2,z3);
Nprv = BPM(n)*4;
figure(1);
stem(BPM(1:n),'r');
legend('true BPM','estimation','location','southoutside');
end
    RES = { 'Result_S01_T01', 'Result_S02_T02', 'Result_S03_T02', 'Result_S04_T02', ...
       'Result_S05_T02', 'Result_S06_T02', 'Result_S07_T01', 'Result_S08_T02',...
       'Result_S09_T02', 'Result_S10_T01','Result_S11_T01','Result_S12_T01'}; 
    save(RES{idnb},'BPM');
   Average_absolute_error = mean(abs(BPM0- BPM(1:size(BPM0))'))
   Average_error_percentage = mean(abs(BPM0- BPM(1:size(BPM0))')./BPM0)
   Estimation_variance = var(BPM0- BPM(1:size(BPM0))')
   Pearson_correlation = corr(BPM0,BPM(1:size(BPM0))')
end