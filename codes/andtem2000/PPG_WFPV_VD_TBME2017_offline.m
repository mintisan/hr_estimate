% Authors: Andriy Temko, INFANT Research Centre
% http://www.infantcentre.ie/
% http://eleceng.ucc.ie/~andreyt/


clear; %close all;

% Dataset IDs
IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

% Overall parameters
srate = 125;    % original sampling rate
FFTres = 1024;   % FFT resolution - user tunable parm1
WFlength = 15; % Wiener filter length - - user tunable parm2
allD = size(IDData,2); % number of recordings

% Filter parameters
CutoffFreqHzHP = 1; % 60 BPM;
CutoffFreqHzLP = 3;% 180 BPM
[b,a] = butter(4, [0.4 4]/(125/2),'bandpass');

% metrics
fullBPM=[];fullBPM0=[]; % for computing the correlation
myError = zeros(1,allD); % Beat per Minute errors
myErrorN = zeros(1,allD); % Beat per Minute errors
myErrorStd = zeros(1,allD); % std error
myRelError = zeros(1,allD); % relative error

% framework rule, 8s window 2s shift, no look into future
window   = 8 * srate;  % window length is 8 seconds
step     = 2 * srate;  % step size is 2 seconds

% loop for each recording
for idnb = 1 : allD
    
    % load the data
    load(['Data/' IDData{idnb}]);
    if idnb>13
        ch1 = 1; ch2 = 2; ch3 = 3; ch4 = 4; ch5 = 5;
    else
        ch1 = 2; ch2 = 3; ch3 = 4; ch4 = 5; ch5 = 6;
    end
    windowNb = floor((length(sig)-window)/step) + 1;  % total number of windows(estimates)
    
    % initialization of variables
    BPM_est = zeros(1,windowNb); % estimated BPM
    BPM_est_N = zeros(1,windowNb); % estimated BPM offline
    rangeIdx = []; % range of search for the next estimates
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN FreqRangePPG;
    
    
    for i =  [1 :  windowNb]
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curData = sig(:,curSegment);
        
        PPG1 = curData(ch1,:); PPG2 = curData(ch2,:);
        ACC_X = curData(ch3,:); ACC_Y = curData(ch4,:); ACC_Z = curData(ch5,:);
        
        % filtering
        PPG1 = filter(b,a,PPG1);
        PPG2 = filter(b,a,PPG2);
        ACC_X = filter(b,a,ACC_X);
        ACC_Y = filter(b,a,ACC_Y);
        ACC_Z = filter(b,a,ACC_Z);
        PPG_ave = 0.5*(PPG1-mean(PPG1))/(std(PPG1))+0.5*(PPG2-mean(PPG2))/(std(PPG2)); % mean overall
        
        % downsampling to 25Hz
        PPG_ave = downsample(PPG_ave,5);
        ACC_X = downsample(ACC_X,5);
        ACC_Y = downsample(ACC_Y,5);
        ACC_Z = downsample(ACC_Z,5);
        srate = 25; % new sampling rate
        
        % Periodogram
        PPG_ave_FFT = fft(PPG_ave,FFTres);
        FreqRange = linspace(0,srate,size(PPG_ave_FFT,2));
        
        % finding the indices for the range of interest
        [extra,lowR] = (min(abs(FreqRange-CutoffFreqHzHP)));
        [extra,highR] = (min(abs(FreqRange-CutoffFreqHzLP)));
        
        %  Getting rid of most spectra outside the range of interest
        FreqRange = FreqRange(lowR:highR);
        PPG_ave_FFT = PPG_ave_FFT(lowR:highR);
        ACC_X_FFT= fft(ACC_X,FFTres); ACC_X_FFT = ACC_X_FFT(lowR:highR);
        ACC_Y_FFT= fft(ACC_Y,FFTres); ACC_Y_FFT = ACC_Y_FFT(lowR:highR);
        ACC_Z_FFT= fft(ACC_Z,FFTres); ACC_Z_FFT = ACC_Z_FFT(lowR:highR);
        
        
        % phase vocoder to refine spectral estimations
        FreqRangePPG(i,:) = FreqRange;
        if i>1 % start phase vocoder for current and previous frames
            for ii=1:size(FreqRangePPG,2)
                curPhase = angle(PPG_ave_FFT(ii));
                prevPhase = angle(PPG_ave_FFTpr(ii)); vocoder = zeros(1,20);
                for n = 1:20
                    vocoder(n) = ((curPhase-prevPhase)+(2*pi*(n-1)))/(2*pi*2);
                end
                difference = vocoder - FreqRange(ii);
                [extra, deltaidx] = min(abs(difference));
                FreqRangePPG(i,ii) = vocoder(deltaidx);
            end
        end
        
        % smooth phase vocoder frequency estimates
        FreqRangePPG(i,:) = moving(FreqRangePPG(i,:),3);
        
        % save previous spectrum for the next phase vocoder call
        PPG_ave_FFTpr = PPG_ave_FFT;
        
        
        % Wiener filtering PPG-ACC, two types
        % W1 filter = clean_signal/(clean_signal+noise) = (all_signal-noise)/all_signal,
        % where all_signal = time_average(ppg_signal), noise = time_average(accelerometer)
        %
        % W2 filter = clean_signal/(clean_signal+noise),
        % where clean_signal = time_average(clean_ppg_signal), noise =
        % time_average(accelerometer)
        %
        % The process of Wiener filtering:
        % 1. Estimate the 'signal'/'clean_signal' level as average over
        %    past 15 spectral envelopes (~30s)
        % 2. Scale the spectrum of noise and signal/clean_signal
        % 3. Apply the filter as weighting of the spectral components
        
        WC1 = WFlength; WC2 = WFlength; % order of the filter for WF1 and WF2
        
        %Wiener 1 / abs & normalised
        W1_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
        if i==1, W1_PPG_ave_FFT_ALL = W1_FFTi(i,:); else W1_PPG_ave_FFT_ALL = mean(W1_FFTi(max(1,i-WC1):i,:),1); end
        W1_PPG_ave_FFT_ALL_norm = (W1_PPG_ave_FFT_ALL)/max(W1_PPG_ave_FFT_ALL);
        W1_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
        W1_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
        W1_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
        WF1 = (1 - 1/3*(W1_ACC_X_FFT_norm+W1_ACC_Y_FFT_norm+W1_ACC_Z_FFT_norm)./(W1_PPG_ave_FFT_ALL_norm));
        WF1 (WF1<0) = -1; % limit negative -inf to -1
        W1_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF1;
        
        % Wiener 1, power2, works better on its own but not in ensemble
        W11_FFTi(i,:) = abs(PPG_ave_FFT).^2;
        if i==1, W11_PPG_ave_FFT_ALL = W11_FFTi(i,:); else W11_PPG_ave_FFT_ALL = mean(W11_FFTi(max(1,i-WC1):i,:),1); end
        W11_PPG_ave_FFT_ALL_norm = (W11_PPG_ave_FFT_ALL)/max(W11_PPG_ave_FFT_ALL);
        W11_ACC_X_FFT_norm = (abs(ACC_X_FFT).^2)/max(abs(ACC_X_FFT.^2));
        W11_ACC_Y_FFT_norm = (abs(ACC_Y_FFT).^2)/max(abs(ACC_Y_FFT).^2);
        W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT).^2)/max(abs(ACC_Z_FFT).^2);
        WF11 = (1 - 1/3*(W11_ACC_X_FFT_norm+W11_ACC_Y_FFT_norm+W11_ACC_Z_FFT_norm)./(W11_PPG_ave_FFT_ALL_norm));
        W11_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*(WF11);
        
        %Wiener 2, abs & normalised
        W2_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
        if i==1, W2_PPG_ave_FFT_ALL = W2_FFTi(i,:); else W2_PPG_ave_FFT_ALL = mean(W2_FFTi(max(1,i-WC2):i,:),1); end
        W2_PPG_ave_FFT_ALL_norm = (W2_PPG_ave_FFT_ALL)/max(W2_PPG_ave_FFT_ALL);
        W2_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
        W2_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
        W2_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
        WF2 = W2_PPG_ave_FFT_ALL_norm./(((W2_ACC_X_FFT_norm+W2_ACC_Y_FFT_norm+W2_ACC_Z_FFT_norm)/3)+W2_PPG_ave_FFT_ALL_norm);
        W2_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF2;
        W2_FFTi(i,:) = (W2_PPG_ave_FFT_Clean(i,:))/max(W2_PPG_ave_FFT_Clean(i,:));
        
        %Wiener 2, power2, work better on its own, but not in ensemble
        W21_FFTi(i,:) = abs(PPG_ave_FFT).^2;
        if i==1, W21_PPG_ave_FFT_ALL = W21_FFTi(i,:); else W21_PPG_ave_FFT_ALL = mean(W21_FFTi(max(1,i-15):i,:),1); end
        W21_PPG_ave_FFT_ALL_norm = W21_PPG_ave_FFT_ALL/max(W21_PPG_ave_FFT_ALL);
        W21_ACC_X_FFT_norm = abs(ACC_X_FFT).^2/max(abs(ACC_X_FFT).^2);
        W21_ACC_Y_FFT_norm = abs(ACC_Y_FFT).^2/max(abs(ACC_Y_FFT).^2);
        W21_ACC_Z_FFT_norm = abs(ACC_Z_FFT).^2/max(abs(ACC_Z_FFT).^2);
        WF21 = W21_PPG_ave_FFT_ALL_norm./(((W21_ACC_X_FFT_norm+W21_ACC_Y_FFT_norm+W21_ACC_Z_FFT_norm)/3)+W21_PPG_ave_FFT_ALL_norm);
        W21_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF21;
        W21_FFTi(i,:) = W21_PPG_ave_FFT_Clean(i,:).^2;
        
        
        W1_PPG_ave_FFT_Clean(i,:) = W1_PPG_ave_FFT_Clean(i,:)/std(W1_PPG_ave_FFT_Clean(i,:)); %4.58BPM
        W11_PPG_ave_FFT_Clean(i,:) = W11_PPG_ave_FFT_Clean(i,:)/std(W11_PPG_ave_FFT_Clean(i,:)); %4.22 BPM, not used
        W2_PPG_ave_FFT_Clean(i,:) = W2_PPG_ave_FFT_Clean(i,:)/std(W2_PPG_ave_FFT_Clean(i,:)); %4.14BPM
        W21_PPG_ave_FFT_Clean(i,:) = W21_PPG_ave_FFT_Clean(i,:)/std(W21_PPG_ave_FFT_Clean(i,:)); %2.42 BPM, not used
        
        PPG_ave_FFT_FIN(i,:) = W1_PPG_ave_FFT_Clean(i,:)+ W2_PPG_ave_FFT_Clean(i,:) ;
        
    end
    
    % comparison with the ground truth
    if idnb>13
        load(['Data/True' IDData{idnb}(5:end)]);           % load groundtruth
    else
        load(['Data/' IDData{idnb} '_BPMtrace']);           % load groundtruth
    end
    
    trans_full = trans_prob(setdiff(1:23,idnb),FreqRange); 
    trans_full = ((moving(trans_full',4))');
    figure; myfS = 8; set(gca,'FontSize',myfS);
    imagesc(FreqRange(end:-1:1)*60,FreqRange(end:-1:1)*60,log(trans_full)); set(gca,'YDir','normal'); set(gca,'FontSize',8); set(gca,'YDir','normal'); 
    xlabel ('Heart Rate (BPM)','FontSize',myfS); ylabel ('Heart Rate (BPM)','FontSize',myfS); 
   
    %figure; imagesc(FreqRange(end:-1:1)*60,FreqRange(end:-1:1)*60,log(trans_full)); set(gca,'YDir','normal');
    %figure; imagesc(FreqRange(end:-1:1)*60,FreqRange(end:-1:1)*60,log(trans_full)); set(gca,'YDir','normal');
    em = PPG_ave_FFT_FIN(:,:)';
    
    %figure; imagesc(1:size(em,2),FreqRange*60,em); set(gca,'YDir','normal'); hold on; plot(BPM0,'w-'); plot(BPM_est,'o')
    [aaa,bbb] = viterbi_path2(trans_full,em);% figure; imagesc(1:size(em,2),FreqRange*60,em); set(gca,'YDir','normal'); hold on; plot(moving(FreqRange(aaa)*60,4),'o-black');
    
    for xxx = 1:size(PPG_ave_FFT_FIN,1),
       BPM_est(xxx) = 60*FreqRangePPG(xxx,aaa(xxx));
    end
    
    save(['Result_' IDData{idnb} '_DATA_WFPV_VD'],'BPM_est', 'BPM0'); 

    % Code to compute the error
    myErrorN(idnb) = mean(abs(BPM0 - moving(BPM_est(1:1:end)',4))); 
    myRelError(idnb) = mean(abs(BPM0 - moving(BPM_est(1:1:end)',4))./BPM0);
    myErrorStd(idnb) = std(abs(BPM0 - moving(BPM_est(1:1:end)',4)));
    
end
fprintf('\nErr12=%2.2f(%2.2f), Err11=%2.2f(%2.2f), ErrAll=%2.2f(%2.2f)', mean(myErrorN(1:12)),mean(myErrorStd(1:12)),mean(myErrorN(13:end)),mean(myErrorStd(13:end)),mean(myErrorN),mean(myErrorStd));

