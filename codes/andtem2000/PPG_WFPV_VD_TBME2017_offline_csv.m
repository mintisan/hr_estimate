% Authors: Andriy Temko, INFANT Research Centre
% http://www.infantcentre.ie/
% http://eleceng.ucc.ie/~andreyt/


clear; close all;

pathroot = '..\..\datasets\心率全场景测试数据\';
file = dir(fullfile(pathroot, '*.csv'));
filenum = size(file,1);
filename = {file.name}';

fs = 25;

window_seconds = 7.2 * fs;
step_second = 1 * fs;

ch_acc_x = 1;
ch_acc_y = 2;
ch_acc_z = 3;
ch_ppg = 4;

CutoffFreqHzHP = 1; % 60 BPM;
CutoffFreqHzLP = 4;% 180 BPM
[b,a] = butter(4, [0.4 4]/(fs/2),'bandpass');

FFTres = 1024;

for k=1:filenum
    filename(k)
    % 1. 加载文件和数据
    fullname = strcat(pathroot, filename(k));
    fid=importdata(fullname{1}, ',', 1);
    polar = fid.data(:,2);
    polar = downsample(polar, fs);
    data=fid.textdata;
    data(1,:) = [];     % remove first line
    total_row=size(data,1);
    data(total_row-100:total_row,:) = [];       % remove last abnormal
    
    % 2. 对每个文件进行划窗处理
    windowNb = floor((size(data,1)-window_seconds)/step_second) + 1;  % total number of windows(estimates)
    clear FFT_PPG_specgram
    for i =  [1 :  windowNb]
        curSegment = (i-1)*step_second+1 : (i-1)*step_second+window_seconds;
        curData = data(curSegment, :);
        % load data to array
        ACC_X = str2double(curData(:,ch_acc_x)); 
        ACC_Y = str2double(curData(:,ch_acc_y)); 
        ACC_Z = str2double(curData(:,ch_acc_z));
        PPG = str2double(curData(:,ch_ppg));
        % filter
        PPG = PPG - mean(PPG);
        PPG = diff(PPG);
%         figure(1);
%         plot(PPG);

        
        % Periodogram
        PPG_FFT = abs(fftshift(fft(PPG,FFTres)));
        
        FreqRange = linspace(0,fs,size(PPG_FFT,1));
        % finding the indices for the range of interest
        [extra,lowR] = (min(abs(FreqRange-CutoffFreqHzHP)));
        [extra,highR] = (min(abs(FreqRange-CutoffFreqHzLP)));
        
        %  Getting rid of most spectra outside the range of interest
        FreqRange = FreqRange(1:highR);

%         figure(2);
        f = (FFTres/2:FFTres-1)*fs/FFTres - fs/2;
        y = 2*PPG_FFT(FFTres/2:FFTres-1)/FFTres;
        bpm = f * 60;
%         plot(bpm, y);
%         xlabel('BPM'); 
%         ylabel('Amplitude'); 
        
%         close all;
        
        FFT_PPG_specgram(:,i) = y(1:highR);
    end 
    
    % 3. 保存和显示结果
    figure; myfS = 8; set(gca,'FontSize',myfS);
    BpmRange = FreqRange(1:highR) * 60;
    SeccondRange = (1:windowNb);
    imagesc(SeccondRange, BpmRange, log(FFT_PPG_specgram)); set(gca,'YDir','normal'); set(gca,'FontSize',8);
    hold on; plot(polar(length(polar)-windowNb:end),'r-');
    xlabel ('Second(s)','FontSize',myfS); ylabel ('Heart Rate (BPM)','FontSize',myfS); 
    title(filename(k))
	saveas(gcf,strcat(fullname{1}, '.jpg'))      % save specgram & result to image
    close all;
end



