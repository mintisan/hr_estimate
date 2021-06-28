% https://www.mathworks.com/matlabcentral/fileexchange/70032-robust-empirical-mode-decomposition-remd?s_tid=FX_rc3_behav

clc; clear; close all;
fs = 10000; % sampling frequency
N = 30000; % data amount
t = (1:N)/fs; % time vector
x = (2+cos(2*pi*0.5*t)).*cos(2*pi*5*t+15*t.^2)+...
     cos(2*pi*2*t);
[imf,ort,fvs,iterNum] = emd_sssc(x,fs,'display',1);
