% https://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform
% https://www.ilovematlab.cn/thread-541151-1-1.html
close all;
[x,Fs] = audioread('Hum.wav');
plot_hht(x(1:6000),1/Fs);
