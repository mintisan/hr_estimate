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
function [ BPM,sample,estimate ] = spectral_analysis( sig , Nprv,sample,estimate,n,x2,x3,z2,z3)
reso = 1:4000;  
spec = zeros(4,4000);
w = (n-1)*250+1:(n+3)*250;
noise = abs(fft(sig([4,5,6],w)-mean(sig([4,5,6],w),2)*ones(1,1000),30000,2)).^2;
noise = noise(:,1:4000);
[~, nsp1] = findpeaks(noise(1,:),'SORTSTR','descend');%finding noise strong peaks
[~, nsp2] = findpeaks(noise(2,:),'SORTSTR','descend');
[~, nsp3] = findpeaks(noise(3,:),'SORTSTR','descend');
nsp = [nsp1(1:2),nsp2(1:2),nsp3(1:2)];
%% using imat
if n > 4
    r = 1;
    r = r + [r == 1] - [r == 2];
    f = IMAT(x2',r*2);
    spec(1,:) = f(1:4000);
    f = IMAT(x3',r*2);
    spec(2,:) = f(1:4000);
else
    f = abs(fft(x2,30000)).^2;
    spec(1,:) = f(1:4000);
    
    f = abs(fft(x3,30000)).^2;
    spec(2,:) = f(1:4000);
end
if n > 4
    r = 2;
    r = r + [r == 1] - [r == 2];
    f = IMAT(z2',r*2);
    spec(3,:) = f(1:4000);
    f = IMAT(z3',r*2);
    spec(4,:) = f(1:4000);
else
    f = abs(fft(z2,30000)).^2;
    spec(3,:) = f(1:4000);
    
    f = abs(fft(z3,30000)).^2;
    spec(4,:) = f(1:4000);
end
%% peak selection
BPM = peakselection(Nprv,nsp,sample,estimate,n,spec);
end

