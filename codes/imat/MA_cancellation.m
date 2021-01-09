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
function [x2, x3] = MA_cancellation(sig,k,deg,n,L,Nprv)
tah = (n+3)*250;
sar = max(tah - L + 1,1);
win = sar:tah;
tool = tah - sar + 1;
x2 = sig(2,win);
x3 = sig(3,win);
y1 = refMA_generation_using_svd(sig(4,:),k,n,tool); %refMA_generation_using_svd
y2 = refMA_generation_using_svd(sig(5,:),k,n,tool);
y3 = refMA_generation_using_svd(sig(6,:),k,n,tool);
s1 = reshape(y1,tool,k);
s2 = reshape(y2,tool,k);
s3 = reshape(y3,tool,k);

x2 = x2 - mean(x2);
x3 = x3 - mean(x3);

en = ((mean(x2.^2)*2)^0.5)*2;
for q = 1:tah-sar+1
    if abs(x2(q)) > en
        x2(q) = (en)*sign(x2(q));
    end
end

en = ((mean(x3.^2)*2)^0.5)*2;
for q = 1:tah-sar+1
    if abs(x3(q)) > en
        x3(q) = (en)*sign(x3(q));
    end
end

x2 = reshape(x2,tool,1);
x3 = reshape(x3,tool,1);
[x2,x3] = adaptive_MA_Cancellation(s1,s2,s3,x2,x3,k,deg,Nprv); %Adaptive_MA_cancellation
end
%%%%%%%%%%%%%%%%%%%%
