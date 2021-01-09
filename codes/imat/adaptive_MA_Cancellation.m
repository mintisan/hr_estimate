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
function [ x2,x3 ] = adaptive_MA_Cancellation( s1,s2,s3,x2,x3,k,deg,Nprv )
ms = 0.005;
[~,frqs1] = max(abs(fft(s1,30000,1)),[],1);
[~,frqs2] = max(abs(fft(s2,30000,1)),[],1);
[~,frqs3] = max(abs(fft(s3,30000,1)),[],1);
ha1 = adaptfilt.lms(deg,ms);
ha2 = adaptfilt.lms(deg,ms);
ha3 = adaptfilt.lms(deg,ms);
k = min(20,k);
r = 1;
% using adaptive filters on signal
for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs1(i) < 250 || frqs1(i) > 800
        [~,x2] = filter(ha1,s1(:,i),x2);
    else
        if frqs1(i) < Nprv+50 && frqs1(i) > Nprv - 50
            [~,x2] = filter(ha3,s1(:,i)*r,x2);
        else
            [~,x2] = filter(ha2,s1(:,i),x2);
        end
    end
end
for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs2(i) < 250 || frqs2(i) > 800
        [~,x2] = filter(ha1,s2(:,i),x2);
    else
        if frqs2(i) < Nprv+50 && frqs2(i) > Nprv - 50
            [~,x2] = filter(ha3,s2(:,i)*r,x2);
        else
            [~,x2] = filter(ha2,s2(:,i),x2);
        end
    end
end
for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs3(i) < 250 || frqs3(i) > 800
        [~,x2] = filter(ha1,s3(:,i),x2);
    else
        if frqs3(i) < Nprv+50 && frqs3(i) > Nprv - 50
            [~,x2] = filter(ha3,s3(:,i)*r,x2);
        else
            [~,x2] = filter(ha2,s3(:,i),x2);
        end
    end
end


for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs1(i) < 250 || frqs1(i) > 800
        [~,x3] = filter(ha1,s1(:,i),x3);
    else
        if frqs1(i) < Nprv+50 && frqs1(i) > Nprv - 50
            [~,x3] = filter(ha3,s1(:,i)*r,x3);
        else
            [~,x3] = filter(ha2,s1(:,i),x3);
        end
    end
end
for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs2(i) < 250 || frqs2(i) > 800
        [~,x3] = filter(ha1,s2(:,i),x3);
    else
        if frqs2(i) < Nprv+50 && frqs2(i) > Nprv - 50
            [~,x3] = filter(ha3,s2(:,i)*r,x3);
        else
            [~,x3] = filter(ha2,s2(:,i),x3);
        end
    end
end
for i=1:k
    reset(ha1);
    reset(ha2);
    reset(ha3);
    if frqs3(i) < 250 || frqs3(i) > 800
        [~,x3] = filter(ha1,s3(:,i),x3);
    else
        if frqs3(i) < Nprv+50 && frqs3(i) > Nprv - 50
            [~,x3] = filter(ha3,s3(:,i)*r,x3);
        else
            [~,x3] = filter(ha2,s3(:,i),x3);
        end
    end
end
end

