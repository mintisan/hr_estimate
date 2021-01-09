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
function[BPM] = peakselection(Nprv,nsp,sample,estimate,n,spec)
[bul,N] = watchdog(spec,Nprv,nsp);
if bul == 1
    estimate(n) = N;
    sample(n) = 1;
else
    e = ceil(Nprv);
    q = findbigpeaks(spec,0.6,e,nsp,60,60);
    q(5:8) = q(5:8)/2;% second harmonic
    q(9:12) = q(9:12)/3; %third harmonic
    jam = 0;
    k = 0;
    for s = 1:8   % mean of first and second haromincs
        if q(s) ~= 0;
            jam = jam + q(s);
            k = k+1;
        end
    end
    estimate(n) = jam /( k + (k==0));
    if jam == 0
        estimate(n) = Nprv;
    end
end
if abs(estimate(n) - Nprv) > 50
    e = ceil(Nprv);
    q = findbigpeaks(spec,0.6,e,nsp,60,60);
    q(5:8) = q(5:8)/2;
    q(9:12) = q(9:12)/3;
    jam = 0;
    k = 0;
    v = zeros(0,0);
    for s = 1:8
        if q(s) ~= 0;
            jam = jam + q(s);
            v = [v,q(s)];
            k = k+1;
        end
    end 
    estimate(n) = jam /( k + (k==0));
    if jam == 0
        estimate(n) = Nprv; % if no dominent peak found 
    end
end
if estimate(n) == 0
    estimate(n) = Nprv;
end
Nprv = estimate(n);
BPM = estimate(n)/4; 
end