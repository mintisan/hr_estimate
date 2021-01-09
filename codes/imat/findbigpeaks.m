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
function [ pik ] = findbigpeaks( spec, trsh, Nprv,nsp,s1,s2)
locs = zeros(4,3);
hgt = zeros(4,3);
locsh = zeros(4,3);
hgth = zeros(4,3);
locs3d = zeros(4,3);
hgt3d = zeros(4,3);
% if min(abs(Nprv - nsp([1,3,5]))) < 40
%     trsh = 0.9; %initializing threshold
% end
for j = 1:4
        b1 = max(Nprv - s1,200);  
        b2 = max(2*Nprv - s2,400); 
        win = b1:Nprv + s1;    %interval around last BPM found
        whar = b2:2*Nprv + s2;  %interval around second harmonic of last BPM found
        w3d = 3*Nprv - 100:3*Nprv + 100;  %interval around third harmonic of last BPM found
        m = max(spec(j,:));
        [h,p]  = findpeaks(spec(j,win),'SORTSTR','descend');
        [h1,p1]  = findpeaks(spec(j,whar),'SORTSTR','descend');
        [h3,p3]  = findpeaks(spec(j,w3d),'SORTSTR','descend');
        
        p = p + b1-1;
        p1 = p1 + b2-1;
        p3 = p3 + 3*Nprv - 101;
        len = min(length(p),3);
        lenh = min(length(p1),3);
        len3d = min(length(p3),3);
        locs(j,1:len) = p(1:len);
        hgt(j,1:len) = h(1:len);
        locsh(j,1:lenh) = p1(1:lenh);
        hgth(j,1:lenh) = h1(1:lenh);
        locs3d(j,1:len3d) = p3(1:len3d);
        hgt3d(j,1:len3d) = h3(1:len3d);
end
    dmn = zeros(1,4);
    ratio = zeros(1,4);
    dmnh = zeros(1,4);
    ratioh = zeros(1,4);
    dmn3d = zeros(1,4);
    ratio3d = zeros(1,4);
    %strong peak selection
    for j = 1:4
        if hgt(j,1)*trsh > hgt(j,2) && hgt(j,1)>m*0.1
            dmn(j) = locs(j,1);
            ratio(j) = hgt(j,1)/(1+hgt(j,2));
        end
    end
    for j = 1:4
        if hgth(j,1)*trsh > hgth(j,2)
            dmnh(j) = locsh(j,1);
            ratioh(j) = hgth(j,1)/(1+hgth(j,2));
        end
    end
    for j = 1:4
        if hgt3d(j,1)*trsh > hgt3d(j,2) 
            dmn3d(j) = locs3d(j,1);
            ratio3d(j) = hgt3d(j,1)/(1+hgt3d(j,2));
        end
    end
    % checking that peaks found are'nt noise peaks
%     if min(abs(Nprv - nsp([1,3,5]))) < 20
%         dmn = zeros(1,4);
%         for k = 1:4
%             if abs(Nprv - dmnh(k)/2) > 30
%                 dmnh(k) = 0;
%             end
%         end
%     end
    pik = [dmn,dmnh,dmn3d];
end
