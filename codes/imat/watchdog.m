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
function [ bul,N ] = watchdog( spec ,Nprv,nsp)
bul = 0;
N = 0;
dstns = 9;
trsh = 0.3; %threshold to find very strong peaks
Nprv = ceil(Nprv);
nd = 5;
y1 = spec(1,250:4000);
y2 = spec(2,250:4000);
[h1,p1] = findpeaks(y1,'SORTSTR','descend');
[h2,p2] = findpeaks(y2,'SORTSTR','descend');
p1 = p1 + 249;
p2 = p2 + 249;
% finding very strong peaks that have strong harmonics too for first ppg channel
if h1(1)*0.33 > h1(3)
if abs(p1(1) - p1(2)/2) < dstns 
    if min(abs(p1(1) - nsp)) < nd || min(abs(p1(2) - nsp)) < nd
    bul = 0;
    else
    N = (p1(1) + p1(2))/3;
    bul = 1;
    end
end
if abs(p1(2) - p1(1)/2) < dstns 
    if min(abs(p1(1) - nsp)) < nd || min(abs(p1(2) - nsp)) < nd
    bul = 0;
    else
    N = (p1(1) + p1(2))/3;
    bul = 1;
    end
end
end
% finding very strong peaks that have strong harmonics too for second ppg chanal
if h2(1)*0.33 > h2(3)
if abs(p2(1) - p2(2)/2) < dstns 
    if min(abs(p2(1) - nsp)) < nd || min(abs(p2(2) - nsp)) < nd
    bul = 0;
    else
    N = (p2(1) + p2(2))/3;
    bul = 1;
    end
end
if abs(p2(2) - p2(1)/2) < dstns 
    if min(abs(p2(1) - nsp)) < nd || min(abs(p2(2) - nsp)) < nd
    bul = 0;
    else
    N = (p2(1) + p2(2))/3;
    bul = 1;
    end
end
end
% if the strong peak is near last bpm the threshold increases 
if abs(p1(1) - Nprv) < 50
    trsh = 0.4;
    if min(abs(p1(1) - p1(2:length(p1))/2)) < 5
        trsh = 0.5;
    end
end

if h1(1)*trsh > h1(2)
    if min(abs(p1(1) - nsp)) < nd % checking that found peak is not very close to noise strong peaks
    bul = 0;
    else
        if p1(1) > 800
        N = p1(1)/2;
        else
        N = p1(1);
        end
         bul = 1;
    end
end
% if the strong peak is near last bpm the threshold increases 
trsh = 0.2;
if abs(p2(1) - Nprv) < 50
    trsh = 0.4;
    if min(abs(p2(1) - p2(2:length(p2))/2)) < 5
        trsh = 0.5;
    end
end


if h2(1)*trsh > h2(2)
    if min(abs(p2(1) - nsp)) < nd   % checking that found peak is not very close to noise strong peaks
    bul = 0;
    else
        if p1(1) > 800
            N = p1(1)/2;
        else
            N = p1(1);
        end
        bul = 1;
    end
end

if min(abs(N - nsp)) < nd
    bul = 0;
end
if abs(N - Nprv) > 100 
    N = Nprv + sign(N - Nprv)*50;
    bul = 0;
end
% if no strong peak is found in 8 second window we check 16 second window 
% the code is similar to 8 second window
if bul == 0
    N = 0;
    dstns = 9;
    trsh = 0.4;
    Nprv = ceil(Nprv);
    nd = 5;
    y3 = spec(3,250:4000);
    y4 = spec(4,250:4000);
    [h3,p3] = findpeaks(y3,'SORTSTR','descend');
    [h4,p4] = findpeaks(y4,'SORTSTR','descend');
    p3 = p3 + 249;
    p4 = p4 + 249;
    if h3(1)*0.33 > h3(3)
        if abs(p3(1) - p3(2)/2) < dstns
            if min(abs(p3(1) - nsp)) < nd || min(abs(p3(2) - nsp)) < nd
                bul = 0;
            else
                N = (p3(1) + p3(2))/3;
                bul = 1;
            end
        end
        if abs(p3(2) - p3(1)/2) < dstns
            if min(abs(p3(1) - nsp)) < nd || min(abs(p3(2) - nsp)) < nd
                bul = 0;
            else
                N = (p3(1) + p3(2))/3;
                bul = 1;
            end
        end
    end
    if h4(1)*0.33 > h4(3)
        if abs(p4(1) - p4(2)/2) < dstns
            if min(abs(p4(1) - nsp)) < nd || min(abs(p4(2) - nsp)) < nd
                bul = 0;
            else
                N = (p4(1) + p4(2))/3;
                bul = 1;
            end
        end
        if abs(p4(2) - p4(1)/2) < dstns
            if min(abs(p4(1) - nsp)) < nd || min(abs(p4(2) - nsp)) < nd
                bul = 0;
            else
                N = (p4(1) + p4(2))/3;
                bul = 1;
            end
        end
    end
    
    if abs(p3(1) - Nprv) < 50
        trsh = 0.4;
        if min(abs(p3(1) - p3(2:length(p3))/2)) < 5
            trsh = 0.5;
        end
    end
    
    if h3(1)*trsh > h3(2)
        if min(abs(p3(1) - nsp)) < nd
            bul = 0;
        else
            if p3(1) > 800
                N = p3(1)/2;
            else
                N = p3(1);
            end
            bul = 1;
        end
    end
    
    trsh = 0.4;
    if abs(p4(1) - Nprv) < 50
        trsh = 0.4;
        if min(abs(p4(1) - p4(2:length(p4))/2)) < 5
            trsh = 0.5;
        end
    end
    
    
    if h4(1)*trsh > h4(2)
        if min(abs(p4(1) - nsp)) < nd
            bul = 0;
        else
            if p3(1) > 800
                N = p3(1)/2;
            else
                N = p3(1);
            end
            bul = 1;
        end
    end
    % checking that the found peak is not  very close to noise strong peaks
    if min(abs(N - nsp)) < nd
        bul = 0;
    end
    if abs(N - Nprv) > 100
        N = Nprv + sign(N - Nprv)*50;
        bul = 0;
    end
end

end