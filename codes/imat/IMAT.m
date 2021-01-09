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
function [f] = IMAT(y,L)
fnum=1;
iternum=5;
length(y);
sig=y;
i=1:(L*250 + [L == 2]*1500);
msk=zeros(1,4000);
msk(L*(i-1)+1)=1;
zsig=zeros(1,4000);
zsig(L*(i-1)+1)=sig(i);
tsig=abs(fft(zsig,4000));
beta=max(tsig)/10;
alpha=0.1;
landa=.2;
xnew=0;
s=0;
for i=1:length(msk)
        if msk(i)==0
            y(i)=s;
        else
            y(i)=zsig(i);
            s=zsig(i);
        end
end
x=zeros(1,length(msk));
for i=1:iternum
    s=0;
    for ii=1:length(msk)
        if msk(ii)==0
            xI(ii)=s;
        else
            xI(ii)=x(ii);
            s=x(ii);
        end
    end
    xnew=x+landa*(y-xI);
    tx=fft(xnew);
    th(i)=beta*exp(-alpha*(i-1));
    ttmsk=zeros(1,length(zsig));
    ttmsk(abs(tx)>th(i))=1;
    x=ifft(ttmsk.*tx);
end
    f = abs(fft(x,120000)).^2;
end