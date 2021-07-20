%================ Reconstruct parameters of the component =================
%============== from its time-frequency support in WFT or WT ==============
% Version 1.01 stable
%-------------------------------Copyright----------------------------------
%
% Author: Dmytro Iatsenko
% Information about these codes (e.g. links to the Video Instructions),
% as well as other MatLab programs and many more can be found at
% http://www.physics.lancs.ac.uk/research/nbmphysics/diats/tfr
% 
% Related articles:
% [1] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part I: Overview, standards of use, related issues and algorithms."
% {preprint:arXiv:1310.7215}
% [2] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part II: Resolution, reconstruction and concentration."
% {preprint:arXiv:1310.7274}
%
%------------------------------Documentation-------------------------------
%
% [iamp,iphi,ifreq,Optional:rtfsupp] = rectfr(tfsupp,TFR,freq,wopt,Optional:method)
% - returns the component's amplitude [iamp], phase [iphi] and frequency
%   [ifreq] as reconstructed from its extracted time-frequency support
%   [tfsupp] in the signal's WFT/WT [TFR] (determines whether TFR is WFT
%   or WT based on the spacings between specified frequencies [freq] -
%   linear (WFT) or logarithmic (WT)). The optional output [rtfsupp]
%   returns the extracted time-frequency support if the input [tfsupp]
%   specifies 1xL frequency profile instead of the full time-frequency
%   support (see below); otherwise returns input [rtfsupp]=[tfsupp].
%
% INPUT:
% tfsupp: 3xL matrix  (or 1xL vector of the desired frequency profile)
%        - extracted time-frequency support of the component, containing
%          frequencies of the TFR amplitude peaks (ridge points) in the
%          first row, support lower bounds (referred as \omega_-(t)/2/pi
%          in [1]) - in the second row, and the upper bounds (referred as
%          \omega_+(t)/2/pi in [1]) - in the third row. Alternatively, one
%          can specify [tfsupp] as 1xL vector of the desired frequency
%          profile, in which case the program will automatically select
%          time-frequency support around it and the corresponding peaks.
% TFR: NFxL matrix (rows correspond to frequencies, columns - to time)
%        - time-frequency representation (WFT or WT), to which [tfsupp]
%          correspond
% freq: NFx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by function wft.m or wt.m
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions wft, wt; [wopt] contains all the
%          needed information, i.e. name of the TFR, sampling frequency,
%          window/wavelet characteristics etc.
% method:'direct'(default)|'ridge'|'both'
%        - the reconstruction method to use for estimating the component's
%          parameters [iamp], [iphi], [ifreq] (see [1]); if set to 'both',
%          all parameters are returned as 2xL matrices with direct and
%          ridge estimates corresponding to 1st and 2nd rows, respectively.
%
% NOTE: in the case of direct reconstruction, if the window/wavelet does
% not allow direct estimation of frequency, i.e. [wopt.wp.omg=Inf] for
% the WFT or [wopt.wp.D=Inf] for the WT (as for the Morlet wavelet),
% corresponding to infinite \bar{\omega}_g or D_\psi (in the notation of
% [1]), then the frequency is reconstructed by hybrid method (see [1]).
%
%-------------------------------Examples-----------------------------------
%
% /Calculate first WFT (see wft function):/
% [WFT,freq,wopt]=wft(sig,fs,...(parameters));
% /Extract somehow the component time-frequency support from it, e.g./
% tfsupp=ecurve(WFT,freq,wopt);
% /Then component's parameters can be reconstructed from its extracted
% WFT support [tfsupp] as/
% [iamp,iphi,ifreq]=rectfr(tfsupp,WFT,freq,wopt);
% /For WT the same procedure is used./
%
%------------------------------Changelog-----------------------------------
% 
% v1.01:
% - some minor changes plus better optimization
% - added possibility to specify [tfsupp] as frequency profile
% - added possibility to specify reconstruction method
% - added some other possibilities
% - for simplicity, only the direct estimates are now returned by default
% 
%--------------------------------------------------------------------------


function [iamp,iphi,ifreq,varargout] = rectfr(tfsupp,TFR,freq,wopt,varargin)

[NF,L]=size(TFR); freq=freq(:); fs=wopt.fs; wp=wopt.wp; method=1;
if nargin>4 && ~isempty(varargin{1})
    if strcmpi(varargin{1},'direct'), method=1;
    elseif strcmpi(varargin{1},'ridge'), method=2;
    else method=0; end
end
idt=1:L; if iscell(tfsupp), idt=tfsupp{2}; tfsupp=tfsupp{1}; end

%Define component parameters and find time-limits
NR=2; if method==1 || method==2, NR=1; end, NC=length(idt);
mm=ones(NR,NC)*NaN; ifreq=mm; iamp=mm; iphi=mm; asig=mm;
tn1=find(~isnan(tfsupp(1,:)),1,'first'); tn2=find(~isnan(tfsupp(1,:)),1,'last');

%Determine the frequency resolution and transform [tfsupp] to indices
if freq(1)<=0 || std(diff(freq))<std(freq(2:end)./freq(1:end-1))
    fres='linear'; fstep=mean(diff(freq));
    tfsupp=1+floor((1/2)+(tfsupp-freq(1))/fstep);
else
    fres='log'; fstep=mean(diff(log(freq)));
    tfsupp=1+floor((1/2)+(log(tfsupp)-log(freq(1)))/fstep);
end
tfsupp(tfsupp<1)=1; tfsupp(tfsupp>NF)=NF;

%Define variables to not use cells or structures
C=wp.C; ompeak=wp.ompeak; fwtmax=wp.fwtmax;
if isfield(wp,'omg'), omg=wp.omg; else D=wp.D; end
if ~iscell(wp.fwt), fwt=wp.fwt; nflag=0;
else
    nflag=1; Lf=length(wp.fwt{2});
    if strcmpi(fres,'linear'), fxi=wp.fwt{2}; fwt=wp.fwt{1};
    else
        fxi=linspace(min(wp.fwt{2}),max(wp.fwt{2}),Lf)';
        fwt=interp1(wp.fwt{2},wp.fwt{1},fxi,'spline');
    end
    wstep=mean(diff(fxi));
end

%If only the frequency profile is specified, extract the full time-frequency support
if min(size(tfsupp))==1
    eind=tfsupp(:)'; tfsupp=zeros(3,L)*NaN;
    for tn=tn1:tn2
        cind=eind(tn); xn=idt(tn); cs=abs(TFR(:,xn));
        
        %Ridge point
        cpeak=cind;
        if cind>1 && cind<NF
            if cs(cind+1)==cs(cind-1)
                cpeak1=cind-1+find(cs(cind:end-1)>=cs(cind-1:end-2) & cs(cind:end-1)>cs(cind+1:end),1,'first'); cpeak1=min([cpeak1,NF]);
                cpeak2=cind+1-find(cs(cind:-1:2)>=cs(cind+1:-1:3) & cs(cind:-1:2)>cs(cind-1:-1:1),1,'first'); cpeak2=max([cpeak2,1]);
                if cs(cpeak1)>0 && cs(cpeak2)>0
                    if cpeak1-cind==cind-cpeak2
                        if cs(cpeak1)>cs(cpeak2), cpeak=cpeak1;
                        else cpeak=cpeak2; end
                    elseif cpeak1-cind<cind-cpeak2, cpeak=cpeak1;
                    elseif cpeak1-cind>cind-cpeak2, cpeak=cpeak2;
                    end
                elseif cs(cpeak1)==0, cpeak=cpeak2;
                elseif cs(cpeak2)==0, cpeak=cpeak1;
                end
            elseif cs(cind+1)>cs(cind-1)
                cpeak=cind-1+find(cs(cind:end-1)>=cs(cind-1:end-2) & cs(cind:end-1)>cs(cind+1:end),1,'first'); cpeak=min([cpeak,NF]);
            elseif cs(cind+1)<cs(cind-1)
                cpeak=cind+1-find(cs(cind:-1:2)>cs(cind-1:-1:1) & cs(cind:-1:2)>=cs(cind+1:-1:3),1,'first'); cpeak=max([cpeak,1]);
            end
        elseif cind==1
            if cs(2)<cs(1), cpeak=cind;
            else
                cpeak=1+find(cs(cind+1:end-1)>=cs(cind:end-2) & cs(cind+1:end-1)>cs(cind+2:end),1,'first'); cpeak=min([cpeak,NF]);
            end
        elseif cind==NF
            if cs(NF-1)<cs(NF), cpeak=cind;
            else
                cpeak=NF-find(cs(cind-1:-1:2)>cs(cind-2:-1:1) & cs(cind-1:-1:2)>=cs(cind:-1:3),1,'first'); cpeak=max([cpeak,1]);
            end
        end
        tfsupp(1,tn)=cpeak;
        
        %Boundaries of time-frequency support
        iup=[]; idown=[];
        if cpeak<NF-1, iup=cpeak+find(cs(cpeak+1:end-1)<=cs(cpeak:end-2) & cs(cpeak+1:end-1)<cs(cpeak+2:end),1,'first'); end
        if cpeak>2, idown=cpeak-find(cs(cpeak-1:-1:2)<=cs(cpeak:-1:3) & cs(cpeak-1:-1:2)<cs(cpeak-2:-1:1),1,'first'); end
        iup=min([iup,NF]); idown=max([idown,1]);
        tfsupp(2,tn)=idown; tfsupp(3,tn)=iup;
    end
end
%If only the boundaries of the frequency profile are specified, extract the peaks
if min(size(tfsupp))==2
    pind=zeros(1,L)*NaN;
    for tn=tn1:tn2
        ii=tfsupp(1,tn):tfsupp(2,tn); xn=idt(tn); cs=abs(TFR(ii,xn));
        [~,mid]=max(cs); pind(tn)=ii(mid);
    end
    tfsupp=[pind;tfsupp];
end
%Return extracted time-frequency support if requested
if nargout>3, varargout{1}=tfsupp*NaN; varargout{1}(:,tn1:tn2)=freq(tfsupp(:,tn1:tn2)); end

%==================================WFT=====================================
if strcmpi(fres,'linear')
    
    %Direct reconstruction-------------------------------------------------
    if method==0 || method==1
        for tn=tn1:tn2
            ii=tfsupp(2,tn):tfsupp(3,tn); xn=idt(tn); cs=TFR(ii,xn);
            if ~isfinite(omg)
                if xn>idt(tn1) && xn<idt(tn2), cw=angle(TFR(ii,xn+1))-angle(TFR(ii,xn-1)); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs/2;
                elseif xn==1, cw=angle(TFR(ii,xn+1))-angle(cs); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs;
                else cw=angle(cs)-angle(TFR(ii,xn-1)); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs; end
                cw=cw/(2*pi);
            end
            
            if isempty(cs(isnan(cs) | ~isfinite(cs)))
                casig=(1/C)*sum(cs*(2*pi*fstep)); asig(1,tn)=casig;
                if isfinite(omg)
                    ifreq(1,tn)=-omg+(1/C)*sum(freq(ii).*cs*(2*pi*fstep))/casig;
                else
                    ifreq(1,tn)=(1/C)*sum(cw.*cs*(2*pi*fstep))/casig;
                end
            end
        end
        ifreq(1,:)=real(ifreq(1,:));
    end
    
    %Ridge reconstruction--------------------------------------------------
    if method==0 || method==2
        rm=2; if method==2, rm=1; end
        for tn=tn1:tn2
            ipeak=tfsupp(1,tn); xn=idt(tn);
            if ipeak>1 && ipeak<NF, cs=TFR([ipeak;ipeak-1;ipeak+1],xn);
            else cs=TFR(ipeak,xn); end
                
            if isfinite(cs(1))
                ifreq(rm,tn)=freq(ipeak)-ompeak/2/pi;
                if ipeak>1 && ipeak<NF %quadratic interpolation
                    a1=abs(cs(2)); a2=abs(cs(1)); a3=abs(cs(3));
                    p=(1/2)*(a1-a3)/(a1-2*a2+a3);
                    if abs(p)<=1, ifreq(rm,tn)=ifreq(rm,tn)+p*fstep; end
                end
                ximax=2*pi*(freq(ipeak)-ifreq(rm,tn));
                if nflag==0 %if window FT is known in analytic form
                    cmax=fwt(ximax); if isnan(cmax), cmax=fwt(ximax+10^(-14)); end
                    if isnan(cmax), cmax=fwtmax; end
                else %if window FT is numerically estimated
                    cid1=1+floor((ximax-fxi(1))/wstep); cid2=cid1+1;
                    cid1=min([max([cid1,1]),Lf]); cid2=min([max([cid2,1]),Lf]);
                    if cid1==cid2, cmax=fwt(cid1);
                    else
                        cmax=fwt(cid1)+(fwt(cid2)-fwt(cid1))*(ximax-fxi(cid1))/(fxi(cid2)-fxi(cid1));
                    end
                end
                casig=2*cs(1)/cmax; asig(rm,tn)=casig;
            end
        end
    end
    
end

%===================================WT=====================================
if strcmpi(fres,'log')
    
    %Direct reconstruction-------------------------------------------------
    if method==0 || method==1
        for tn=tn1:tn2
            ii=tfsupp(2,tn):tfsupp(3,tn); xn=idt(tn); cs=TFR(ii,xn);
            if ~isfinite(D)
                if xn>idt(tn1) && xn<idt(tn2), cw=angle(TFR(ii,xn+1))-angle(TFR(ii,xn-1)); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs/2;
                elseif xn==1, cw=angle(TFR(ii,xn+1))-angle(cs); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs;
                else cw=angle(cs)-angle(TFR(ii,xn-1)); cw(cw<0)=cw(cw<0)+2*pi; cw=cw*fs; end
                cw=cw/(2*pi);
            end
            
            if isempty(cs(isnan(cs) | ~isfinite(cs)))
                casig=(1/C)*sum(cs*fstep); asig(1,tn)=casig;
                if isfinite(D)
                    ifreq(1,tn)=(1/D)*sum(freq(ii).*cs*fstep)/casig;
                else
                    ifreq(1,tn)=(1/C)*sum(cw.*cs*fstep)/casig;
                end
            end
        end
        ifreq(1,:)=real(ifreq(1,:));
    end
    
    %Ridge reconstruction--------------------------------------------------
    if method==0 || method==2
        rm=2; if method==2, rm=1; end
        for tn=tn1:tn2
            ipeak=tfsupp(1,tn); xn=idt(tn);
            if ipeak>1 && ipeak<NF, cs=TFR([ipeak;ipeak-1;ipeak+1],xn);
            else cs=TFR(ipeak,xn); end
            
            if isfinite(cs(1))
                ifreq(rm,tn)=freq(ipeak);
                if ipeak>1 && ipeak<NF %quadratic interpolation
                    a1=abs(cs(2)); a2=abs(cs(1)); a3=abs(cs(3));
                    p=(1/2)*(a1-a3)/(a1-2*a2+a3);
                    if abs(p)<=1, ifreq(rm,tn)=exp(log(ifreq(rm,tn))+p*fstep); end
                end
                ximax=ompeak*ifreq(rm,tn)/freq(ipeak);
                if nflag==0 %if wavelet FT is known in analytic form
                    cmax=fwt(ximax); if isnan(cmax), cmax=fwt(ximax*(1+10^(-14))); end
                    if isnan(cmax), cmax=fwtmax; end
                else %if wavelet FT is numerically estimated
                    cid1=1+floor((ximax-fxi(1))/wstep); cid2=cid1+1;
                    cid1=min([max([cid1,1]),Lf]); cid2=min([max([cid2,1]),Lf]);
                    if cid1==cid2, cmax=fwt(cid1);
                    else
                        cmax=fwt(cid1)+(fwt(cid2)-fwt(cid1))*(ximax-fxi(cid1))/(fxi(cid2)-fxi(cid1));
                    end
                end
                casig=2*cs(1)/cmax; asig(rm,tn)=casig;
            end
        end
    end
    
end

%Estimate amplitude and phase (faster to do all at once)
iamp=abs(asig); iphi=angle(asig);
%Unwrap phases at the end
for sn=1:size(iphi,1), iphi(sn,:)=unwrap(iphi(sn,:)); end

end

