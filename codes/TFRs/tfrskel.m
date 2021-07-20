% Version 1.00 stable
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
% RTFR=tfrskel(TFR,freq,wopt,Optional:method)
% constructs the skeleton [RTFR] of the original [TFR], i.e. partitions the
% latter on time-frequency supports and reconstructs each component's
% amplitude, phase and frequency, assigning the corresponding values to a
% particular frequency bins; determines whether TFR is WFT or WT based on
% the spacings between specified frequencies [freq] - linear (WFT) or
% logarithmic (WT).
%
% INPUT: ({{...}} denotes default)
% TFR: FNxL matrix (rows correspond to frequencies, columns - to time)
%        - WFT or WT for which to calculate skeleton
% freq: FNx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by wft,wt,sswft or sswt
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions wft, wt, sswft and sswt; [wopt]
%          contains all the needed information, i.e. name of the TFR,
%          sampling frequency, window/wavelet characteristics etc.
% method:{{'ridge'}}|'direct'
%        - method for amplitude, phase and frequency reconstruction
%
% NOTE: if the window/wavelet does not allow direct estimation of frequency
%       from WFT/WT, i.e. [wopt.omg=Inf] for WFT or [wopt.D=Inf] for WT (as 
%       for the Morlet wavelet), corresponding to not finite \bar{\omega}_g
%       and D_\psi in [1], then uses hybrid frequency reconstruction
%       instead of the direct one, for which calculates the instantaneous
%       TFR frequency (referred as \nu_{G,W} in [1]).
%
%-------------------------------Examples-----------------------------------
%
% /Calculate first WFT (see wft function):/
% [WFT,freq,wopt]=wft(sig,fs,...(parameters));
% /Then to obtain ridge-based WFT skeleton, use/
% WFTrskel=tfrskel(WFT,freq,wopt,'ridge');
% /while to obtain direct-based WFT skeleton, use/
% WFTdskel=tfrskel(WFT,freq,wopt,'direct');
% /The same applies for WT (just replace everywhere wft->wt, WFT->WT)./
%
%--------------------------------------------------------------------------


function  RTFR = tfrskel(TFR,freq,wopt,varargin)

[FN,L]=size(TFR); freq=freq(:); fs=wopt.fs; wp=wopt.wp;
method='ridge'; if nargin>3, method=varargin{1}; end
RTFR=zeros(FN,L); RTFR(isnan(TFR))=NaN;

%Determine the frequency resolution and transform [tfsupp] from frequencies to indices
if freq(1)<=0 || std(diff(freq))<std(freq(2:end)./freq(1:end-1))
    fres='linear'; fstep=mean(diff(freq));
else
    fres='log'; fstep=mean(diff(log(freq)));
end

%Calculate TFR frequency (only if one is forced to use hybrid frequency estimation)
IFR=[];
if (isfield(wp,'omg') && wp.omg==Inf) || (isfield(wp,'D') && wp.D==Inf)
    if strcmpi(method,'direct')
        IFR=angle(TFR); IFR=diff(IFR,1,2);
        unwr=find(IFR>pi); IFR(unwr)=IFR(unwr)-2*pi; %unwrapping
        unwr=find(IFR<-pi); IFR(unwr)=IFR(unwr)+2*pi; %unwrapping
        IFR=horzcat(IFR(:,1)*fs/(2*pi),(IFR(:,1:end-1)+IFR(:,2:end))*fs/(2*2*pi),IFR(:,end)*fs/(2*pi));
    end
end

%For JIT acceleration, define variables to not use cells or structures
C=wp.C; ompeak=wp.ompeak;
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

if strcmpi(fres,'linear') %for WFT
    tn1=find(~isnan(TFR(end,:)),1,'first'); tn2=find(~isnan(TFR(end,:)),1,'last');
    for tn=tn1:tn2
        cs=TFR(:,tn); acs=abs(cs);
        
        if strcmpi(method,'ridge')
            idr=1+find(acs(2:end-1)>=acs(1:end-2) & acs(2:end-1)>acs(3:end));
            for rn=1:length(idr)
                ipeak=idr(rn);
                ifreq=freq(ipeak)-ompeak/2/pi;
                if ipeak>1 && ipeak<FN %quadratic interpolation
                    a1=acs(ipeak-1); a2=acs(ipeak); a3=acs(ipeak+1);
                    p=(1/2)*(a1-a3)/(a1-2*a2+a3);
                    if abs(p)<=1, ifreq=ifreq+p*fstep; end
                end
                ximax=2*pi*(freq(ipeak)-ifreq);
                if nflag==0 %if window FT is known in analytic form
                    cmax=fwt(ximax);
                else %if window FT is numerically estimated
                    cid1=1+floor((ximax-fxi(1))/wstep); cid2=cid1+1;
                    cid1=min([max([cid1,1]),Lf]); cid2=min([max([cid2,1]),Lf]);
                    if cid1==cid2, cmax=fwt(cid1);
                    else
                        cmax=fwt(cid1)+(fwt(cid2)-fwt(cid1))*(ximax-fxi(cid1))/(fxi(cid2)-fxi(cid1));
                    end
                end
                if ifreq>freq(1)-fstep/2 && ifreq<freq(end)+fstep/2
                    RTFR(1+floor((1/2)+(ifreq-freq(1))/fstep),tn)=2*cs(ipeak)/cmax;
                end
            end
        else
            idm=[find(~isnan(acs),1,'first');1+find(acs(2:end-1)<=acs(1:end-2) & acs(2:end-1)<acs(3:end));find(~isnan(acs),1,'last')];
            for sn=1:(length(idm)-1)
                ii=idm(sn):idm(sn+1);
                ansig=(1/C)*sum(cs(ii)*(2*pi*fstep));
                if isfinite(omg)
                    ifreq=real(-omg+(1/C)*sum(freq(ii).*cs(ii)*(2*pi*fstep))/ansig);
                else
                    ifreq=real((1/C)*sum(IFR(ii,tn).*cs(ii)*(2*pi*fstep))/ansig);
                end
                if ifreq>freq(1)-fstep/2 && ifreq<freq(end)+fstep/2
                    RTFR(1+floor((1/2)+(ifreq-freq(1))/fstep),tn)=ansig;
                end
            end
        end
    end
    
else %for WT
    tn1=find(~isnan(TFR(end,:)),1,'first'); tn2=find(~isnan(TFR(end,:)),1,'last');
    for tn=tn1:tn2
        cs=TFR(:,tn); acs=abs(cs);
        
        if strcmpi(method,'ridge')
            idr=1+find(acs(2:end-1)>=acs(1:end-2) & acs(2:end-1)>acs(3:end));
            for rn=1:length(idr)
                ipeak=idr(rn);
                ifreq=freq(ipeak);
                if ipeak>1 && ipeak<FN %quadratic interpolation
                    a1=acs(ipeak-1); a2=acs(ipeak); a3=acs(ipeak+1);
                    p=(1/2)*(a1-a3)/(a1-2*a2+a3);
                    if abs(p)<=1, ifreq=exp(log(ifreq)+p*fstep); end
                end
                ximax=ompeak*ifreq/freq(ipeak);
                if nflag==0 %if wavelet FT is known in analytic form
                    cmax=fwt(ximax);
                else %if wavelet FT is numerically estimated
                    cid1=1+floor((ximax-fxi(1))/wstep); cid2=cid1+1;
                    cid1=min([max([cid1,1]),Lf]); cid2=min([max([cid2,1]),Lf]);
                    if cid1==cid2, cmax=fwt(cid1);
                    else
                        cmax=fwt(cid1)+(fwt(cid2)-fwt(cid1))*(ximax-fxi(cid1))/(fxi(cid2)-fxi(cid1));
                    end
                end
                if log(ifreq)>log(freq(1))-fstep/2 && log(ifreq)<log(freq(end))+fstep/2
                    RTFR(1+floor((1/2)+log(ifreq/freq(1))/fstep),tn)=2*cs(ipeak)/cmax;
                end
            end
        else
            idm=[find(~isnan(acs),1,'first');1+find(acs(2:end-1)<=acs(1:end-2) & acs(2:end-1)<acs(3:end));find(~isnan(acs),1,'last')];
            for sn=1:(length(idm)-1)
                ii=idm(sn):idm(sn+1);
                ansig=(1/C)*sum(cs(ii)*fstep);
                if isfinite(D)
                    ifreq=real((1/D)*sum(freq(ii).*cs(ii)*fstep)/ansig);
                else
                    ifreq=real((1/C)*sum(IFR(ii,tn).*cs(ii)*fstep)/ansig);
                end
                if log(ifreq)>log(freq(1))-fstep/2 && log(ifreq)<log(freq(end))+fstep/2
                    RTFR(1+floor((1/2)+log(ifreq/freq(1))/fstep),tn)=ansig;
                end
            end
        end
    end


end

