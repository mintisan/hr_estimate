%================ Reconstruct parameters of the component =================
%============= from its time-frequency support in SWFT or SWT =============
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
% [iamp,iphi,ifreq,Optional:rtfsupp] = ssrectfr(tfsupp,TFR,freq,wopt,Optional:method)
% - returns the component's amplitude [iamp], phase [iphi] and frequency
%   [ifreq] as reconstructed from its extracted time-frequency support
%   [tfsupp] in the signal's SWFT/SWT [TFR] (determines whether TFR is SWFT
%   or SWT based on the spacings between specified frequencies [freq] -
%   linear (SWFT) or logarithmic (SWT)). The optional output [rtfsupp]
%   returns the extracted time-frequency support if the input [tfsupp]
%   specifies 1xL frequency profile instead of the full time-frequency
%   support (see below); otherwise returns input [rtfsupp]=[tfsupp].
%
% INPUT:
% tfsupp: 3xL matrix
%        - extracted time-frequency support of the component, containing
%          frequencies of the TFR amplitude peaks (ridge points) in the
%          first row, support lower bounds (referred as \omega_-(t)/2/pi
%          in [1]) - in the second row, and the upper bounds (referred as
%          \omega_+(t)/2/pi in [1]) - in the third row. Alternatively, one
%          can specify [tfsupp] as 1xL vector of the desired frequency
%          profile, in which case the program will automatically select
%          time-frequency support around it and the corresponding peaks.
% TFR: NFxL matrix (rows correspond to frequencies, columns - to time)
%        - time-frequency representation (SWFT or SWT), to which [tfsupp]
%          correspond
% freq: NFx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by function sswft.m or sswt.m
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions sswft, sswt; [wopt] contains all the
%          needed information, i.e. name of the TFR, sampling frequency,
%          window/wavelet characteristics etc.
% method:'direct'(default)|'ridge'|'both'
%        - the reconstruction method to use for estimating the component's
%          parameters [iamp], [iphi], [ifreq] (see [1]); if set to 'both',
%          all parameters are returned as 2xL matrices with direct and
%          ridge estimates corresponding to 1st and 2nd rows, respectively.
%
% NOTE: ridge reconstruction of amplitude is inappropriate for the
% syncrosqueezed TFRs, so that the resultant amplitude estimates should
% not be used if reconstructed by ridge method (while 'direct' is fine).
%
%-------------------------------Examples-----------------------------------
%
% /Calculate first SWFT (see sswft function):/
% [SWFT,freq,wopt]=sswft(sig,fs,...(parameters));
% /Extract somehow the component time-frequency support from it, e.g./
% tfsupp=ssecurve(SWFT,freq,wopt);
% /Then component's parameters can be reconstructed from its extracted
% SWFT support [tfsupp] as/
% [iamp,iphi,ifreq]=ssrectfr(tfsupp,SWFT,freq,wopt);
% /For SWT the same procedure is used./
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


function [iamp,iphi,ifreq,varargout] = ssrectfr(tfsupp,TFR,freq,wopt,varargin)

[NF,L]=size(TFR); freq=freq(:); method=1;
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

%If only the frequency profile is specified, extract the full time-frequency support
if min(size(tfsupp))==1
    eind=tfsupp(:)'; tfsupp=zeros(3,L)*NaN;
    for tn=tn1:tn2
        cind=eind(tn); xn=idt(tn); cs=abs(TFR(:,xn));
        
        %Assure that the point is in the support
        cpeak=cind;
        if cs(cpeak)==0
            cpeak1=cind-1+find(cs(cind:end)>0,1,'first');
            cpeak2=cind+1-find(cs(cind:-1:1)>0,1,'first');
            if ~isempty(cpeak1) && ~isempty(cpeak2)
                if cpeak1-cind==cind-cpeak2
                    if cs(cpeak1)>cs(cpeak2), cpeak=cpeak1;
                    else cpeak=cpeak2; end
                elseif cpeak1-cind<cind-cpeak2, cpeak=cpeak1;
                elseif cpeak1-cind>cind-cpeak2, cpeak=cpeak2;
                end
            elseif isempty(cpeak1), cpeak=cpeak2;
            elseif isempty(cpeak2), cpeak=cpeak1;
            end
        end
        
        %Boundaries of time-frequency support
        iup=[]; idown=[];
        if cpeak<NF-1, iup=cpeak+find(cs(cpeak+1:end)==0,1,'first'); end
        if cpeak>2, idown=cpeak-find(cs(cpeak-1:-1:1)==0,1,'first'); end
        iup=min([iup,NF]); idown=max([idown,1]);
        tfsupp(2,tn)=idown; tfsupp(3,tn)=iup;
        
        %Find the maximum peak in the current support
        ii=idown:iup; [~,idm]=max(abs(cs(ii))); tfsupp(1,tn)=ii(idm);
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

%==================================SWFT====================================
if strcmpi(fres,'linear')
    
    %Direct reconstruction-------------------------------------------------
    if method==0 || method==1
        for tn=tn1:tn2
            ii=tfsupp(2,tn):tfsupp(3,tn); xn=idt(tn); cs=TFR(ii,xn);
            if isempty(cs(isnan(cs) | ~isfinite(cs)))
                casig=sum(cs); asig(1,tn)=casig; cfsig=sum(freq(ii).*cs)/casig;
                if cfsig>=freq(ii(1)) && cfsig<=freq(ii(end)), ifreq(1,tn)=cfsig;
                else ifreq(1,tn)=sum(freq(ii).*abs(cs))/sum(abs(cs)); end
            end
        end
        ifreq(1,:)=real(ifreq(1,:));
    end
    
    %Ridge reconstruction--------------------------------------------------
    if method==0 || method==2
        rm=2; if method==2, rm=1; end
        for tn=tn1:tn2
            ipeak=tfsupp(1,tn); xn=idt(tn); cs=TFR(ipeak,xn);
            if isfinite(cs(1))
                asig(rm,tn)=cs(1); ifreq(rm,tn)=freq(ipeak);
            end
        end
    end
    
end

%===================================SWT====================================
if strcmpi(fres,'log')
    
    %Direct reconstruction-------------------------------------------------
    if method==0 || method==1
        for tn=tn1:tn2
            ii=tfsupp(2,tn):tfsupp(3,tn); xn=idt(tn); cs=TFR(ii,xn);
            if isempty(cs(isnan(cs) | ~isfinite(cs)))
                casig=sum(cs); asig(1,tn)=casig; cfsig=sum(freq(ii).*cs)/casig;
                if cfsig>=freq(ii(1)) && cfsig<=freq(ii(end)), ifreq(1,tn)=cfsig;
                else ifreq(1,tn)=sum(freq(ii).*abs(cs))/sum(abs(cs)); end
            end
        end
        ifreq(1,:)=real(ifreq(1,:));
    end
    
    %Ridge reconstruction--------------------------------------------------
    if method==0 || method==2
        rm=2; if method==2, rm=1; end
        for tn=tn1:tn2
            ipeak=tfsupp(1,tn); xn=idt(tn); cs=TFR(ipeak,xn);
            if isfinite(cs(1))
                asig(rm,tn)=cs(1); ifreq(rm,tn)=freq(ipeak);
            end
        end
    end
    
end

%Estimate amplitude and phase (faster to do all at once)
iamp=abs(asig); iphi=angle(asig);
%Unwrap phases at the end
for sn=1:size(iphi,1), iphi(sn,:)=unwrap(iphi(sn,:)); end

end

