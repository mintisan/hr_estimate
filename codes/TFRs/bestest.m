%========== Determine the best among direct and ridge estimates ===========
% Version 1.00 stable
%------------------------------Requirements--------------------------------
% requires wft.m, wt.m, rectfr.m
%-------------------------------Copyright----------------------------------
%
% Author: Dmytro Iatsenko
% Information about these codes (e.g. links to the Video Instructions),
% as well as other MatLab programs and many more can be found at
% http://www.physics.lancs.ac.uk/research/nbmphysics/diats/tfr
% 
% Related articles:
% [1] D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part I: Overview, standards of use, related issues and algorithms."
% {preprint - arXiv:1310.7215}
% [2] D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part II: Resolution, reconstruction and concentration."
% {preprint - arXiv:1310.7274}
% [3] D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "On the extraction of instantaneous frequencies from ridges in
%  time-frequency representations of signals."
% {preprint - arXiv:1310.7276}
% [4] D. Iatsenko, D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "Nonlinear Mode Decomposition: a new noise-robust, adaptive
%  decomposition method."
% {preprint - arXiv:1207.5567}
%
%---------------------------- Documentation -------------------------------
% [eamp,ephi,efreq,Optional:allerr,allrec]=bestest(tfsupp,TFR,freq,wopt,Optional:'PropertyName',PropertyValue);
% - given the extracted component's support [tfsupp] in the current [TFR],
%   reconstructs its instantaneous parameters (amplitude, phase and
%   frequency) by both direct and ridge methods, between which it then
%   chooses the best estimates based on the criteria proposed in [2];
%   these 'best' estimates are then returned in [eamp] (amplitude), [ephi]
%   (phase) and [efreq] (frequency).
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
%          corresponds
% freq: NFx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by function wft.m or wt.m
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions wft.m, wt.m; [wopt] contains all the
%          needed information, i.e. name of the TFR, sampling frequency,
%          window/wavelet characteristics etc.
% 
% PROPERTIES: ({{...}} denotes default)
% ################################ BASIC ##################################
% 'Display':{{'on'}}|'on-'|'off'
%        - to display or not the progress information; if set to 'on-',
%          then displays the information but more compactly.
% ############################## ADVANCED #################################
% 'Kappa':2x3 matrix (default = [[3;4;2],[1;1;1]])
%        - the multipliers before the errors for each estimate (columns -
%          direct and ridge, rows - amplitude, phase and frequency),
%          required to make these errors for different methods comparable
%          (see [2,4]); the default values were found empirically.
% 'MinSupp':{{'on'}}|'off'|SuppAcc|[SuppAcc,SuppPerc]
%        - use or not the minimal support of the component, which is
%          needed solely to improve the speed of the method. Thus, when
%          the TFR needs to be recalculated, it is recalculated not in
%          the whole frequency range (which would be quite expensive
%          computationally), but only in the range where the component is
%          expected to reside. To estimate this range, the original
%          component's support [tfsupp] is first squeezed at each time to
%          the region where the TFR amplitude is higher than some portion
%          (specified by [SuppAcc], default = 0.001) of the peak amplitude
%          at that time; then the respective frequency range is estimated
%          as some percentage (specified by [SuppPerc], default = 0.95) of
%          the lowest lower and highest upper values of such squeezed
%          support; see [4] for more details. Setting 'MinSupp' to 'off'
%          is equivalent to [SuppAcc]=0, [SuppPerc]=1, while specifying it
%          as a vector allows to manually define [SuppAcc] and [SuppPerc].
% 
% OUTPUT:
% eamp,ephi,efreq: 1xL vectors
%   - chosen estimates ([eamp] - amplitude, [ephi] - phase, [efreq] -
%     frequency), i.e. those having the smallest re-extraction errors
% allerr: 3x2 vectors
%   - the corresponding re-extraction errors (multiplied on [kappa]):
%     amplitude - 1 row, phase - 2 row, frequency - 3 row;
%     direct estimates - 1 column, ridge estimates - 2 column.
%     For example, if allerr(1,1)<=allerr(1,2), then the direct estimate
%     of amplitude is preferred and returned in [eamp]. Note, that
%     allerr(3,:) are the discrepancies for the usual frequency (in Hz),
%     i.e. the ones defined in [2,4] but divided on 2\pi.
% allrec: 3x2 cell of 3xL matrices
%   - all reconstructed parameters, with 3xL matrix in each cell containing
%     the estimated amplitude, phase and frequency in 1, 2 and 3 rows,
%     respectively; allerr{1,1} and allerr{1,2} are the original direct
%     and ridge estimates, respectively; allerr{2,1} and allerr{2,2}
%     are the same estimates but obtained from the time-frequency support
%     'cutted' to the minimal range; and allerr{3,1} and allerr{3,2} are
%     the refined direct and ridge estimates, respectively.
% 
% WARNING: the calculated discrepancies [allerr] can be used ONLY to
% compare the quality of direct and ridge estimates; [allerr] CANNOT be
% used to determine the overal accuracy of the obtained estimates, or to
% judge about anything else.
% 
% NOTE: One can alternatively pass the structure with the properties as the
% 5th argument, e.g. /opt.MinSupp='off'; bestest(tfsupp,TFR,freq,wopt,opt);/.
% If the other properties are specified next, they override those in the
% structure, e.g. /bestest(tfsupp,TFR,freq,wopt,opt,'MinSupp','on');/ will
% always use 'MinSupp'='on', irrespectively to what is specified in [opt].
% 
%-------------------------------Examples-----------------------------------
% 
% /Calculate TFR, WFT in the present case ([fs] denote sampling frequency)/
% [WFT,freq,wopt]=wft(signal,fs);
% /Extract the curve from it/
% tfsupp=ecurve(WFT,freq,wopt);
% /Reconstruct the parameters and return the better estimates/
% [eamp,ephi,efreq]=recerr(tfsupp,WFT,freq,wopt);
% 
%--------------------------------------------------------------------------


function [eamp,ephi,efreq,varargout]=bestest(tfsupp,TFR,freq,wopt,varargin)

%Default parameters
Kappa=[[3;4;2],[1;1;1]]; MinSupp='on'; DispMode='on';
%Update if user-defined
vst=1;
if nargin>4 && isstruct(varargin{1})
    copt=varargin{1}; vst=2;
    if isfield(copt,'Kappa'), cvv=copt.Kappa; if ~isempty(cvv), Kappa=cvv; end, end
    if isfield(copt,'MinSupp'), cvv=copt.MinSupp; if ~isempty(cvv), MinSupp=cvv; end, end
    if isfield(copt,'Display'), cvv=copt.Display; if ~isempty(cvv), DispMode=cvv; end, end
end
for vn=vst:2:nargin-4
    if strcmpi(varargin{vn},'Kappa'), if ~isempty(varargin{vn+1}), Kappa=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'MinSupp'), if ~isempty(varargin{vn+1}), MinSupp=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Display'), if ~isempty(varargin{vn+1}), DispMode=varargin{vn+1}; end
    else error(['There is no Property ''',varargin{vn},'''']);
    end
end
if ~strcmpi(DispMode,'off'), fprintf('Selecting the optimal method for reconstructing parameters of the component:\n'); end

%Determine the frequency resolution
if min(freq)<=0 || std(diff(freq))<std(diff(log(freq))), fres=1;
else fres=2; end

%Reconstruct the parameters from a full time-frequency support
[eamp,ephi,efreq]=rectfr(tfsupp,TFR,freq,wopt,'both');

%Determine the frequency range and reconstruct the parameters from it
[~,mrange]=minsupp(tfsupp,TFR,freq,wopt,MinSupp);
tfsupp(2,:)=max(tfsupp(2,:),mrange(1));
tfsupp(3,:)=min(tfsupp(3,:),mrange(2));
[iamp,iphi,ifreq]=rectfr(tfsupp,TFR,freq,wopt,'both');
wopt.fmin=mrange(1); wopt.fmax=mrange(2);

sig=iamp.*cos(iphi); ramp=iamp; rphi=iphi; rfreq=ifreq;
%Direct refinement
if strcmpi(DispMode,'on'), fprintf('-- calculating discrepancies of the direct estimates;\n'); end
if fres==1, [TFR,freq]=wft(sig(1,:),wopt.fs,wopt,'Display','off','Plot','off');
else [TFR,freq]=wt(sig(1,:),wopt.fs,wopt,'Display','off','Plot','off'); end
tfsupp=ecurve_max(TFR,freq); [camp,cphi,cfreq]=rectfr(tfsupp,TFR,freq,wopt,'both');
ramp(1,:)=camp(1,:); rphi(1,:)=cphi(1,:); rfreq(1,:)=cfreq(1,:);
%Ridge refinement
if strcmpi(DispMode,'on'), fprintf('-- calculating discrepancies of the ridge estimates;\n'); end
if fres==1, [TFR,freq]=wft(sig(2,:),wopt.fs,wopt,'Display','off','Plot','off');
else [TFR,freq]=wt(sig(2,:),wopt.fs,wopt,'Display','off','Plot','off'); end
tfsupp=ecurve_max(TFR,freq); [camp,cphi,cfreq]=rectfr(tfsupp,TFR,freq,wopt,'both');
ramp(2,:)=camp(2,:); rphi(2,:)=cphi(2,:); rfreq(2,:)=cfreq(2,:);

[aerr,perr,ferr]=errcalc(iamp,iphi,ifreq,ramp,rphi,rfreq);
allerr=[aerr(:)';perr(:)';ferr(:)']; allerr=Kappa.*allerr;
allrec={[eamp(1,:);ephi(1,:);efreq(1,:)],[eamp(2,:);ephi(2,:);efreq(2,:)];...
    [iamp(1,:);iphi(1,:);ifreq(1,:)],[iamp(2,:);iphi(2,:);ifreq(2,:)];...
    [ramp(1,:);rphi(1,:);rfreq(1,:)],[ramp(2,:);rphi(2,:);rfreq(2,:)]};

aest='direct'; pest='direct'; fest='direct';
if allerr(1,1)<=allerr(1,2), eamp=eamp(1,:); else eamp=eamp(2,:); aest='ridge'; end
if allerr(2,1)<=allerr(2,2), ephi=ephi(1,:); else ephi=ephi(2,:); pest='ridge'; end
if allerr(3,1)<=allerr(3,2), efreq=efreq(1,:); else efreq=efreq(2,:); fest='ridge'; end
if strcmpi(DispMode,'on'), fprintf(['-- chosen estimates: amplitude - ',aest,', phase - ',pest,', frequency - ',fest,'.\n']); end
if strcmpi(DispMode,'on-'), fprintf(['\b amplitude - ',aest,', phase - ',pest,', frequency - ',fest,'.\n']); end

if nargout>3, varargout{1}=allerr; end
if nargout>4, varargout{2}=allrec; end

end


%------------- Function for maximum-based curve extraction ----------------
function tfsupp=ecurve_max(TFR,freq)

[NF,L]=size(TFR); TFR=abs(TFR);
tn1=find(~isnan(TFR(end,:)),1,'first'); tn2=find(~isnan(TFR(end,:)),1,'last');

ridges=zeros(L,1)*NaN; maxv=zeros(L,1)*NaN;
for tn=tn1:tn2, [maxv(tn),ridges(tn)]=max(abs(TFR(:,tn))); end
idz=tn1-1+find(maxv(tn1:tn2)==0);
if ~isempty(idz)
    idnz=tn1:tn2; idnz=idnz(~ismember(idnz,idz));
    ridges(idz)=interp1(idnz,ridges(idnz),idz,'linear','extrap');
    ridges(idz)=round(ridges(idz));
end
tfsupp(1,:)=ridges(:)';

for tn=tn1:tn2
    cs=abs(TFR(:,tn)); cpeak=tfsupp(1,tn);
    iup=[]; idown=[];
    if cpeak<NF-1, iup=cpeak+find(cs(cpeak+1:end-1)<=cs(cpeak:end-2) & cs(cpeak+1:end-1)<cs(cpeak+2:end),1,'first'); end
    if cpeak>2, idown=cpeak-find(cs(cpeak-1:-1:2)<=cs(cpeak:-1:3) & cs(cpeak-1:-1:2)<cs(cpeak-2:-1:1),1,'first'); end
    iup=min([iup,NF]); idown=max([idown,1]);
    tfsupp(2,tn)=idown; tfsupp(3,tn)=iup;
end

tfsupp(:,tn1:tn2)=freq(tfsupp(:,tn1:tn2));

end

%----------------- Function for calculating the errors --------------------
function [aerr,perr,ferr]=errcalc(iamp,iphi,ifreq,ramp,rphi,rfreq)

NE1=size(iamp,1); NE2=size(ramp,1);
if NE1==NE2
    NE=NE1; aerr=zeros(NE,1); perr=zeros(NE,1); ferr=zeros(NE,1);
    for en=1:NE
        aerr(en)=sqrt(mean((ramp(en,:)-iamp(en,:)).^2));
        perr(en)=sqrt(1-abs(mean(exp(1i*(rphi(en,:)-iphi(en,:)))))^2);
        ferr(en)=sqrt(mean((rfreq(en,:)-ifreq(en,:)).^2));
    end
else
    aerr=zeros(NE1,NE2); perr=zeros(NE1,NE2); ferr=zeros(NE1,NE2);
    for en1=1:NE1
        for en2=1:NE2
            aerr(en1,en2)=sqrt(mean((ramp(en2,:)-iamp(en1,:)).^2));
            perr(en1,en2)=sqrt(1-abs(mean(exp(1i*(rphi(en2,:)-iphi(en1,:)))))^2);
            ferr(en1,en2)=sqrt(mean((rfreq(en2,:)-ifreq(en1,:)).^2));
        end
    end
end

end

%--------------- Function for finding the minimal support -----------------
function [msupp,varargout] = minsupp(tfsupp,TFR,freq,wopt,MinSupp)

msupp=tfsupp; if nargout>1, varargout{1}=[min(tfsupp(2,:));max(tfsupp(3,:))]; end
SuppAcc=0.001; SuppPerc=0.95;
if strcmpi(MinSupp,'off'), return; end
if ~ischar(MinSupp)
    SuppAcc=MinSupp(1); if length(MinSupp)>1, SuppPerc=MinSupp(2); end
end

tn1=find(~isnan(tfsupp(1,:)),1,'first'); tn2=find(~isnan(tfsupp(1,:)),1,'last');
if min(freq)<=0 || std(diff(freq))<std(diff(log(freq))), fres=1; else fres=2; end

if fres==1, idsupp=1+floor((1/2)+(tfsupp-freq(1))/mean(diff(freq)));
else idsupp=1+floor((1/2)+(log(tfsupp)-log(freq(1)))/mean(diff(log(freq)))); end

for tn=tn1:tn2
    ii=idsupp(2,tn):idsupp(3,tn); cs=abs(TFR(ii,tn)); [vmax,imax]=max(cs);
    idown=imax+1-find(cs(imax:-1:1)>=SuppAcc*vmax,1,'last'); idown=max([1,idown]);
    iup=imax-1+find(cs(imax:end)>=SuppAcc*vmax,1,'last'); iup=min([iup,length(ii)]);
    msupp(2,tn)=freq(ii(idown)); msupp(3,tn)=freq(ii(iup));
end

if nargout>1
    mrange=[sort(msupp(2,:),'descend');sort(msupp(3,:))];
    mrange=mrange(:,round(SuppPerc*(tn2-tn1+1)));
    varargout{1}=[min([mrange(1),min(tfsupp(1,:))]);max([mrange(2),max(tfsupp(1,:))])];
end

end

