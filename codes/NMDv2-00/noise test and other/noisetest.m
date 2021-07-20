%============ Test the extracted component for being not noise ============
% Version 1.00 stable
%------------------------------Requirements--------------------------------
% requires wft.m, wt.m, ecurve.m, rectfr.m
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
% [3] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "On the extraction of instantaneous frequencies from ridges in
%  time-frequency representations of signals."
% {preprint:arXiv:1310.7276}
% [4] D. Iatsenko, D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "Nonlinear Mode Decomposition: a new noise-robust, adaptive
%  decomposition method."
% {preprint - arXiv:1207.5567}
% 
%------------------------------Documentation-------------------------------
% [p,Optional:signif,DA,DF]=noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo,Optional:'PropertyName',PropertyValue)
% - tests whether the extracted component (given by its time-frequency
%   support [tfsupp] in the WFT or WT [TFR]) is genuine (p=1) or is just
%   formed from the noise pciked in the time-frequency plane (p=0).
% 
% INPUT:
% sig: 1xL OR Lx1 vector
%        - the original signal from which the component was extracted
% TFR: NFxL matrix (rows correspond to frequencies, columns - to time)
%        - the signal's time-frequency representation (WFT or WT)
% freq: NFx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by function wft.m or wt.m
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions wft.m, wt.m; [wopt] contains all the
%          needed information, i.e. name of the TFR, sampling frequency,
%          window/wavelet characteristics etc.
% tfsupp: 3xL matrix  (or 1xL vector of the desired frequency profile)
%        - extracted time-frequency support of the component, containing
%          frequencies of the TFR amplitude peaks (ridge points) in the
%          first row, support lower bounds (referred as \omega_-(t)/2/pi
%          in [1]) - in the second row, and the upper bounds (referred as
%          \omega_+(t)/2/pi in [1]) - in the third row. Alternatively, one
%          can specify [tfsupp] as 1xL vector of the desired frequency
%          profile, in which case the program will automatically select
%          time-frequency support around it and the corresponding peaks.
% ecinfo: structure returned by function ecurve.m
%        - parameters of the procedure used for component extraction,
%          returned as a second output by function ecurve.m
% 
% OUTPUT:
% p: 0 or 1
%        - the outcome of the test: p=0 means that the test accepted the
%          null hypothesis of noise, i.e. the extracted component might
%          have no physical meaning, while p=1 indicates that the null
%          hypothesis was rejected, i.e. the extracted component is
%          most probably physically meaningful.
% signif: NAx1 vector
%        - significance of the surrogate test for each of the [NA]
%          discrimination statistics that were used (see 'Alpha' property);
%          the overall significance of the test is taken as maximum among
%          all values of [signif].
% DA: (1+[SurrNum])x1 vector
%        - original (first) and surrogate ([SurrNum] other) values of the
%          component's amplitude spectral entropy (see [4]).
% DF: (1+[SurrNum])x1 vector
%        - original (first) and surrogate ([SurrNum] other) values of the
%          component's frequency spectral entropy (see [4]).
% 
% PROPERTIES: ({{...}} denotes default)
% ################################ BASIC ##################################
% 'Display':{{'on'}}|'on-'|'off'
%        - to display or not the progress information; if set to 'on-',
%          then displays the information but more compactly.
% ############################## ADVANCED #################################
% 'Alpha':NAx2 matrix (default = [1,0;0,1;1,1])
%        - multipliers that determine the discriminating statistics in use;
%          thus, the [NA] discriminating statisctics are calculated as
%          Alpha(:,1)*DA+Alpha(:,2)*DF, where [DA] and [DF] are amplitude
%          and frequency spectral entropies, respectively; the overall
%          significance of the surrogate test is then taken as the maximum
%          among significances for each of the [NA] measures.
% 'Signif':value from 0 to 1 (default = 0.95)
%        - significance level for rejecting the null hypothesis of noise,
%          so that the component is regarded as genuine if the overall
%          significance of the test exceeds this level.
% 'SurrNum':positive integer value (default = 40)
%        - number of surrogates to use for testing
% 'SurrType':{{'FT'}}|'AAFT'|'IAAFT'
%        - type of the surrogates to use for the test: 'FT' - Fourier
%          Transform surrogates, obtained by randomizing the phases of
%          the original Fourier coefficients, preserving the power
%          spectrum (addresses the null hypothesis of linear noise);
%          'AAFT' - amplitude adjusted FT surrogates, which preserve the
%          distribution of signal values and also try to preserve the
%          power spectrum (addresses the null hypothesis of invertibly
%          rescaled linear noise, but does not always do this well due to
%          not exactly preserving the power spectrum); 'IAAFT' - iterative
%          AAFT, which preserve both the power spectrum and, as much as
%          possible, the distribution of values (addresses the null
%          hypothesis of rescaled linear noise).
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
% NOTE: One can alternatively pass the structure with the properties as the
% 5th argument, e.g. /opt.SurrNum=20; noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo,opt);/.
% If the other properties are specified next, they override those in the
% structure, e.g. /noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo,opt,'SurrNum',100);/
% will always use 100 surrogates for testing, irrespectively to what is
% specified in [opt].
%
%-------------------------------Examples-----------------------------------
% 
% /Calculate TFR, WFT in the present case ([fs] denote sampling frequency)/
% [TFR,freq,wopt]=wft(signal,fs);
% /Extract the component's time-frequency support from the signal's TFR/
% [tfsupp,ecinfo]=ecurve(TFR,freq,wopt);
% /Test whether the extracted component have physical meaning/
% p=noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo);
% 
%--------------------------------------------------------------------------


function [p,varargout]=noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo,varargin)

%Default parameters
DispMode='on'; Alpha=[1,0;0,1;1,1]; SLevel=0.95; NS=40; SurrType='FT'; MinSupp='on';
%Update if user-defined
vst=1;
if nargin>6 && isstruct(varargin{1})
    copt=varargin{1}; vst=2;
    if isfield(copt,'Display'), cvv=copt.Display; if ~isempty(cvv), DispMode=cvv; end, end
    if isfield(copt,'Alpha'), cvv=copt.Alpha; if ~isempty(cvv), Alpha=cvv; end, end
    if isfield(copt,'Signif'), cvv=copt.Signif; if ~isempty(cvv), SLevel=cvv; end, end
    if isfield(copt,'SurrNum'), cvv=copt.SurrNum; if ~isempty(cvv), NS=cvv; end, end
    if isfield(copt,'SurrType'), cvv=copt.SurrType; if ~isempty(cvv), SurrType=cvv; end, end
    if isfield(copt,'MinSupp'), cvv=copt.MinSupp; if ~isempty(cvv), MinSupp=cvv; end, end
end
for vn=vst:2:nargin-6
    if strcmpi(varargin{vn},'Display'), if ~isempty(varargin{vn+1}), DispMode=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Alpha'), if ~isempty(varargin{vn+1}), Alpha=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Signif'), if ~isempty(varargin{vn+1}), SLevel=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'SurrNum'), if ~isempty(varargin{vn+1}), NS=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'SurrType'), if ~isempty(varargin{vn+1}), SurrType=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'MinSupp'), if ~isempty(varargin{vn+1}), MinSupp=varargin{vn+1}; end
    else error(['There is no Property ''',varargin{vn},'''']);
    end
end
if ~strcmpi(DispMode,'off'), fprintf('Testing the component against noise: '); end
if strcmpi(DispMode,'on-'), cstr='0%%'; fprintf(cstr); cslen=length(cstr); end

%Determine the frequency resolution and detrend the signal
if min(freq)<=0 || std(diff(freq))<std(diff(log(freq))), fres=1; else fres=2; end
X=(1:length(sig))'/wopt.fs; XM=ones(length(X),4); for pn=1:3, CX=X.^pn; XM(:,pn+1)=(CX-mean(CX))/std(CX); end
w=warning('off','all'); sig=sig(:)-XM*(pinv(XM)*sig(:)); warning(w);

%Determine the minimal frequency range and re-extract the component
[~,mrange]=minsupp(tfsupp,TFR,freq,wopt,MinSupp); wopt.fmin=mrange(1); wopt.fmax=mrange(2);
if strcmpi(DispMode,'on'), fprintf('\n-- re-extracting component from the minimal range [%0.3f,%0.3f] Hz',mrange(1),mrange(2)); end
if fres==1, [TFR,freq,wopt]=wft(sig,wopt.fs,wopt,'Padding',0,'Display','off','Plot','off');
else [TFR,freq,wopt]=wt(sig,wopt.fs,wopt,'Padding',0,'Display','off','Plot','off'); end
[tfsupp,ecinfo]=ecurve(TFR,freq,wopt,ecinfo,'Display','off','Plot','off','Skel',[]);
[iamp,iphi,ifreq]=rectfr(tfsupp,TFR,freq,wopt,'ridge');

%Construct the surrogates
if strcmpi(DispMode,'on'), fprintf(';\n-- creating surrogates: '); end
surr=createsurr(sig,NS,SurrType,DispMode);

%Calculate the original and surrogate spectral entropies
if strcmpi(DispMode,'on'), cstr=['0 (out of ',num2str(NS),')']; fprintf([';\n-- calculating surrogate DS: ',cstr]); cslen=length(cstr); end
if strcmpi(DispMode,'on-'), cstr=[num2str(80/(NS+1)+20,'%0.1f'),'%%']; cslen=length(cstr)-1; end
DA=zeros(1+NS,1); DA(1)=spentropy(iamp);
DF=zeros(1+NS,1); DF(1)=spentropy(ifreq);
for sn=1:NS
    if fres==1, [TFR,freq,wopt]=wft(surr(sn,:),wopt.fs,wopt);
    else [TFR,freq,wopt]=wt(surr(sn,:),wopt.fs,wopt); end
    tfsupp=ecurve(TFR,freq,wopt,ecinfo);
    [iamp,iphi,ifreq]=rectfr(tfsupp,TFR,freq,wopt,'ridge');
    DA(sn+1)=spentropy(iamp); DF(sn+1)=spentropy(ifreq);
    if strcmpi(DispMode,'on'), cstr=[num2str(sn),' (out of ',num2str(NS),')']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr); end
    if strcmpi(DispMode,'on-'), cstr=[num2str(80*(sn+1)/(NS+1)+20,'%0.1f'),'%%']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr)-1; end
end

%Calculate the significances
NA=length(Alpha); signif=zeros(NA,1);
for an=1:NA
    cval=Alpha(an,1)*DA+Alpha(an,2)*DF;
    signif(an)=length(find(cval(2:end)>cval(1)))/NS;
end
p=0; if max(signif)>=SLevel, p=1; end

%Final display
if strcmpi(DispMode,'on')
    cstr=[]; for an=1:NA, cstr=[cstr,num2str(signif(an),'%0.3f'),',']; end
    if p==0, fprintf([';\n-- component was identified as noise (significance=',num2str(max(signif),'%0.3f'),'=max[',cstr,'\b]<',num2str(SLevel,'%0.3f'),').\n']);
    else fprintf([';\n-- component was identified as genuine (significance=',num2str(max(signif),'%0.3f'),'=max[',cstr,'\b]>',num2str(SLevel,'%0.3f'),').\n']); end
end
if strcmpi(DispMode,'on-')
    if p==0, fprintf([' (identified as noise, significance=',num2str(max(signif),'%0.3f'),'<',num2str(SLevel,'%0.3f'),').\n']);
    else fprintf([' (identified as genuine, significance=',num2str(max(signif),'%0.3f'),'>',num2str(SLevel,'%0.3f'),').\n']); end
end

%Assign optional output
if nargout>1, varargout{1}=signif; end
if nargout>2, varargout{2}=DA; end
if nargout>3, varargout{3}=DF; end

end


%==========================================================================
%======================== Supporting functions ============================
%==========================================================================

%----------------- Function for creating the surrogates -------------------
function surr = createsurr(sig,NS,SurrType,varargin)

DispMode='off'; if nargin>3 && ~isempty(varargin{1}), DispMode=varargin{1}; end
L=length(sig); surr=zeros(NS,L); sig=sig(:)';

if strcmpi(DispMode,'on'), cstr=['0 (out of ',num2str(NS),')']; fprintf(cstr); cslen=length(cstr); end
if strcmpi(DispMode,'on-'), cstr=[num2str(80/(NS+1),'%0.1f'),'%%']; fprintf(['\b\b',cstr]); cslen=length(cstr)-1; end

if strcmpi(SurrType,'FT')
    ll=ceil((L+1)/2)-1; ftsig=fft(sig,L);
    for sn=1:NS
        surr(sn,1)=ftsig(1); randph=2*pi*rand(1,ll-1);
        surr(sn,2:ll)=ftsig(2:ll).*exp(1i*randph);
        surr(sn,2+L-ll:L)=conj(fliplr(surr(sn,2:ll)));
        
        surr(sn,:)=real(ifft(surr(sn,:),L));
        if strcmpi(DispMode,'on'), cstr=[num2str(sn),' (out of ',num2str(NS),')']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr); end
        if strcmpi(DispMode,'on-'), cstr=[num2str(80/(NS+1)+20*sn/NS,'%0.1f'),'%%']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr)-1; end
    end
end

if strcmpi(SurrType,'AAFT')
    ll=ceil((L+1)/2)-1; [sortsig,sortind]=sort(sig);
    rankind=zeros(1,L); rankind(sortind)=1:L;
    rescrank=zeros(1,L);
    for sn=1:NS
        rgs=sort(randn(1,L));
        rescsig=rgs(rankind);
        
        ftresc=fft(rescsig,L);
        surr(sn,1)=ftresc(1); randph=2*pi*rand(1,ll-1);
        surr(sn,2:ll)=ftresc(2:ll).*exp(1i*randph);
        surr(sn,2+L-ll:L)=conj(fliplr(surr(sn,2:ll)));
        
        surr(sn,:)=real(ifft(surr(sn,:),L));
        
        [~,rescind]=sort(surr(sn,:)); rescrank(rescind)=1:L;
        surr(sn,:)=sortsig(rescrank);
        if strcmpi(DispMode,'on'), cstr=[num2str(sn),' (out of ',num2str(NS),')']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr); end
        if strcmpi(DispMode,'on-'), cstr=[num2str(80/(NS+1)+20*sn/NS,'%0.1f'),'%%']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr)-1; end
    end
end

if strcmpi(SurrType,'IAAFT')
    ftsig=fft(sig,L); [sortsig,sortind]=sort(sig);
    rankind=zeros(1,L); rankind(sortind)=1:L;
    ovitn=zeros(NS,1); maxitn=1000; %maximum number of iterations
    for sn=1:NS
        surr(sn,:)=sig(randperm(L));       
        
        itn=1; iterrank=rankind; olditrank=zeros(1,L);
        while (max(abs(olditrank-iterrank))~=0 && itn<maxitn)
            olditrank=iterrank;
            iterf=real(ifft(abs(ftsig).*exp(1i*angle(fft(surr(sn,:),L))))); %replace Fourier amplitudes (take real() because makes mistakes of order \epsilon)
            [~,iterind]=sort(iterf); iterrank(iterind)=1:L;
            surr(sn,:)=sortsig(iterrank);
            itn=itn+1;
        end
        ovitn(sn)=itn-1;
        surr(sn,:)=iterf;
        if strcmpi(DispMode,'on'), cstr=[num2str(sn),' (out of ',num2str(NS),')']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr); end
        if strcmpi(DispMode,'on-'), cstr=[num2str(80/(NS+1)+20*sn/NS,'%0.1f'),'%%']; fprintf([repmat('\b',1,cslen),cstr]); cslen=length(cstr)-1; end
    end
end

end

%-------------- Function for calculating spectral entropies ---------------
function DS = spentropy(sig)

pp=abs(fft(sig)).^2;
pp=pp(pp>0)/sum(pp);
DS=sum(-pp.*log(pp));

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

