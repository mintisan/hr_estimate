%========================= Harmonics  extraction ==========================
% Version 1.00 stable
%------------------------------Requirements--------------------------------
% wft.m, wt.m, rectfr.m, bestest.m
%-------------------------------Copyright----------------------------------
%
% Author: Dmytro Iatsenko
% Information about these codes (e.g. links to the Video Instructions),
% as well as other MatLab programs and many more can be found at
% http://www.physics.lancs.ac.uk/research/nbmphysics/diats/nmd
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
%
% [hid,hamp,hphi,hfreq,Optional:hinfo]=eharm(sig,tfsupp,TFR,freq,wopt,Optional:'PropertyName',PropertyValue);
% - extracts the harmonics and returns their numbers in [hid], amplitudes
%   in [hamp], phases in [hphi], and instantaneous frequencies in [hfreq];
%   if specified, returns also the structure [hinfo] with all information
%   about the process of harmonics extraction.
%
% INPUT:
% sig: 1xL or Lx1 vector
%        - the signal from which to extract the harmonics
% tfsupp: 3xL matrix
%        - extracted time-frequency support of the component, containing
%          frequencies of the TFR amplitude peaks (ridge points) in the
%          first row, support lower bounds (referred as \omega_-(t)/2/pi
%          in [1]) - in the second row, and the upper bounds (referred as
%          \omega_+(t)/2/pi in [1]) - in the third row.
% TFR: FNxL matrix (rows correspond to frequencies, columns - to time)
%        - time-frequency representation (WFT or WT), to which [tfsupp]
%          corresponds
% freq: FNx1 vector
%        - the frequencies corresponding to the rows of [TFR]
% wopt: structure returned by function [wft.m] or [wt.m]
%        - parameters of the window/wavelet and the simulation, returned as
%          a third output by functions wft, wt, sswft and sswt; [wopt]
%          contains all the needed information, i.e. name of the TFR,
%          sampling frequency, window/wavelet characteristics etc.
%
%
% Properties: ({{...}} denotes default)
% ################################ BASIC ##################################
% 'Display':{{'on'}}|'off'|'plot'|'plot+'|'plot++'
%        - defines level of displayed information; when set to 'plot',
%          shows the information about adapting the resolution parameter
%          for each harmonic by plotting the corresponding amplitude-phase
%          consistency values in dependence on [f0]; for 'plot+' also
%          plots the TFR (resampled to the dimensions of the plot to not
%          consume too much memory), the extracted harmonic parameters and
%          the distributions of surrogate consistencies for an (adapted)
%          optimal [f0]; for 'plot++' shows the same as for 'plot+', but
%          for all tested [f0], and not only for the optimal one.
% 'S':value(default = 3)|1i*(1xNH integer vector)
%        - stopping factor: the program stops extracting harmonics
%          when it encounters [value] consecutive false harmonics;
%          alternatively, one can specify S as an imaginary vector, in
%          which case it will extracts only the specified harmonics
%          (e.g. setting 'S' to 1i*[2,5] tells the program to extract
%          only the 2nd and 5th harmonics); note, however, that as a first
%          row in [hid], [hamp] etc. is always returned the first
%          (fundamental) harmonic, so that e.g. setting 'S' to 1i*[2,5] and
%          to 1i*[1,2,5] will give the same result.
% 'S0':value(default = 3)|1i*(1xNH integer vector)
%        - a very similar property as the previous, but this one applies to
%          'fractional' harmonics 1/h (1/2,1/3,...): to check whether the
%          extracted component (defined by a specified time-frequency
%          support [tfsupp] in [TFR]) really represents the first (and not
%          higher) harmonic of some NM, the program checks its 1/h
%          harmonics until [value] consecutive ones will not be identified
%          as false; if some of 1/h harmonics were identified as true,
%          switches to the lowest of them as a new first harmonic;
%          similarly to previous property, setting 'S0' to imaginary
%          vector will manually set the harmonics to check (e.g. 1i*[2,5]
%          will check 1/2 and 1/5 harmonics only); to not apply any checks
%          and regard the specified component as first harmonic from the
%          beginning, just set 'S0' to 0.
% 'AdaptRange':{{'auto'}}|[val1,val2]|@(h)[r1(h),r2(h)]
%        - the range in which the optimal resolution parameter [f0] should
%          be initially searched for: 'auto' sets it automatically for each
%          harmonic, as described in [4]; [val1,val2] (e.g.\ [0.5,2]) sets
%          it manually for all harmonics, but if [val1] or [val2] is
%          imaginary, then the corresponding limit is estimated as in
%          'auto' and then multiplied on [imag(val1)] or [imag(val2)],
%          respectively (e.g. if [val1,val2]=[0.5,1i*2], then the lower
%          limit is 0.5 for all harmonics, while upper is determined
%          automatically for each harmonic and then multiplied on 2);
%          passing 1x2 vector function @(h)[r1(h),r2(h)] allows to manually
%          set range [val1,val2] (including imaginary, see above) for each
%          harmonic [h], e.g. passing @(h)[1./h,h] will set for harmonic
%          h=2 the range [1/2,2], for h=3 -- [1/3,3], for h=1/2 --
%          sort([2,1/2])=[1/2,2] and so on, where the fractional h=1/n
%          might be used to test initial component for being the first
%          harmonic, see 'S0' property). Note that, even when the range
%          is set manually, if at some end of this range the consistency
%          [rho] grows, then program continues to check [f0] in the
%          corresponding direction until the maximum appears (unless the
%          [erflag] is set to 0, see 'AdaptRes' property).
% ############################## ADVANCED #################################
% 'AdaptRes':{{'on'}}|'off'|NR|[NR,erflag]|[NR,erflag,acc]|[NR,eflag,acc,mit]
%        - to adapt or not the resolution parameter for each harmonic; if
%          'on', then searches for [f0] in WFT/WT which maximizes the
%          consistency [rho] of the current harmonic with the first one;
%          passing the vector of up to 4 values (e.g. [25,1] or
%          [25,1,0.01,10]), allows to specify additional parameters:
%          [NR] - the number of points onto which the initial search range
%          of the optimal [f0] is broken (default = 10), [erflag] -
%          specifies whether to search further for optimal [f0] if at the
%          ends of the search interval consistency [rho] increases
%          (erflag = 1, default) or not (erflag = 0), [acc] - the accuracy
%          with which to determine the precise optimal [f0] by golden
%          section search (default = wopt.f0/100, i.e. one hundredth of the
%          resolution parameter used originally for [TFR]; this iterative
%          search is performed only if significant dependence was found
%          at some of the tested [f0]; to turn it off, set [acc]=Inf);
%          [mit] - the maximum number of iterations for the iterative
%          search discussed above (default = Inf, i.e.\ not limited).
% 'SNumber':positive integer value (default = 100)
%        - number of surrogates for significance testing
% 'SLevel':value from 0 to 1 (default = 0.95)
%        - the minimal significance of the surrogate test to regard
%          harmonic as true.
% 'WCons':[wa,wp,wf] (default = [1,1,0])
%        - the weighting factors of each consistency ([wa] for amplitude,
%          [wp] for phase, [wf] for frequency) in the overall consistency
%          [rho]=([rho_a]^wa)*([rho_p]^wp)*([rho_f]^wf).
% 'CLevel':value from 0 to 1 (default = 0.5^(wa+wp))
%        - the minimal value of the consistency [rho] to regard the
%          harmonic as true
% 'RecMethod':{{'auto'}}|'direct'|'ridge'
%        - what reconstruction method to use for estimating amplitude,
%          phase and frequency of the first and higher harmonics;
%          if 'auto', then determines it automatically as that giving best
%          consistency [rho] with the first harmonic (while for the first
%          harmonic uses function bestest.m to select optimal method).
% 'MinSupp':{{'on'}}|'off'|SuppAcc|[SuppAcc,SuppPerc]
%        - to use or not the minimal support of the component, which is
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
%          support; see [4] for more details. The frequency ranges in
%          which the TFR is calculated for harmonics are all based on the
%          initial frequency range of the component according to the
%          formulas in [4]. Setting 'MinSupp' to 'off' is equivalent to
%          [SuppAcc]=0, [SuppPerc]=1, while specifying it as a vector
%          allows to manually define [SuppAcc] and [SuppPerc].
%          
%
% NOTE: One can alternatively pass the structure with the properties as the
% 6th argument, e.g. /opt.ResNum=25; eharm(sig,tfsupp,TFR,freq,wopt,opt);/.
% If the other properties are specified next, they override those in the
% structure, e.g. /eharm(sig,tfsupp,TFR,freq,wopt,opt,'ResNum',10);/ will
% always use 'ResNum'=10, irrespectively to what is specified in [opt].
%
%-------------------------------Examples-----------------------------------
%
% /Calculate first WFT of the reference signal [rsig] (see wft function):/
% [WFT,freq,wopt]=wft(rsig,fs);
% /Extract somehow the component time-frequency support from it, e.g./
% tfsupp=ecurve(WFT,freq);
% /Use it to extract all associated harmonics and their characteristics:/
% [hid,hamp,hphi,hfreq]=eharm(sig,tfsupp,TFR,freq,wopt);
% /For WT the same procedure is used./
%
%--------------------------------------------------------------------------

function [hid,hamp,hphi,hfreq,varargout] = eharm(sig,tfsupp,TFR,freq,wopt,varargin)

[NF,L]=size(TFR); fs=wopt.fs; sig=sig(:)';
tn1=find(~isnan(tfsupp(1,:)),1,'first'); tn2=find(~isnan(tfsupp(1,:)),1,'last');
if min(freq)<=0 || std(diff(freq))<std(diff(log(freq))), fres=1; else fres=2; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default parameters
S=3; S0=3;
NS=100; slev=0.95; ww=[1,1,0]; rhomin=[];
AdaptRes='on'; NR=10; erflag=1; ResAcc=wopt.f0/100; mit=Inf; AdaptRange='auto';
RecMethod='auto'; MinSupp='on'; DispMode='on';
%Update if user-defined
vst=1;
if nargin>5 && isstruct(varargin{1})
    copt=varargin{1}; vst=2;
    if isfield(copt,'S'), cvv=copt.S; if ~isempty(cvv), S=cvv; end, end
    if isfield(copt,'S0'), cvv=copt.S0; if ~isempty(cvv), S0=cvv; end, end
    if isfield(copt,'SNumber'), cvv=copt.SNumber; if ~isempty(cvv), NS=cvv; end, end
    if isfield(copt,'SLevel'), cvv=copt.SLevel; if ~isempty(cvv), slev=cvv; end, end
    if isfield(copt,'CLevel'), cvv=copt.CLevel; if ~isempty(cvv), rhomin=cvv; end, end
    if isfield(copt,'AdaptRes'), cvv=copt.AdaptRes; if ~isempty(cvv), AdaptRes=cvv; end, end
    if isfield(copt,'AdaptRange'), cvv=copt.AdaptRange; if ~isempty(cvv), AdaptRange=cvv; end, end
    if isfield(copt,'RecMethod'), cvv=copt.RecMethod; if ~isempty(cvv), RecMethod=cvv; end, end
    if isfield(copt,'MinSupp'), cvv=copt.MinSupp; if ~isempty(cvv), MinSupp=cvv; end, end
    if isfield(copt,'Display'), cvv=copt.Display; if ~isempty(cvv), DispMode=cvv; end, end
    if isfield(copt,'WCons'), cvv=copt.WCons; if ~isempty(cvv), ww=cvv; end, end
end
for vn=vst:2:nargin-5
    if strcmpi(varargin{vn},'S'), if ~isempty(varargin{vn+1}), S=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'S0'), if ~isempty(varargin{vn+1}), S0=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'SNumber'), if ~isempty(varargin{vn+1}), NS=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'SLevel'), if ~isempty(varargin{vn+1}), slev=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'CLevel'), if ~isempty(varargin{vn+1}), rhomin=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'AdaptRes'), if ~isempty(varargin{vn+1}), AdaptRes=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'AdaptRange'), if ~isempty(varargin{vn+1}), AdaptRange=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'RecMethod'), if ~isempty(varargin{vn+1}), RecMethod=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'MinSupp'), if ~isempty(varargin{vn+1}), MinSupp=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Display'), if ~isempty(varargin{vn+1}), DispMode=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'WCons'), if ~isempty(varargin{vn+1}), ww=varargin{vn+1}; end
    else error(['There is no Property ''',varargin{vn},'''']);
    end
end
if isempty(rhomin), rhomin=0.5^(ww(1)+ww(2)); end
if ~ischar(AdaptRes),
    NR=AdaptRes(1); ARL=length(AdaptRes);
    if ARL>1, erflag=AdaptRes(2); end
    if ARL>2, ResAcc=AdaptRes(3); if imag(ResAcc)~=0, ResAcc=abs(ResAcc)*wopt.f0; end, end
    if ARL>3, mit=AdaptRes(4); end
elseif strcmpi(AdaptRes,'off')
    NR=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove the trend from the signal
X=(1:length(sig))'/fs; XM=ones(length(X),4); for pn=1:3, CX=X.^pn; XM(:,pn+1)=(CX-mean(CX))/std(CX); end
w=warning('off','all'); trend=XM*(pinv(XM)*sig(:)); sig=sig(:)'-trend(:)'; warning(w);

iflag=0; %return or not the structure [hinfo] with all information
if nargout>3
    iflag=1;
    hinfo=struct('hnum',[],'signif',[],'cons',[],'f0',[],'ac',[],'pc',[],'fc',[],'rm',[],'iamp',[],'iphi',[],'ifreq',[],'tfsupp',[]);
    hinfo.adapt_all=struct('mid',[],'signif',[],'cons',[],'f0',[],'ac',[],'pc',[],'fc',[]);
    hinfo.check_first=hinfo;
    hinfo.sim=struct('S',S,'S0',S0,'SNumber',NS,'SLevel',slev,'CLevel',rhomin,'WCons',ww,'AdaptRes',[NR,erflag,ResAcc],'AdaptRange',AdaptRange,'RecMethod',RecMethod,'MinSupp',MinSupp,'Display',DispMode);
end

%First harmonic (expected)
BDisp=DispMode; if ~isempty(strfind(lower(DispMode),'plot')), BDisp='on'; end
if strcmpi(RecMethod,'auto'), [hamp,hphi,hfreq]=bestest(tfsupp,TFR,freq,wopt,'MinSupp',MinSupp,'Display',BDisp);
else [hamp,hphi,hfreq]=rectfr(tfsupp,TFR,freq,wopt,RecMethod); end

%Find the minimal support around the first harmonic curve (to know frequency range for harmonics)
[~,brange]=minsupp(tfsupp,TFR,freq,wopt,MinSupp); brange=[brange(:);brange(:)];
sridgef=sort(tfsupp(1,:)); brange(3:4)=[sridgef(round(0.05*(tn2-tn1+1))),sridgef(round(0.95*(tn2-tn1+1)))];
sfreq=sort(hfreq(1,~isnan(hfreq(1,:)))); mfreq=sfreq(round([0.05,0.95]*length(sfreq)));

%Check does the current harmonic really represent the first harmonic
if max(real(S0))>0 || max(imag(S0))>1
    if ~strcmpi(DispMode,'off'), fprintf('Checking does the extracted component represent first harmonic:\n'); end
    if iflag==1, hinfo.check_first=upd_info(hinfo.check_first,1,NaN,1,[NaN;NaN],[NaN;NaN],[NaN;NaN],[NaN;NaN],[NaN;NaN],wopt.f0,hamp(1,:),hphi(1,:),hfreq(1,:),tfsupp); end
    
    hid=1; itsig=sig-hamp(1,:).*cos(hphi(1,:)); fcounter=0; gwopt=wopt;
    hnmax=round(mean(hfreq(1,:))/(fs/L)); hnums=2:hnmax; if max(imag(S0))~=0, hnums=imag(S0); hnums=hnums(hnums>1); end
    for k=1:length(hnums)
        hn=hnums(k); if max(imag(S0))==0 && fcounter>=S0, break; end %stop if last [S0] harmonics are false (if this is the criterion)
        if ~strcmpi(DispMode,'off')
            [range1,range2]=harmfrange(1/hn,brange,wopt.f0,wopt.f0,fres,fs/L);
            fprintf('-- extracting harmonic 1/%d (basic frequency range [%0.3f,%0.3f] Hz)',hn,range2(1),range2(2));
            if ~strcmpi(DispMode,'on-'), fprintf(':\n'); end
        end
        
        [nrm,nmid,nf0,ncons,nac,npc,nfc,nslev,niamp,niphi,nifreq,ntfsupp,nTFR,nfreq,nwopt]=harmextract(1/hn,itsig,hamp(1,:),hphi(1,:),hfreq(1,:));
        if iflag==1, hinfo.check_first=upd_info(hinfo.check_first,hn,nrm,nmid,ncons,nac,npc,nfc,nslev,nf0,niamp,niphi,nifreq,ntfsupp); end
        
        %Check is harmonic true and if yes, subtract it and assign
        if nslev(nrm,nmid)>=slev && ncons(nrm,nmid)>=rhomin %assign the harmonic if true
            iamp=niamp; iphi=niphi; ifreq=nifreq;
            itsig=itsig(:)'-niamp(nrm,:).*cos(niphi(nrm,:)); %subtract the harmonic
            hid=[hid;hn]; hamp=vertcat(hamp,niamp(nrm,:)); hphi=vertcat(hphi,niphi(nrm,:)); hfreq=vertcat(hfreq,nifreq(nrm,:));
            sridgef=sort(ntfsupp(1,:)); brange(3:4)=[sridgef(round(0.05*(tn2-tn1+1))),sridgef(round(0.95*(tn2-tn1+1)))];
            fcounter=0; %reset the counter
            tfsupp=ntfsupp; TFR=nTFR; freq=nfreq; gwopt=nwopt; %update the first harmonic parameters
        else %update the counter of consecutive false harmonics
            fcounter=fcounter+1;
        end
        
        %Display if needed
        if ~strcmpi(DispMode,'off')
            if strcmpi(DispMode,'on-'), fprintf('\n'); end
            sstr1='<'; if ncons(nrm,nmid)>=rhomin, sstr1='>'; end
            sstr2='<'; if nslev(nrm,nmid)>=slev, sstr2='>'; end
            bstr='false'; if strcmpi(sstr1,'>') && strcmpi(sstr2,'>'), bstr='true'; end
            fprintf(['-- >> harmonic 1/',num2str(hn),' was identified as ',bstr,' (rho=',num2str(ncons(nrm,nmid)),...
                sstr1,num2str(rhomin),', signif=',num2str(nslev(nrm,nmid)),sstr2,num2str(slev),', f0=',num2str(nf0(nmid)),...
                ', amp. ratio = ',num2str(mean(niamp(nrm,:))/mean(hamp(1,:))),', phase shift = ',num2str(angle(mean(exp(1i*(hn*niphi(nrm,:)-hphi(1,:)))))/pi),'*Pi)\n']);
            if max(imag(S0))==0 && fcounter<S0 && hn==hnmax
                fprintf('-- next harmonic (1/%d) lies too high in frequency; aborting harmonic extraction.\n',hn+1);
            elseif max(imag(S0))~=0 && k<length(hnums) && hnums(k+1)>hnmax
                fprintf('-- next harmonic (1/%d) lies too high in frequency; aborting harmonic extraction.\n',hnums(k+1));
                break;
            end
        end
    end
    
    %Switch to the first harmonic if needed
    wopt=gwopt; clear ntfsupp nTFR nfreq nwopt;
    if ~strcmpi(DispMode,'off')
        if hid(end)==1, fprintf('-- it was confirmed that the specified component is the fundamental (first) harmonic.\n');
        else fprintf('-- the specified component is a %d harmonic; switching to the fundamental (first) harmonic.\n',hid(end)); end
    end
    hamp=hamp(end,:); hphi=hphi(end,:); hfreq=hfreq(end,:);
    
    %Find the minimal support around the first harmonic curve (to know frequency range for harmonics)
    [~,brange]=minsupp(tfsupp,TFR,freq,wopt,MinSupp); brange=[brange(:);brange(:)];
    sridgef=sort(tfsupp(1,:)); brange(3:4)=[sridgef(round(0.05*(tn2-tn1+1))),sridgef(round(0.95*(tn2-tn1+1)))];
end

%Extract all other harmonics
if ~strcmpi(DispMode,'off'), fprintf('Identifying and reconstructing higher harmonics:\n'); end
if iflag==1, hinfo=upd_info(hinfo,1,NaN,1,[NaN;NaN],[NaN;NaN],[NaN;NaN],[NaN;NaN],[NaN;NaN],wopt.f0,hamp(1,:),hphi(1,:),hfreq(1,:),tfsupp); end
hid=1; itsig=sig-hamp(1,:).*cos(hphi(1,:)); fcounter=0;
hnmax=round((fs/2)/mean(hfreq(1,:))); hnums=2:hnmax; if max(imag(S))~=0, hnums=imag(S); hnums=hnums(hnums>1); end
for k=1:length(hnums)
    hn=hnums(k); if max(imag(S))==0 && fcounter>=S, break; end %stop if last [S] harmonics are false (if this is the criterion)
    if ~strcmpi(DispMode,'off')
        [range1,range2]=harmfrange(hn,brange,wopt.f0,wopt.f0,fres,fs/L);
        fprintf('-- extracting harmonic %d (basic frequency range [%0.3f,%0.3f] Hz)',hn,range2(1),range2(2));
        if ~strcmpi(DispMode,'on-'), fprintf(':\n'); end
    end
    
    [nrm,nmid,nf0,ncons,nac,npc,nfc,nslev,niamp,niphi,nifreq,ntfsupp]=harmextract(hn,itsig,hamp(1,:),hphi(1,:),hfreq(1,:));
    if iflag==1, hinfo=upd_info(hinfo,hn,nrm,nmid,ncons,nac,npc,nfc,nslev,nf0,niamp,niphi,nifreq,ntfsupp); end
    
    %Check is harmonic true and if yes, subtract it and assign
    if nslev(nrm,nmid)>=slev && ncons(nrm,nmid)>=rhomin %assign the harmonic if true
        itsig=itsig(:)'-niamp(nrm,:).*cos(niphi(nrm,:)); %subtract the harmonic
        hid=[hid;hn]; hamp=vertcat(hamp,niamp(nrm,:)); hphi=vertcat(hphi,niphi(nrm,:)); hfreq=vertcat(hfreq,nifreq(nrm,:));
        sridgef=sort(ntfsupp(1,:)); brange(3:4)=[sridgef(round(0.05*(tn2-tn1+1))),sridgef(round(0.95*(tn2-tn1+1)))];
        fcounter=0; %reset the counter
    else %update the counter of consecutive false harmonics
        fcounter=fcounter+1;
    end
    
    %Display if needed
    if ~strcmpi(DispMode,'off')
        if strcmpi(DispMode,'on-'), fprintf('\n'); end
        sstr1='<'; if ncons(nrm,nmid)>=rhomin, sstr1='>'; end
        sstr2='<'; if nslev(nrm,nmid)>=slev, sstr2='>'; end
        bstr='false'; if strcmpi(sstr1,'>') && strcmpi(sstr2,'>'), bstr='true'; end
        fprintf(['-- >> harmonic ',num2str(hn),' was identified as ',bstr,' (rho=',num2str(ncons(nrm,nmid)),...
            sstr1,num2str(rhomin),', signif=',num2str(nslev(nrm,nmid)),sstr2,num2str(slev),', f0=',num2str(nf0(nmid)),...
            ', amp. ratio = ',num2str(mean(niamp(nrm,:))/mean(hamp(1,:))),', phase shift = ',num2str(angle(mean(exp(1i*(niphi(nrm,:)-hn*hphi(1,:)))))/pi),'*Pi)\n']);
        if max(imag(S))==0 && fcounter<S && hn==hnmax
            fprintf('-- next harmonic (%d) lies too high in frequency; aborting harmonic extraction.\n',hn+1);
        elseif max(imag(S))~=0 && k<length(hnums) && hnums(k+1)>hnmax
            fprintf('-- next harmonic (%d) lies too high in frequency; aborting harmonic extraction.\n',hnums(k+1));
            break;
        end
    end
end

if iflag==1, varargout{1}=hinfo; end



%============================= (Nested) ==========================================================================================================================================================================
%============= Function for extracting the specified harmonic ====================================================================================================================================================
%=================================================================================================================================================================================================================
    function [crm,cmid,cf0,ccons,cac,cpc,cfc,cslev,ciamp,ciphi,cifreq,ctfsupp,varargout]=harmextract(ch,csig,fhamp,fhphi,fhfreq)
        
        if NR>1
            czz=zeros(2,2*(NR+min([round(ch/ResAcc),mit])))*NaN;
            ccons=czz; cac=czz; cpc=czz; cfc=czz; cslev=czz; cf0=czz(1,:);
            cf0lim=wopt.f0*sort([1/ch,1]); if fres==2, cf0lim=ch*cf0lim; end
            if ~ischar(AdaptRange)
                car=AdaptRange; if isa(car,'function_handle'), car=car(ch); end
                if imag(car(1))~=0, cf0lim(1)=cf0lim(1)*imag(car(1)); else cf0lim(1)=car(1); end
                if imag(car(2))~=0, cf0lim(2)=cf0lim(2)*imag(car(2)); else cf0lim(2)=car(2); end
            end
            cf0(1:NR)=sort(linspace(cf0lim(1),cf0lim(2),NR));
            
            %--------- Search for the best resolution parameter -----------
            CDisp='off'; if strcmpi(DispMode,'plot++'), CDisp=DispMode; end
            if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-')
                astr='-- >> adapting the resolution parameter f0: '; estr=['-- >>',repmat(' ',1,length(astr)-5)];
                fprintf([astr,'testing ',num2str(NR),' values (from ',num2str(cf0(1)),' to ',num2str(cf0(NR)),') - ']); lstr=0;
            end
            %Crude search
            for rn=1:NR
                if strcmpi(DispMode,'on-'), fprintf('.');
                elseif ~strcmpi(DispMode,'off'), cstr=num2str(rn); fprintf([repmat('\b',1,lstr),cstr]); lstr=length(cstr); end
                [ccons(:,rn),cac(:,rn),cpc(:,rn),cfc(:,rn),cslev(:,rn),~,~,~,~]=hcalcall(ch,csig,cf0(rn),CDisp);
            end
            
            [~,cmid]=max(ccons,[],2);
            if ccons(1,cmid(1))>=ccons(2,cmid(2)), crm=1; else crm=2; end
            if strcmpi(RecMethod,'direct'), crm=1; end
            if strcmpi(RecMethod,'ridge'), crm=2; end
            cmid=cmid(crm);
            
            %Investigate further if consistency grows near the ends
            if erflag==1
                if ccons(crm,1)>=ccons(crm,2) || (cslev(crm,1)>=slev && cslev(crm,2)<slev) || NR==2 || (NR>2 && ccons(crm,2)>=ccons(crm,3))
                    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-'), fprintf(['\n',estr,'exploring f0 below ',num2str(min(cf0)),', f0 = ']); end
                    idnn=find(~isnan(cf0)); nnid=idnn(end:-1:1); cf0(idnn)=cf0(nnid);
                    ccons(:,idnn)=ccons(:,nnid); cac(:,idnn)=cac(:,nnid); cpc(:,idnn)=cpc(:,nnid); cfc(:,idnn)=cfc(:,nnid); cslev(:,idnn)=cslev(:,nnid);
                    while 1==1
                        cn=length(cf0(~isnan(cf0)))+1; cf0(cn)=cf0(cn-1)/(1+1/NR);
                        if strcmpi(DispMode,'on-'), fprintf('.');
                        elseif ~strcmpi(DispMode,'off'), fprintf([num2str(cf0(cn)),'; ']); end
                        [ccons(:,cn),cac(:,cn),cpc(:,cn),cfc(:,cn),cslev(:,cn),~,~,~]=hcalcall(ch,csig,cf0(cn),CDisp);
                        if ccons(crm,cn-2)>ccons(crm,cn-1) && ccons(crm,cn-2)>ccons(crm,cn), break; end
                    end
                    idnn=find(~isnan(cf0)); nnid=idnn(end:-1:1); cf0(idnn)=cf0(nnid);
                    ccons(:,idnn)=ccons(:,nnid); cac(:,idnn)=cac(:,nnid); cpc(:,idnn)=cpc(:,nnid); cfc(:,idnn)=cfc(:,nnid); cslev(:,idnn)=cslev(:,nnid);
                end
                idnn=find(~isnan(cf0));
                if ccons(crm,idnn(end))>=ccons(crm,idnn(end-1)) || (cslev(crm,idnn(end))>=slev && cslev(crm,idnn(end-1))<slev) || length(idnn)==2 || (length(idnn)>2 && ccons(crm,idnn(end-1))>=ccons(crm,idnn(end-2)))
                    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-'), fprintf(['\n',estr,'exploring f0 above ',num2str(max(cf0)),', f0 = ']); end
                    while 1==1
                        cn=length(cf0(~isnan(cf0)))+1; cf0(cn)=cf0(cn-1)*(1+1/NR);
                        if strcmpi(DispMode,'on-'), fprintf('.');
                        elseif ~strcmpi(DispMode,'off'), fprintf([num2str(cf0(cn)),'; ']); end
                        [ccons(:,cn),cac(:,cn),cpc(:,cn),cfc(:,cn),cslev(:,cn),~,~,~]=hcalcall(ch,csig,cf0(cn),CDisp);
                        if ccons(crm,cn-2)>ccons(crm,cn-1) && ccons(crm,cn-2)>ccons(crm,cn), break; end
                    end
                end
            end
            
            %Find the highest significant peak
            idnn=find(~isnan(cf0)); [qcons,qid]=sort(ccons(:,idnn),2,'descend');
            cmid1=find(cslev(1,qid(1,:))>=slev & ccons(1,qid(1,:))>=rhomin,1,'first'); cmv1=qcons(1,cmid1);
            cmid2=find(cslev(2,qid(2,:))>=slev & ccons(2,qid(2,:))>=rhomin,1,'first'); cmv2=qcons(2,cmid2);
            if strcmpi(RecMethod,'direct'), crm=1; cmid=qid(1,cmid1);
            elseif strcmpi(RecMethod,'ridge'), crm=2; cmid=qid(2,cmid2);
            else %automatically determine the method
                if isempty(cmv1) || isempty(cmv2)
                    if isempty(cmv1) && isempty(cmv2)
                        cmid=[]; if qcons(1,1)>=qcons(2,1), crm=1; else crm=2; end
                    elseif isempty(cmv1), crm=2; cmid=qid(2,cmid2);
                    elseif isempty(cmv2), crm=1; cmid=qid(1,cmid1);
                    end
                elseif cmv1>=cmv2, crm=1; cmid=qid(1,cmid1);
                else crm=2; cmid=qid(2,cmid2);
                end
            end
            if isempty(cmid), cmid=qid(crm,1); end
            if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-')
                mstr='direct'; if crm==2, mstr='ridge'; end
                fprintf(['\n',estr,'maximum found at f0=',num2str(cf0(cmid)),' (rho=',num2str(ccons(crm,cmid)),...
                    ', signif=',num2str(cslev(crm,cmid)),', method=',mstr,')\n']);
            end
            
            %Refine search
            if cslev(crm,cmid)>=slev
                if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-')
                    fprintf([estr,'approaching precise maximum (accuracy=',num2str(ResAcc),', maxiter=',num2str(mit),'), iterations - ']); lstr=0;
                end
                
                %Golden section search
                r=(1+sqrt(5))/2; cn0=length(cf0(~isnan(cf0))); cn=cn0+1;
                x=zeros(1,3); y=zeros(1,3); %initialization
                x(2)=cf0(cmid); y(2)=ccons(crm,cmid);
                if cmid>1, x(1)=cf0(cmid-1); y(1)=ccons(crm,cmid-1); end
                if cmid<idnn(end), x(3)=cf0(cmid+1); y(3)=ccons(crm,cmid+1); end
                if cmid==1, x(1)=x(2)-(x(3)-x(2)); y(1)=y(3); end
                if cmid==idnn(end), x(3)=x(2)+(x(2)-x(1)); y(3)=y(1); end
                itn=0;
                while (x(3)-x(1)>ResAcc && itn<mit)
                    if strcmpi(DispMode,'on-'), fprintf('.');
                    elseif ~strcmpi(DispMode,'off'), cstr=num2str(cn-cn0); fprintf([repmat('\b',1,lstr),cstr]); lstr=length(cstr); end
                    if x(3)-x(2)>x(2)-x(1)
                        if r*(x(2)-x(1))<0.9*(x(3)-x(2)), xx=x(2)+r*(x(2)-x(1));
                        else xx=x(2)+(x(2)-x(1))/r; end
                    else
                        if r*(x(3)-x(2))<0.9*(x(2)-x(1)), xx=x(2)-r*(x(3)-x(2));
                        else xx=x(2)-(x(3)-x(2))/r; end
                    end
                    [ccons(:,cn),cac(:,cn),cpc(:,cn),cfc(:,cn),cslev(:,cn),~,~,~]=hcalcall(ch,csig,xx,CDisp);
                    yy=ccons(crm,cn);
                    if yy>y(2)
                        if xx<x(2), x=[x(1),xx,x(2)]; y=[y(1),yy,y(2)];
                        else x=[x(2),xx,x(3)]; y=[y(2),yy,y(3)]; end
                    else
                        if xx<x(2), x=[xx,x(2),x(3)]; y=[yy,y(2),y(3)];
                        else x=[x(1),x(2),xx]; y=[y(1),y(2),yy]; end
                    end
                    cf0(cn)=xx; cn=cn+1; itn=itn+1;
                end
                
                %Final display
                if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'on-')
                    mstr='direct'; if crm==2, mstr='ridge'; end
                    idnn=find(~isnan(cf0)); idsf=idnn(cslev(crm,idnn)>=slev); [~,cmid]=max(ccons(crm,idsf)); cmid=idsf(cmid);
                    fprintf(['\n',estr,'precise maximum found at f0=',num2str(cf0(cmid)),' (rho=',num2str(ccons(crm,cmid)),...
                        ', signif=',num2str(cslev(crm,cmid)),', method=',mstr,')\n']);
                end
            end
            
            %Final estimates
            idnn=find(~isnan(cf0)); cf0=cf0(idnn); ccons=ccons(:,idnn); cac=cac(:,idnn); cpc=cpc(:,idnn); cfc=cfc(:,idnn); cslev=cslev(:,idnn);
            [cf0,sid]=sort(cf0); ccons=ccons(:,sid); cac=cac(:,sid); cpc=cpc(:,sid); cfc=cfc(:,sid); cslev=cslev(:,sid);
            idsf=find(cslev(crm,:)>=slev); if isempty(idsf), idsf=1:length(cf0); end
            [~,cmid]=max(ccons(crm,idsf)); cmid=idsf(cmid);
            [~,~,~,~,~,ciamp,ciphi,cifreq,ctfsupp,cTFR,cfreq,cwopt]=hcalcall(ch,csig,cf0(cmid),DispMode);
            
            %Plot if needed
            if ~isempty(strfind(lower(DispMode),'plot'))
                figure; axes('YLim',[0,1],'NextPlot','Add');
                if ch>1, title(['Harmonic ',num2str(ch)]); else title(['Harmonic 1/',num2str(round(1/ch))]); end
                lines=zeros(20,1); tcons=ccons;
                lines(1)=plot(cf0,tcons(1,:),'--b','LineWidth',2,'DisplayName','direct (nonsignificant)');
                lines(2)=plot(cf0,tcons(2,:),'--r','LineWidth',2,'DisplayName','ridge (nonsignificant)');
                tcons(cslev<slev)=NaN;
                lines(3)=plot(cf0,tcons(1,:),'-b','Marker','.','MarkerSize',15,'LineWidth',2,'DisplayName','direct (significant)');
                lines(4)=plot(cf0,tcons(2,:),'-r','Marker','.','MarkerSize',15,'LineWidth',2,'DisplayName','ridge (significant)');
                lines(5)=plot(cf0,cslev(1,:),':b','LineWidth',1,'Visible','off','DisplayName','significance (direct reconstruction)');
                lines(6)=plot(cf0,cslev(2,:),':r','LineWidth',1,'Visible','off','DisplayName','significance (ridge reconstruction)');
                lines(7)=plot(cf0,cac(1,:),'-b','LineWidth',1,'Visible','off','DisplayName','amplitude consistency (direct)');
                lines(8)=plot(cf0,cac(2,:),'-r','LineWidth',1,'Visible','off','DisplayName','amplitude consistency (ridge)');
                lines(9)=plot(cf0,cpc(1,:),'-b','LineWidth',1,'Visible','off','DisplayName','phase consistency (direct)');
                lines(10)=plot(cf0,cpc(2,:),'-r','LineWidth',1,'Visible','off','DisplayName','phase consistency (ridge)');
                lines(11)=plot(cf0,cfc(1,:),'-b','LineWidth',1,'Visible','off','DisplayName','frequency consistency (direct)');
                lines(12)=plot(cf0,cfc(2,:),'-r','LineWidth',1,'Visible','off','DisplayName','frequency consistency (ridge)');
                lines(13)=plot([min(cf0),max(cf0)],[rhomin,rhomin],':k','LineWidth',1,'DisplayName','threshold');
                legend(lines([3,4,1,2,13]));
            end
            
        else
            cmid=1; cf0=wopt.f0;
            [ccons,cac,cpc,cfc,cslev,ciamp,ciphi,cifreq,ctfsupp,cTFR,cfreq,cwopt]=hcalcall(ch,csig,cf0,DispMode);
            if strcmpi(RecMethod,'direct'), crm=1;
            elseif strcmpi(RecMethod,'ridge'), crm=2;
            else %automatically determine the method
                if cslev(1)>=slev && cslev(2)<slev, crm=1;
                elseif cslev(1)<slev && cslev(2)>=slev, crm=2;
                else
                    if ccons(1)>=ccons(2), crm=1; else crm=2; end
                end
            end
        end
        
        if nargout>12, varargout{1}=cTFR; end
        if nargout>13, varargout{2}=cfreq; end
        if nargout>14, varargout{3}=cwopt; end
        
        
%========================= (Nested in Nested) ====================================================================================================================================================================
%============= Function for calculating all for a given f0 =======================================================================================================================================================
%=================================================================================================================================================================================================================
        function [xcons,xac,xpc,xfc,xslev,xiamp,xiphi,xifreq,varargout]=hcalcall(xh,xsig,xf0,varargin)
            
            %Calculate TFR
            [xrange,yrange]=harmfrange(xh,brange,wopt.f0,xf0,fres,fs/L);
            xsig=bandpass(xsig,fs,xrange);
            if fres==1
                xfstep=wopt.fstepsim; if isempty(strfind(xfstep,'auto')), xfstep=xfstep*wopt.f0/xf0; end
                [xTFR,xfreq,xwopt]=wft(xsig,fs,wopt,'f0',xf0,'fmin',yrange(1),'fmax',yrange(2),'fstep',xfstep,'Display','off','Plot','off');
                xNF=length(xfreq); xhidf=1+floor(0.5+(xh*fhfreq-xfreq(1))/xwopt.fstep); %frequency indices for harmonic
            else
                xnv=wopt.nvsim; if isempty(strfind(xnv,'auto')), xnv=xnv*wopt.f0/xf0; end
                [xTFR,xfreq,xwopt]=wt(xsig,fs,wopt,'f0',xf0,'fmin',yrange(1),'fmax',yrange(2),'nv',xnv,'Display','off','Plot','off');
                xNF=length(xfreq); xhidf=1+floor(0.5+xwopt.nv*log(xh*fhfreq/xfreq(1))/log(2)); %frequency indices for harmonic
            end
            xhidf=max(min(xhidf,xNF),1);
            
            %Perform an initial preprocessing (to use fast extraction and reconstruction in the following)
            [Zsupp,Ztn,Zind,Zall]=partsupp(xTFR); Zsupp=xfreq(Zsupp); %partition the current TFR onto supports
            [Zamp,Zphi,Zfreq,Ztfsupp]=rectfr({Zsupp,Ztn},xTFR,xfreq,xwopt,'both'); %find the corresponding parameters
            
            %Original harmonic candidate
            xidz=Zall(xNF*(0:L-1)+xhidf); %find the z-indices for possible harmonic candidate
            xiamp=Zamp(:,xidz); xiphi=unwrap(Zphi(:,xidz),[],2); xifreq=Zfreq(:,xidz); xtfsupp=Ztfsupp(:,xidz); %reconstruct the parameters
            [xac,xpc,xfc]=acpcfc(xh,xiamp,xiphi,xifreq,fhamp,fhphi,fhfreq); xcons=(xac.^ww(1)).*(xpc.^ww(2)).*(xfc.^ww(3)); %calculate original consistencies
            
            %Surrogates and independence test
            xsac=zeros(2,NS+1); xspc=zeros(2,NS+1); xsfc=zeros(2,NS+1);
            MS=ceil(L/4); MS=MS+mod(MS,2); xid0=1+MS/2:L-MS/2;
            xtsh=[0,round([linspace(-MS,-max([MS/NS/2,1]),floor(NS/2)),linspace(max([MS/NS/2,1]),MS,ceil(NS/2))])]; %time shifts
            [xsac(:,1),xspc(:,1),xsfc(:,1)]=acpcfc(xh,xiamp(:,xid0),xiphi(:,xid0),xifreq(:,xid0),fhamp(xid0),fhphi(xid0),fhfreq(xid0));
            for sn=1:NS
                xid1=xid0-floor(xtsh(sn+1)/2); xid2=xid0+ceil(xtsh(sn+1)/2);
                xsidz=Zall(xNF*(xid2-1)+xhidf(xid1));
                [xsac(:,sn+1),xspc(:,sn+1),xsfc(:,sn+1)]=acpcfc(xh,Zamp(:,xsidz),Zphi(:,xsidz),Zfreq(:,xsidz),fhamp(xid1),fhphi(xid1),fhfreq(xid1));
            end
            xscons=(xsac.^ww(1)).*(xspc.^ww(2)).*(xsfc.^ww(3)); %overall consistencies
            xslev=zeros(2,1); for xrn=1:2, xslev(xrn)=length(find(xscons(xrn,2:end)<=xscons(xrn,1)))/NS; end
            
            %Optional outputs
            if nargout>8, varargout{1}=xtfsupp; end
            if nargout>9, varargout{2}=xTFR; end
            if nargout>10, varargout{3}=xfreq; end
            if nargout>11, varargout{4}=xwopt; end
            
            %Plot if needed
            if nargin>3 && ~isempty(strfind(lower(varargin{1}),'plot+'))
                scrsz=get(0,'ScreenSize'); sz1=0.014*scrsz(4); sz2=0.013*scrsz(4); 
                %Plot TFR and the parameter estimates
                figure('Position',[scrsz(3)/16,scrsz(4)/8,7*scrsz(3)/8,3*scrsz(4)/4]);
                axes('Position',[0.06,0.075,0.38,0.275],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add');
                xtt=(0:(L-1))/fs; if fres==2, set(gca,'YScale','log'); end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                XX=xtt; YY=xfreq; ZZ=abs(xTFR);
                if size(ZZ,2)>round(scrsz(3)), XX=linspace(XX(1),XX(end),round(scrsz(3))); end
                if size(ZZ,1)>round(scrsz(4)), YY=linspace(YY(1),YY(end),round(scrsz(4))); end
                ZZ=aminterp(xtt,xfreq,ZZ,XX,YY,'max'); XX=XX(:); YY=YY(:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                xhh=pcolor(XX,YY,ZZ); set(xhh,'LineStyle','none');
                xlabel('Time (s)','FontSize',sz1); ylabel('Frequency (Hz)','FontSize',sz1);
                xtstr='FALSE'; if ~isempty(find(xslev>=slev & xcons>=rhomin,1)), xtstr='TRUE'; end
                title({['TFR (f0 = ',num2str(xf0),') for harmonic ',num2str(xh),' (',xtstr,')'],...
                    ['(direct method: rho=',num2str(xcons(1),'%0.3f'),', signif=',num2str(xslev(1),'%0.3f'),', a_h=',num2str(mean(xiamp(1,:)./fhamp),'%0.3f'),', \phi_h/\pi=',num2str(angle(mean(exp(1i*(xiphi(1,:)-xh*fhphi))))/pi,'%0.3f'),')'],...
                    ['( ridge method: rho=',num2str(xcons(2),'%0.3f'),', signif=',num2str(xslev(2),'%0.3f'),', a_h=',num2str(mean(xiamp(2,:)./fhamp),'%0.3f'),', \phi_h/\pi=',num2str(angle(mean(exp(1i*(xiphi(2,:)-xh*fhphi))))/pi,'%0.3f'),')']},...
                    'FontSize',sz1);
                xlim([xtt(1),xtt(end)]); ylim([xfreq(1),xfreq(end)]);
                plot(xtt,xtfsupp(1,:),'-k','LineWidth',2,'DisplayName','ridge curve');
                plot(xtt,xtfsupp(2,:),'Color',[0.6,0.6,0.6],'LineWidth',1,'DisplayName','lower bound of time-frequency support');
                plot(xtt,xtfsupp(3,:),'Color',[0.6,0.6,0.6],'LineWidth',1,'DisplayName','upper bound of time-frequency support');
                plot(xtt,xh*fhfreq,'--k','LineWidth',2,'Color',[0.4,0.4,0.4],'DisplayName','expected frequency profile');
                
                axes('Position',[0.06,0.5,0.38,0.125],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add','XTickLabel',{});
                ylabel('Amplitude','FontSize',sz1); xlines=zeros(10,1);
                xlines(1)=plot(xtt,fhamp*mean(xiamp(1,:)./fhamp),'-','Color',[0.6,0.6,0.6],'LineWidth',3,'DisplayName','direct amplitude (expected)');
                xlines(2)=plot(xtt,fhamp*mean(xiamp(2,:)./fhamp),'-','Color',[0.6,0.6,0.6],'LineWidth',3,'DisplayName','ridge amplitude (expected)');
                xlines(3)=plot(xtt,xiamp(1,:),'-b','LineWidth',2,'DisplayName','direct amplitude (obtained)');
                xlines(4)=plot(xtt,xiamp(2,:),'-r','LineWidth',2,'DisplayName','ridge amplitude (obtained)');
                axis tight; if xcons(1)<xcons(2), set(xlines(1),'Visible','off'); else set(xlines(2),'Visible','off'); end
                
                axes('Position',[0.06,0.65,0.38,0.125],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add','XTickLabel',{});
                ylabel({'Detrended','phase'},'FontSize',sz1); xlines=zeros(10,1);
                xlines(1)=plot(xtt,detrend(xh*unwrap(fhphi)+angle(mean(exp(1i*(xiphi(1,:)-xh*fhphi))))),'-','Color',[0.6,0.6,0.6],'LineWidth',3,'DisplayName','direct phase (expected)');
                xlines(2)=plot(xtt,detrend(xh*unwrap(fhphi)+angle(mean(exp(1i*(xiphi(2,:)-xh*fhphi))))),'-','Color',[0.6,0.6,0.6],'LineWidth',3,'DisplayName','ridge phase (expected)');
                xlines(3)=plot(xtt,detrend(unwrap(xiphi(1,:))),'-b','LineWidth',2,'DisplayName','direct phase (obtained)');
                xlines(4)=plot(xtt,detrend(unwrap(xiphi(2,:))),'-r','LineWidth',2,'DisplayName','ridge phase (obtained)');
                xwphi=mod(xiphi-ones(2,1)*xh*fhphi,2*pi); if abs(angle(mean(exp(1i*(xiphi(2,:)-xh*fhphi)))))<pi/2, xwphi=-pi+mod(xiphi-ones(2,1)*xh*fhphi+pi,2*pi); end
                xlines(5)=plot(xtt,xwphi(1,:),'--b','LineWidth',1,'DisplayName','wrapped direct phase difference');
                xlines(6)=plot(xtt,xwphi(2,:),'--r','LineWidth',1,'DisplayName','wrapped ridge phase difference');
                axis tight; if xcons(1)<xcons(2), set(xlines(1),'Visible','off'); else set(xlines(2),'Visible','off'); end
                
                axes('Position',[0.06,0.8,0.38,0.125],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add','XTickLabel',{});
                ylabel({'Instantaneous','Frequency'},'FontSize',sz1); xlines=zeros(10,1);
                xlines(2)=plot(xtt,xh*fhfreq,'-','Color',[0.6,0.6,0.6],'LineWidth',3,'DisplayName','direct frequency (expected)');
                xlines(1)=plot(xtt,xifreq(1,:),'-b','LineWidth',2,'DisplayName','direct frequency (obtained)');
                xlines(3)=plot(xtt,xifreq(2,:),'-r','LineWidth',2,'DisplayName','ridge frequency (obtained)');
                xlines(4)=plot(xtt,ones(size(xtt))*NaN,'--k','LineWidth',1,'Visible','off');
                axis tight;
                xhl=legend(xlines([1,3,2,4]),'direct','ridge','expected','phase difference (wrapped)','Orientation','horizontal','Location','NorthOutside');
                xpos=get(xhl,'Position'); set(xhl,'Position',[xpos(1),0.95,xpos(3),xpos(4)]);
                
                %Plot surrogate values and histograms
                annotation(gcf,'line',[0.5,0.5],[0,1],'LineWidth',3);
                for xrn=1:2
                    for yrn=1:4
                        xax1=axes('Position',[0.57+0.235*(xrn-1),0.75-0.225*(yrn-1),0.13,0.175],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add');
                        if xrn==1 && yrn==1, title({'Direct consistencies',['(\alpha=',num2str(ww(1)),', \beta=',num2str(ww(2)),', \gamma=',num2str(ww(3)),')']},'FontSize',sz1); end
                        if xrn==2 && yrn==1, title({'Ridge consistencies',['(\alpha=',num2str(ww(1)),', \beta=',num2str(ww(2)),', \gamma=',num2str(ww(3)),')']},'FontSize',sz1); end
                        if xrn==1
                            if yrn==2, ylabel('\rho_A (amplitude)','FontSize',sz1);
                            elseif yrn==3, ylabel('\rho_\phi (phase)','FontSize',sz1);
                            elseif yrn==4, ylabel('\rho_\nu (frequency)','FontSize',sz1);
                            else ylabel({'\rho=\rho_A^\alpha\rho_\phi^\beta\rho_\nu^\gamma','(overall)'},'FontSize',sz1);
                            end
                        end
                        xax2=axes('Position',[0.71+0.235*(xrn-1),0.75-0.225*(yrn-1),0.04,0.175],'Layer','top','Box','on','FontSize',sz2,'NextPlot','Add');
                        
                        if yrn==2, xxsval=xsac(xrn,:); elseif yrn==3, xxsval=xspc(xrn,:);
                        elseif yrn==4, xxsval=xsfc(xrn,:); else xxsval=xscons(xrn,:); end
                        xxslev=length(find(xxsval(2:end)<=xxsval(1)))/NS;
                        xxval=xxsval(1); [xxtsh,xxid]=sort(xtsh); xxsval=xxsval(xxid);
                        
                        axes(xax1);
                        plot(xxtsh/fs,xxsval,'-k','Marker','.','MarkerSize',12);
                        plot(0,xxval,'-r','Marker','.','MarkerSize',30); plot([xxtsh(1)/fs,xxtsh(end)/fs],[xxval,xxval],'--r','LineWidth',2);
                        axis tight; ymym=get(gca,'YLim'); set(gca,'XTickLabel',{});
                        
                        axes(xax2);
                        title(['p=',num2str(xxslev),''],'FontSize',sz1); view(-90,90); set(gca,'YDir','reverse');
                        hist(xxsval(xxtsh~=0),ceil(sqrt(NS))); axis tight; set(gca,'XLim',ymym);
                        plot([xxval,xxval],get(gca,'YLim'),'--r','LineWidth',2,'DisplayName','original value');
                        set(gca,'XTickLabel',{},'YTickLabel',{});
                    end
                    axes(xax1); set(gca,'XTickLabelMode','auto'); xlabel('Time (s)','FontSize',sz1);
                end
            end
            
        end
        
    end


end



%=================================================================================================================================================================================================================
%========================== Other functions ======================================================================================================================================================================
%=================================================================================================================================================================================================================

%Function for partitioning the given TFR onto supports:
%Zsupp: 2xSL matrix with frequency indices of the supports
%Ztn: 1xSL vector with time-indices corresponding to [Zsupp]
%Zind: 2xL matrix containing range indices for [Zsupp] corresponding to each time
%(SL - total number of supports, L - length of a signal)
function [Zsupp,Ztn,Zind,Zall]=partsupp(TFR)

[NF,L]=size(TFR); TFR=vertcat(0*ones(1,L),-1*ones(1,L),abs(TFR),-1*ones(1,L),0*ones(1,L));
idft=1+find(TFR(2:end-1)<=TFR(1:end-2) & TFR(2:end-1)<TFR(3:end)); %find linear indices of the minimums
[idf,idt]=ind2sub(size(TFR),idft); idf=idf(:)-2; idt=idt(:); %find frequency and time indices of the peaks
idf(idf<1)=1; idf(idf>NF)=NF; %correct the "wrong" indices
SL=length(idf); dind=[0;find(diff(idt)>0);SL];

Ztn=idt(1:end-1)'; Ztn(dind(2:end-1))=[]; ZL=length(Ztn);
Zsupp=[idf(1:end-1)';idf(2:end)']; Zsupp(:,dind(2:end-1))=[];
Zind=[0,find(diff(Ztn)>0),ZL]; Zind=[Zind(1:end-1)+1;Zind(2:end)];

%Form the matrix of indices
Zall=ones(NF,L); Zall(1,:)=Zind(1,:); nadd=Zsupp(2,:)-Zsupp(1,:);
n2=Ztn+cumsum(nadd); n1=n2-nadd+1; for zn=1:ZL, Zall(n1(zn):n2(zn))=zn; end

end

%Function for calculating the frequency ranges for harmonics:
%xrange - basic range where the harmonic FT resides;
%yrange - range in which to calculate TFR (taking f0 into account);
function [xrange,yrange]=harmfrange(hn,brange,af0,bf0,fres,T)

if fres==1
    xrange=hn*mean(brange(1:2))+[-1/2,1/2]*diff(brange(1:2))*max([1,hn]);
    xmean=mean(xrange); xd=[xrange(1)-xmean,xrange(2)-xmean];
    if hn>1 && xrange(1)<brange(3), xd(1)=brange(3)-xmean; end
    if hn<1 && xrange(2)>brange(4), xd(2)=brange(4)-xmean; end
    yrange=xmean+xd*max(1,af0/bf0/max([1,hn]));
else
    xrange=hn*mean(brange(1:2))+[-1/2,1/2]*diff(brange(1:2))*max([1,hn]);
    xmean=mean(xrange); xd=[xrange(1)-xmean,xrange(2)-xmean];
    if hn>1 && xrange(1)<brange(3), xd(1)=brange(3)-xmean; end
    if hn<1 && xrange(2)>brange(4), xd(2)=brange(4)-xmean; end
    if xmean+xd(1)<1/T, xd(1)=min([1/T,xmean/2])-xmean; end
    yrange=(xmean+xd).*((xmean+xd)/xmean).^(max([1,af0/bf0])-1);
    if yrange(1)<1/T/max([1,1/hn]), yrange(1)=min([1/T/max([1,1/hn]),xmean/2]); end
end

end

%Function for calculating all consistencies
function [ac,pc,fc,varargout]=acpcfc(h,hamp,hphi,hfreq,iamp,iphi,ifreq,varargin)

[NR,L]=size(hamp);
if size(iamp,2)~=L, iamp=iamp'; end
if size(iphi,2)~=L, iphi=iphi'; end
if size(ifreq,2)~=L, ifreq=ifreq'; end

NS=0; if nargin>7, NS=varargin{1}; end
MS=ceil(L/4); MS=MS+mod(MS,2); id0=1+MS/2:L-MS/2;
tsh=[0,round([linspace(-MS,-max([MS/NS/2,1]),floor(NS/2)),linspace(max([MS/NS/2,1]),MS,ceil(NS/2))])]; %time shifts
ac=zeros(NR,NS+1); pc=zeros(NR,NS+1); fc=zeros(NR,NS+1); id1=1:L; id2=1:L;
for rn=1:NR
    for sn=1:NS+1
        if NS>0, id1=id0-floor(tsh(sn)/2); id2=id0+ceil(tsh(sn)/2); end
        %ac(rn,sn)=exp(-sqrt(mean((-1+(hamp(rn,id2)./iamp(1,id1))/mean(hamp(rn,id2)./iamp(1,id1))).^2)));
        %ac(rn,sn)=exp(-sqrt(mean((-1+(iamp(1,id1)./hamp(rn,id2))/mean(iamp(1,id1)./hamp(rn,id2))).^2)));
        ac(rn,sn)=exp(-sqrt(mean((hamp(rn,id2)*mean(iamp(1,id1))-iamp(1,id1)*mean(hamp(rn,id2))).^2))/mean(iamp(1,id1).*hamp(rn,id2)));
        pc(rn,sn)=abs(mean(exp(1i*(max([1,1/h])*hphi(rn,id2)-max([1,h])*iphi(1,id1)))));
        %if h>=1, fc(rn,sn)=exp(-sqrt(mean((-1+h*ifreq(1,id1)./hfreq(rn,id2)).^2)));
        %else fc(rn,sn)=exp(-sqrt(mean((-1+hfreq(rn,id2)./(h*ifreq(1,id1))).^2))); end
        fc(rn,sn)=exp(-sqrt(mean((hfreq(rn,id2)-h*ifreq(1,id1)).^2))/mean(hfreq(rn,id2)));
    end
end

if nargout>3, varargout{1}=tsh; end

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

%------ Function for updating the structure with full information ---------
function hinfo=upd_info(hinfo,hn,nrm,nmid,ncons,nac,npc,nfc,nslev,nf0,niamp,niphi,nifreq,ntfsupp)

CH=length(hinfo.rm);

hinfo.hnum=[hinfo.hnum,hn];
hinfo.cons=[hinfo.cons,ncons(:,nmid)];
hinfo.signif=[hinfo.signif,nslev(:,nmid)];
hinfo.f0=[hinfo.f0,nf0(nmid)];
hinfo.ac=[hinfo.ac,nac(:,nmid)];
hinfo.pc=[hinfo.pc,npc(:,nmid)];
hinfo.fc=[hinfo.fc,nfc(:,nmid)];
hinfo.rm=[hinfo.rm,nrm];

hinfo.iamp{CH+1}=niamp;
hinfo.iphi{CH+1}=niphi;
hinfo.ifreq{CH+1}=nifreq;
hinfo.tfsupp{CH+1}=ntfsupp;

hinfo.adapt_all.mid=[hinfo.adapt_all.mid,nmid];
hinfo.adapt_all.cons{CH+1}=ncons;
hinfo.adapt_all.signif{CH+1}=nslev;
hinfo.adapt_all.f0{CH+1}=nf0;
hinfo.adapt_all.ac{CH+1}=nac;
hinfo.adapt_all.pc{CH+1}=npc;
hinfo.adapt_all.fc{CH+1}=nfc;

end

%=============== Interpolation function for plotting ======================
% ZI=aminterp(X,Y,Z,XI,YI,method)
% - the same as interp2, but uses different interpolation types,
% maximum-based ([method]='max') or average-based ([method]='avg'), where
% [ZI] at each point [XI,YI] represents maximum or average among values of
% [Z] corresponding to the respective quadrant in [X,Y]-space which
% includes point [XI,YI];
% X,Y,XI,YI are all linearly spaced vectors, and the lengths of XI,YI
% should be the same or smaller than that of X,Y; Z should be real,
% ideally positive; X and Y correspond to the 2 and 1 dimension of Z, as
% always.
function ZI = aminterp(X,Y,Z,XI,YI,method)

%Interpolation over X
ZI=zeros(size(Z,1),length(XI))*NaN;
xstep=mean(diff(XI));
xind=1+floor((1/2)+(X-XI(1))/xstep); xind=xind(:);
xpnt=[0;find(xind(2:end)>xind(1:end-1));length(xind)];
if strcmpi(method,'max')
    for xn=1:length(xpnt)-1
        xid1=xpnt(xn)+1; xid2=xpnt(xn+1);
        ZI(:,xind(xid1))=max(Z(:,xid1:xid2),[],2);
    end
else
    for xn=1:length(xpnt)-1
        xid1=xpnt(xn)+1; xid2=xpnt(xn+1);
        ZI(:,xind(xid1))=mean(Z(:,xid1:xid2),2);
    end
end

Z=ZI;

%Interpolation over Y
ZI=zeros(length(YI),size(Z,2))*NaN;
ystep=mean(diff(YI));
yind=1+floor((1/2)+(Y-YI(1))/ystep); yind=yind(:);
ypnt=[0;find(yind(2:end)>yind(1:end-1));length(yind)];
if strcmpi(method,'max')
    for yn=1:length(ypnt)-1
        yid1=ypnt(yn)+1; yid2=ypnt(yn+1);
        ZI(yind(yid1),:)=max(Z(yid1:yid2,:),[],1);
    end
else
    for yn=1:length(ypnt)-1
        yid1=ypnt(yn)+1; yid2=ypnt(yn+1);
        ZI(yind(yid1),:)=mean(Z(yid1:yid2,:),1);
    end
end

end

%------------------------- Bandpass filtering -----------------------------
function bsig = bandpass(sig,fs,fband)

L=length(sig); Nq=ceil((L+1)/2); %signal's length and the Nyquist index
ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; %frequencies in Fourier transform

fx=fft(sig); %Fourier transform of the signal
fx(abs(ff)<fband(1) | abs(ff)>fband(2))=0; %filtering
bsig=ifft(fx); %recover the filteres signal

end