%====================== Nonlinear Mode Decomposition ======================
% Version 1.00 stable
%------------------------------Requirements--------------------------------
% wft.m, wt.m, ecurve.m, rectfr.m, bestest.m, eharm.m, recnm.m, noisetest.m
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
% {preprint - arXiv:1310.7215}
% [2] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part II: Resolution, reconstruction and concentration."
% {preprint - arXiv:1310.7274}
% [3] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "On the extraction of instantaneous frequencies from ridges in
%  time-frequency representations of signals."
% {preprint - arXiv:1310.7276}
% [4] D. Iatsenko, D. Iatsenko, P.V.E. McClintock and A. Stefanovska,
% "Nonlinear Mode Decomposition: a new noise-robust, adaptive
%  decomposition method."
% {preprint - arXiv:1207.5567}
%
%------------------------------Documentation-------------------------------
%
% [NM,Optional:info,hamp,hphi,hfreq,hid,ah,ph]=nmd(sig,fs,Optional:'PropertyName',PropertyValue);
% - applies Nonlinear Mode Decomposition procedure to signal [sig] which
%   is sampled at [fs] Hz. Returns all physically meaningful Nonlinear
%   Modes in the rows of [NM] (the residual is then [sig(:)'-sum(NM,1)]).
%
% INPUT:
% sig: 1xL or Lx1 vector
%        - the signal from which to extract the harmonics
% fs: value
%        - sampling frequency of the signal
%
% OUTPUT:
% NM: MxL matrix
%        - matrix containing in its rows all extracted Nonlinear Modes that
%          were identified as physically meaningful (i.e. not noise)
% info: structure
%        - structure with all parameters and outcome of the decomposition
% hamp,hphi,hfreq: Mx1 cell (with [m]th cell containing [H(m) x L] matrix)
%        - estimated instantaneous amplitudes [hamp], phases [hphi] and
%          frequencies [hfreq] of each of the [H(m)] harmonics
%          (corresponding to rows of [hamp{m}], [hphi{m}] and [hfreq{m}])
%          associated with each [m=1,...,M] Nonlinear Mode in [NM] (so
%          that [NM(m,:)=sum(hamp{m}.*cos(hphi{m}),1)]).
% hid,ah,ph: Mx1 cell (with [m]th cell containing [H(m) x 1] vector)
%        - harmonic numbers [hid] (e.g. [1;3;4;7]), amplitude ratios [ah]
%          and phase shifts [ph] for each Nonlinear Mode in [NM].
%
% PROPERTIES: ({{...}} denotes default)
% ################################ BASIC ##################################
% 'Display':{{'on'}}|'on+'|'off'
%        - defines the level of displayed information.
% 'Template':{{'Accurate'}}|'Fast'|'VeryAccurate'|'VeryFast'
%        - the template determining set of parameters which control the
%          speed/accuracy tradeoff (the parameters and their values
%          corresponding to each template can be seen below).
% 'ModeNum':{{'auto'}}|number
%        - the number of modes to extract; the default auto determines it
%          automatically, stopping decomposition when the residual is
%          regarded as noise (i.e. does not pass the associate surrogate
%          test implemented in <noisetest.m>).
% 'TFRtype':{{'auto-WT'}}|'auto-WFT'|'WT'|'WFT'
%        - the type of TFR to be used in all calculations; for 'auto-WT'
%          or 'auto-WFT' determines it automatically, but at the first
%          step uses the wavelet transform (WT) or windowed Fourier
%          transform (WFT), respectively.
% ############################## ADVANCED #################################
% All properties of <wft.m>, <wt.m>, <ecurve.m>, <bestest.m>, <eharm.m> and
% <noisetest.m>, with the same defaults, except:
% - 'Display' (all functions) is always determined by the level of display
%   set by the corresponding property of the NMD function (<nmd.m>);
% - 'Plot' (<wft.m>, <wt.m>, <ecurve.m>) is always set to 'off';
% - 'CutEdges' (<wft.m>, <wt.m>) is always set to 'off';
% - 'Skel' (<ecurve.m>) is always set to [] (empty matrix);
% - 'Normalize' (<ecurve.m>) by default is set to 'on'.
%
%
% NOTE: One can alternatively pass the structure with the properties as the
% 3rd argument, e.g. /opt.ModeNum='auto'; nmd(sig,fs,opt);/. If the other
% properties are specified next, they override those in the structure,
% e.g. /nmd(sig,fs,opt,'ModeNum',4);/ will always extract M=4 modes,
% irrespectively to what is specified in [opt].
%
%-------------------------------Examples-----------------------------------
%
% NM=nmd(sig,fs);
% - extracts all meaningful modes using default parameters.
% NM=nmd(sig,fs,'Template','Fast');
% - extracts all meaningful modes using 'Fast' parameter template.
% NM=nmd(sig,fs,'Template','Fast','ModeNum',1);
% - extracts one mode (not necessarily meaningful) using 'Fast' template.
% NM=nmd(sig,fs,'S',4);
% - extracts all meaningful modes using default parameters, except setting
%   the 'S' property of <eharm.m> to 4, so that for each NM the harmonics
%   are extracted until 4 of them are identified as false in a row.
%
%==========================================================================

function [NM,varargout] = nmd(sig,fs,varargin)

L=length(sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default parameters
opts=struct('Display','on','ModeNum','auto','Template','Accurate','TFRtype','auto-WT');
opts.MinSupp='on'; opts.Normalize='on';
%Update if user-defined
vst=1;
while (nargin>=vst+2 && isstruct(varargin{vst}))
    copt=varargin{vst}; fnames=filednames(copt); vst=vst+1;
    for fn=1:length(fnames), opts.(fnames{fn})=copt.(fnames{fn}); end
end
for vn=vst:2:nargin-2
    opts.(varargin{vn})=varargin{vn+1};
end
%Apply the template
fnames={'S','S0','AdaptRes','AdaptRange','SurrNum'};
if strcmpi(opts.Template,'VeryFast'), fvals={2,2,[1,0,0,0],'auto',20}; end
if strcmpi(opts.Template,'Fast'), fvals={2,2,[3,1,0,3],'auto',40}; end
if strcmpi(opts.Template,'Accurate'), fvals={3,3,[10,1,1i/100,Inf],'auto',40}; end
if strcmpi(opts.Template,'VeryAccurate'), fvals={3,3,[20,1,1i/100,Inf],[1i*1,1i*2],80}; end
for fn=1:length(fnames), opts.(fnames{fn})=fvals{fn}; end
%Make some adjustments of properties
opts.CutEdges='off'; opts.Plot='off'; opts.Skel=[];
CDisp='off'; if strcmpi(opts.Display,'on'), CDisp='on-'; elseif ~strcmpi(opts.Display,'off'), CDisp='on'; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove the trend from the signal
X=(1:length(sig))'/fs; XM=ones(length(X),4); for pn=1:3, CX=X.^pn; XM(:,pn+1)=(CX-mean(CX))/std(CX); end
w=warning('off','all'); trend=XM*(pinv(XM)*sig(:)); sig=sig(:)'-trend(:)'; warning(w);

%Extract the Nonlinear Modes
miter=opts.ModeNum; if strcmpi(miter,'auto'), miter=Inf; end
M=ceil(log(L)); cc=cell(M,1); hamp=cc; hphi=cc; hfreq=cc; hid=cc; ah=cc; ph=cc;
NM=zeros(M,L); info=struct; info.opts=opts; info.harm=cc; info.signif=[];
itn=1;
while itn<=miter
    if ~strcmpi(opts.Display,'off')
        fprintf('=========================== Extracting NM %d ===========================\n',itn);
    end
    
    %Calculate the TFR
    if ~isempty(strfind(opts.TFRtype,'WT')), [TFR,freq,wopt]=wt(sig,fs,opts,'Display',CDisp); end
    if ~isempty(strfind(opts.TFRtype,'WFT')), [TFR,freq,wopt]=wft(sig,fs,opts,'Display',CDisp); end
    
    %Extract the component
    if ~strcmpi(opts.Display,'off'), fprintf('Extracting the component...'); end
    [tfsupp,ecinfo]=ecurve(TFR,freq,wopt,opts,'Display','off'); [ramp,rphi,rfreq]=rectfr(tfsupp,TFR,freq,wopt,'ridge');
    if ~strcmpi(opts.Display,'off'), fprintf(' done (ridge frequencies = %0.3f+-%0.3f, ridge amplitudes = %0.3f+-%0.3f)\n',mean(rfreq),std(rfreq),mean(ramp),std(ramp)); end
    
    %Test the extracted component against noise
    if ~isfinite(miter)
        [p,signif]=noisetest(sig,TFR,freq,wopt,tfsupp,ecinfo,opts,'Display',CDisp);
        info.signif=[info.signif;signif(:)'];
        if p==0
            if ~strcmpi(opts.Display,'off'), fprintf('Stopping the decomposition (number of extracted modes: %d).\n',itn-1); end
            break;
        end
    end
    
    %Determine the optimal TFR type for the extracted component
    tfrtype=opts.TFRtype; if ~isempty(strfind(tfrtype,'auto')), tfrtype=checktype(fs,ramp,rfreq,CDisp); end
    
    %Recalculate all (if optimal TFR type differs from the calculated one)
    if isempty(strfind(opts.TFRtype,tfrtype))
        [~,mrange]=minsupp(tfsupp,TFR,freq,wopt,opts.MinSupp); fmean=mean(rfreq);
        if strcmpi(tfrtype,'WT')
            frange=fmean*exp((mrange-fmean)/fmean); %frequency range for the WT
            [TFR,freq,wopt]=wt(bandpass(sig,fs,mrange),fs,opts,'f0',wopt.f0*fmean,'fmin',frange(1),'fmax',frange(2),'Display',CDisp);
        else
            frange=fmean*(1+log(mrange/fmean)); %frequency range for the WFT
            [TFR,freq,wopt]=wft(bandpass(sig,fs,mrange),fs,opts,'f0',wopt.f0/fmean,'fmin',frange(1),'fmax',frange(2),'Display',CDisp);
        end
        if ~strcmpi(opts.Display,'off'), fprintf('Extracting the component...'); end
        [tfsupp,ecinfo]=ecurve(TFR,freq,wopt,opts,'Display','off'); [ramp,rphi,rfreq]=rectfr(tfsupp,TFR,freq,wopt,'ridge');
        if ~strcmpi(opts.Display,'off'), fprintf(' done (ridge frequencies = %0.3f+-%0.3f, ridge amplitudes = %0.3f+-%0.3f)\n',mean(rfreq),std(rfreq),mean(ramp),std(ramp)); end
    end
    
    %Extract the harmonics
    [hid{itn},hamp{itn},hphi{itn},hfreq{itn},hinfo]=eharm(sig,tfsupp,TFR,freq,wopt,opts,'Display',CDisp);
    
    %Reconstruct the Nonlinear Mode
    [NM(itn,:),ah{itn},ph{itn}]=recnm(hid{itn},hamp{itn},hphi{itn},hfreq{itn});
    
    %Subtract NM from the signal and update the counter
    sig=sig-NM(itn,:); info.harm{itn}=hinfo; itn=itn+1;
end
ii=1:itn-1; NM=NM(ii,:); info.harm=info.harm(ii); hamp=hamp(ii); hphi=hphi(ii); hfreq=hfreq(ii); hid=hid(ii); ah=ah(ii); ph=ph(ii);

if nargout>1, varargout{1}=info; end
if nargout>2, varargout{2}=hamp; end
if nargout>3, varargout{3}=hphi; end
if nargout>4, varargout{4}=hfreq; end
if nargout>5, varargout{5}=hid; end
if nargout>6, varargout{6}=ah; end
if nargout>7, varargout{7}=ph; end

end


%==========================================================================
%======================== Supporting functions ============================
%==========================================================================

%------------ Function for determining the optimal TFR type ---------------
function tfrtype = checktype(fs,iamp,ifreq,DispMode)

%Different parameters
Perc=0.75; DLev=1.1; %percentile range and the discrimination threshold
L=length(iamp); i1=round((0.5-Perc/2)*L); i2=round((0.5+Perc/2)*L);

%Estimate localized frequencies/amplitudes and their time-derivatives
tamp=(iamp(1:end-2)+iamp(2:end-1)+iamp(3:end))/3; dtamp=fs*(iamp(3:end)-iamp(1:end-2))/2;
tfreq=(ifreq(1:end-2)+ifreq(2:end-1)+ifreq(3:end))/3; dtfreq=fs*(ifreq(3:end)-ifreq(1:end-2))/2;

%Calculate the values of [V]'s and the overall functional [U]
gamp1=sort(abs(hilbert(dtamp./tfreq-mean(dtamp./tfreq)))); gamp2=sort(abs(hilbert((dtamp-mean(dtamp))./mean(tfreq))));
gfreq1=sort(abs(hilbert(dtfreq./tfreq-mean(dtfreq./tfreq)))); gfreq2=sort(abs(hilbert((dtfreq-mean(dtfreq))./mean(tfreq))));
Vamp=(gamp1(i2)-gamp1(i1))/(gamp2(i2)-gamp2(i1)); Vfreq=(gfreq1(i2)-gfreq1(i1))/(gfreq2(i2)-gfreq2(i1));
if isnan(Vamp), Vamp=1; end, if isnan(Vfreq), Vfreq=1; end
U=1/(1+Vamp)+1/(1+Vfreq);

%Determine the optimal TFR type
if U<DLev, tfrtype='WFT'; cstr='<'; else tfrtype='WT'; cstr='>'; end
if ~strcmpi(DispMode,'off')
    fprintf(['Optimal TFR type was determined to be ',tfrtype,' (']);
    fprintf('Va=%0.3f, Vf=%0.3f, 1/(1+Va)+1/(1+Vf)=%0.3f',Vamp,Vfreq,U)
    fprintf([cstr,num2str(DLev),')\n']);
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

%------------------------- Bandpass filtering -----------------------------
function bsig = bandpass(sig,fs,fband)

L=length(sig); Nq=ceil((L+1)/2); %signal's length and the Nyquist index
ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L; %frequencies in Fourier transform

fx=fft(sig); %Fourier transform of the signal
fx(abs(ff)<fband(1) | abs(ff)>fband(2))=0; %filtering
bsig=ifft(fx); %recover the filteres signal

end
