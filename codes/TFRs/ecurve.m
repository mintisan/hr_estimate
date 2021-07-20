%================== Extract ridge curve from WFT or WT ====================
% Version 2.20 stable
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
% {preprint - arXiv:1310.7215}
% [2] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "Linear and synchrosqueezed time-frequency representations revisited.
%  Part II: Resolution, reconstruction and concentration."
% {preprint - arXiv:1310.7274}
% [3] D. Iatsenko, A. Stefanovska and P.V.E. McClintock,
% "On the extraction of instantaneous frequencies from ridges in
%  time-frequency representations of signals."
% {preprint - arXiv:1310.7276}
%
%------------------------------Documentation-------------------------------
%
% [tfsupp,Optional:ecinfo,Skel] = ecurve(TFR,freq,wopt,'PropertyName',PropertyValue)
% - extracts the curve (i.e. the sequence of the amplitude ridge points)
%   and its full support (the widest region of unimodal TFR amplitude
%   around them) from the given time-frequency representation [TFR] of the
%   signal. TFR can be either WFT or WT.
% 
% OUTPUT:
% tfsupp: 3xL matrix
%    - extracted time-frequency support of the component, containing
%      frequencies of the TFR amplitude peaks (ridge points) in the
%      first row (referred as \omega_p(t)/2\pi in [3]), support lower
%      bounds (referred as \omega_-(t)/2/pi in [1]) - in the second row,
%      and the upper bounds (referred as \omega_+(t)/2/pi in [1]) - in
%      the third row.
% ecinfo: structure
%    - contains all the relevant information about the process of curve
%      extraction.
% Skel: 4x1 cell (returns empty matrix [] if 'Method' property is not 1,2,3 or 'nearest')
%    - contains the number of peaks N_p(t) in [Skel{1}], their frequency
%      indices m(t) in [Skel{2}], the corresponding frequencies
%      \nu_m(t)/2\pi in [Skel{3}], and the respective amplitudes Q_m(t)
%      in [Skel{4}] (in notations of [3]).
% 
% INPUT:
% TFR: NFxL matrix (rows correspond to frequencies, columns - to time)
%    - WFT or WT from which to extract [tfsupp]
% freq: NFx1 vector
%    - the frequencies corresponding to the rows of [TFR]
% wopt: structure | value (except if 'Method' property is 1 or 3) | 1x2 vector
%    - structure with parameters of the window/wavelet and the simulation,
%      returned as a third output by functions wft.m and wt.m;
%      alternatively, one can set wopt=[fs], where [fs] is the signal
%      sampling frequency (except methods 1 and 3); for methods 1 and 3
%      one can set wopt=[fs,D], where [D] is particular parameter of the
%      method (see [3]): (method 1) the characteristic growth rate of the
%      frequency - df/dt (in Hz/s) - for the WFT, or of the log-frequency
%      - d\log(f)/dt - for the WT; (method 2) the minimal distinguishable
%      frequency difference (in Hz) for WFT, or log-difference for WT.
%
% PROPERTIES: ({{...}} denotes default)
% ################################ BASIC ##################################
% 'Method': 1|{{2}}|3|'nearest'|'max'|efreq|1i*efreq
%    - method by which to extract the curve, as described in [3]; to use
%      frequency-based extraction, specify an L-length vector [efreq] of
%      the frequencies (in Hz), in which case the program will form the
%      ridge curve from the peaks lying at each time in the same support
%      (i.e. region of unimodal TFR amplitude) as the specified frequency
%      profile [efreq]; to select just ridge points nearest to [efreq]
%      (on a logarithmic scale for WT), use imaginary profile 1i*[efreq],
%      but this approach is more susceptible to noise and other effects.
% 'Param':[alpha,beta] - for method 2 (default = [1,1])
%                alpha - for method 1 (default = 1)
%                   [] - for all other methods (no parameters)
%    - parameters for each method, as described in [3].
% 'Normalize':{{'off'}}|'on'
%    - noise power in time-frequency domain can depend on frequency, so
%      that at lowest or highest frequencies the noise-induced amplitude
%      peaks might overgrow the peaks associated with a genuine components,
%      thus being selected instead of the latter. To avoid this, the curve
%      can be extracted using the normalized amplitude peaks, that are
%      non-uniformly reduced in dependence on their frequencies, with a
%      suitable normalization being determined based on the dependence of
%      the mean TFR amplitude on frequency; setting 'Normalize' to 'on'
%      applies such a normalization. Note that this property does not apply
%      if 'Method' property is set to [efreq] or [1i*efreq] (see above).
% 'Display': {{'on'}}|'off'|'notify'|'plot'
%    - to display the progress information or not; if set to 'notify',
%      displays only notifications if something went wrong; if set to
%      'plot', additionally plots all the obtained curves and their
%      characteristics (such as the averages and standard deviations of
%      ridge frequencies and their differences for scheme II iterations).
% 'Plot': 'on'|'on-wr'|{{'off'}}
%    - if set to 'on', plots additionally the full [TFR] and shows on it
%      the extracted frequency profile and the boundaries of its support at
%      each time (do not confuse with 'Display'='plot', which plots the
%      progress information and not only the final result). To avoid
%      unnecessary plotting of the huge data, by default the TFR plot is
%      resampled to have no more than few data points displayed per pixel;
%      to turn this option off, use addition '-wr' (i.e.\ 'Plot'='on-wr').
% ############################## ADVANCED #################################
% 'AmpFunc':function (default AmpFunc=@(x)log(x))
%    - the functional of the ridge amplitudes which to use in the
%      optimization, so that one maximizes (over all ridge sequences) the
%      path functional $\sum_n[AmpFunc(Q_p(t_n))+(penalization terms)]$,
%      where Q_p(t_n) denotes the ridge amplitudes at time t_n.
% 'PenalFunc': 1x2 cell with two functions {func1,func2}
%    - penalization functions used for curve extraction. For the WFT:
%      Method 1 - {@(drf,Xdrf),@(rf,Xrf)}, where [rf] and [drf] denote the
%      ridge frequency (in Hz) and its time derivative (difference between
%      two consecutive ridge frequencies times signal sampling frequency),
%      while [Xrf] and [Xdrf] are minimal resolvable frequency difference
%      and the characteristic rate of frequency change ($\Delta\xi_g/2\pi$
%      and $(\Delta\xi_g/2\pi)/\Delta\tau_g$ in notation of [1]); default =
%      = {@(drf,Xdrf)-p(1)*abs(drf./Xdrf),[]}, where p(1) is the parameter
%      of the method (see 'Param' property), while empty matrix [] means
%      no discrimination over [rf]).
%      Method 2 - {@(drf,m,s),@(rf,m,s)}, where [m] and [s] are (adaptively
%      chosen and refined) medians and quantiles of [drf] and [rf]; default
%      = {@(drf,m,s)-p(1)*abs((drf-m)./s),@(rf,m,s)-p(2)*abs((rf-m)./s)}).
%      For the WT all is the same but all frequency variables are taken
%      on a logarithmic scale, so that e.g. [rf] and [drf] now denote the
%      logarithms of the ridge frequencies and their time derivatives.
%      Note, that 'PenalFunc' property overrides the 'Param' property, so
%      the functions should be specified with all the parameters.
% 'PathOpt':{{'on'}}|'off'
%    - optimize the ridge curve over all possible trajectories ('on') or
%      use the one-step approach ('off'), see [3]; the path optimization
%      GREATLY improves the performance of all methods and is not
%      computationally expensive (is performed in O(N) operations using
%      fast algorithm of [3]), so DO NOT CHANGE THIS PROPERTY unless you
%      want to just play and see the advantages of the path optimization
%      over the one-step approach; note, that this property applies only
%      to methods 1 and 2, as the others are simple one-step approaches.
% 'Skel': {{[]}} | 4x1 cell returned as the third output of ecurve.m
%    - can be specified to avoid performing search of the peaks, their
%      numbers and corresponding amplitudes each time if the procedures
%      are applied to the same TFR (e.g. one compares the performance of
%      the different schemes).
% 'MaxIter': value (default = 20) /method 2 only/
%    - maximum number of iterations allowed for methods 2 and 3 to converge
%
%----------------------- Additional possibilities -------------------------
%
% NOTE: One can alternatively pass the structure with the properties as the
% 4th argument, e.g. /opt.Method=1; ecurve(TFR,freq,wopt,opt);/.
% If the other properties are specified next, they override those in the
% structure, e.g. /ecurve(TFR,freq,wopt,opt,'Method',2);/ will
% always use 2nd method, irrespectively to what is specified in [opt].
% 
%-------------------------------Examples-----------------------------------
%
% [WFT,freq,wopt]=wft(sig,fs); tfsupp=ecurve(WFT,freq,wopt);
% - extracts the ridge curve from the windowed Fourier transform [WFT] of
%   the signal [sig] sampled at [fs] Hz.
%
%------------------------------Changelog-----------------------------------
%
% v2.20:
% - some minor changes and fixes
% - changed slightly specifications of some properties
% - to make methods slightly more universal, changed penalization over both
%   frequency difference and frequency to the 1st order in both methods,
%   plus in method 2 changed mean/std to median/range
% - removed third scheme, as it always works worse than the other ones
% v2.10:
% - in methods 1 and 2, the frequency differences are now penalized to
%   the 2nd order (instead of the first order, as was before)
% - added possibility of normalization (see Properties)
% - some minor changes
% v2.00:
% - the algorithms are reprogrammed in accordance with the second version of [3]
% - some additional minor changes (Optimization/Display/Plotting)
% - added new possibilities (see PROPERTIES documentation)
% - added corrections for discretization effects, e.g. cubic
%   interpolation to better locate the peaks etc.
%
%--------------------------------------------------------------------------

function [tfsupp,varargout] = ecurve(TFR,freq,wopt,varargin)

[NF,L]=size(TFR); freq=freq(:);
tfsupp=zeros(3,L)*NaN; pind=zeros(1,L)*NaN; pamp=zeros(1,L)*NaN; idr=zeros(1,L)*NaN;

%Default parameters
method=2; pars=[]; NormMode='off'; DispMode='on'; PlotMode='off';
Skel=[]; PathOpt='on'; AmpFunc=@(x)log(x); PenalFunc={[],[]}; MaxIter=20;
%Update if user-defined
vst=1;
if nargin>3 && isstruct(varargin{1})
    copt=varargin{1}; vst=2;
    if isfield(copt,'Display'), cvv=copt.Display; if ~isempty(cvv), DispMode=cvv; end, end
    if isfield(copt,'Method'), cvv=copt.Method; if ~isempty(cvv), method=cvv; end, end
    if isfield(copt,'Param'), cvv=copt.Param; if ~isempty(cvv), pars=cvv; end, end
    if isfield(copt,'Normalize'), cvv=copt.Normalize; if ~isempty(cvv), NormMode=cvv; end, end
    if isfield(copt,'Plot'), cvv=copt.Plot; if ~isempty(cvv), PlotMode=cvv; end, end
    if isfield(copt,'Skel'), cvv=copt.Skel; if ~isempty(cvv), Skel=cvv; end, end
    if isfield(copt,'PathOpt'), cvv=copt.PathOpt; if ~isempty(cvv), PathOpt=cvv; end, end
    if isfield(copt,'AmpFunc'), cvv=copt.AmpFunc; if ~isempty(cvv), AmpFunc=cvv; end, end
    if isfield(copt,'PenalFunc'), cvv=copt.PenalFunc; if ~isempty(cvv), PenalFunc=cvv; end, end
    if isfield(copt,'MaxIter'), cvv=copt.MaxIter; if ~isempty(cvv), MaxIter=cvv; end, end
end
for vn=vst:2:nargin-3
    if strcmpi(varargin{vn},'Display'), if ~isempty(varargin{vn+1}), DispMode=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Method'), if ~isempty(varargin{vn+1}), method=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Param'), if ~isempty(varargin{vn+1}), pars=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Normalize'), if ~isempty(varargin{vn+1}), NormMode=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Plot'), if ~isempty(varargin{vn+1}), PlotMode=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'Skel'), if ~isempty(varargin{vn+1}), Skel=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'PathOpt'), if ~isempty(varargin{vn+1}), PathOpt=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'AmpFunc'), if ~isempty(varargin{vn+1}), AmpFunc=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'PenalFunc'), if ~isempty(varargin{vn+1}), PenalFunc=varargin{vn+1}; end
    elseif strcmpi(varargin{vn},'MaxIter'), if ~isempty(varargin{vn+1}), MaxIter=varargin{vn+1}; end
    else error(['There is no Property ''',varargin{vn},'''']);
    end
end
if isempty(pars) %if parameters are not specified, assign the defaults
    if ischar(method) || length(method)>1, pars=[];
    elseif method==1, pars=1;
    else pars=[1,1];
    end
else %if parameters are specified, check are they specified correctly
    if (ischar(method) || length(method)>1) && ~isempty(pars)
        error('Wrong number of parameters for the chosen method: there should be no parameters.');
    elseif method==1 && length(pars)~=1
        error('Wrong number of parameters for the chosen method: there should be one parameter.');
    elseif (method==2 || method==3) && length(pars)~=2
        error('Wrong number of parameters for the chosen method: there should be two parameters.');
    end
end
if nargout>1
    ec=struct;
    ec.Method=method; ec.Param=pars; ec.Display=DispMode; ec.Plot=PlotMode;
    ec.Skel=Skel; ec.PathOpt=PathOpt; ec.AmpFunc=AmpFunc;
end

%Determine the frequency resolution
if min(freq)<=0 || std(diff(freq))<std(diff(log(freq)))
    fres=1; fstep=mean(diff(freq));
    dfreq=[freq(1)-freq(end:-1:2);freq(end)-freq(end:-1:1)];    
else
    fres=2; fstep=mean(diff(log(freq)));
    dfreq=[log(freq(1))-log(freq(end:-1:2));log(freq(end))-log(freq(end:-1:1))];
end
%Assign numerical parameters
if isstruct(wopt)
    fs=wopt.fs; DT=(wopt.wp.t2h-wopt.wp.t1h);
    if fres==1, DF=(wopt.wp.xi2h-wopt.wp.xi1h)/2/pi;
    else DF=log(wopt.wp.xi2h/wopt.wp.xi1h); end
    if method==1, DD=DF/DT; end
    if method==3, DD=DF; end
else
    fs=wopt(1);
    if method==1 || method==3, DD=wopt(2); end
end

%//////////////////////////////////////////////////////////////////////////
TFR=abs(TFR); %convert to absolute values, since we need only them; also improves computational speed as TFR is no more complex and is positive
nfunc=ones(NF,1); tn1=find(~isnan(TFR(end,:)),1,'first'); tn2=find(~isnan(TFR(end,:)),1,'last'); sflag=0;
if (ischar(method) && ~strcmpi(method,'max')) || length(method)==1 %if not frequency-based or maximum-based extraction
    sflag=1;
    %----------------------------------------------------------------------
    %Construct matrices of ridge indices, frequencies and amplitudes:
    %[Ip],[Fp],[Qp], respectively; [Np] - number of peaks at each time.
    if ~isempty(Skel)
        Np=Skel{1}; Ip=Skel{2}; Fp=Skel{3}; Qp=Skel{4}; Mp=max(Np);
    else
        if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
            fprintf('Locating the amplitude peaks in TFR... ');
        end
        TFR=vertcat(zeros(1,L),TFR,zeros(1,L)); %pad TFR with zeros
        idft=1+find(TFR(2:end-1)>=TFR(1:end-2) & TFR(2:end-1)>TFR(3:end)); %find linear indices of the peaks
        [idf,idt]=ind2sub(size(TFR),idft); idf=idf-1; %find frequency and time indices of the peaks
        idb=find(idf==1 | idf==NF); idft(idb)=[]; idf(idb)=[]; idt(idb)=[]; %remove the border peaks
        dind=[0;find(diff(idt(:))>0);length(idt)]; Mp=max([max(diff(dind)),2]);
        Np=zeros(1,L); idn=zeros(length(idt),1);
        for dn=1:length(dind)-1,
            ii=dind(dn)+1:dind(dn+1); idn(ii)=1:length(ii);
            Np(idt(ii(1)))=length(ii);
        end
        idnt=sub2ind([Mp,L],idn(:),idt(:));
        %Quadratic interpolation to better locate the peaks
        a1=TFR(idft-1); a2=TFR(idft); a3=TFR(idft+1);
        dp=(1/2)*(a1-a3)./(a1-2*a2+a3);
        %Assign all
        Ip=ones(Mp,L)*NaN; Fp=ones(Mp,L)*NaN; Qp=ones(Mp,L)*NaN;
        Ip(idnt)=idf+dp; Qp(idnt)=a2-(1/4)*(a1-a3).*dp;
        if fres==1, Fp(idnt)=freq(idf)+dp(:)*fstep;
        else Fp(idnt)=freq(idf).*exp(dp(:)*fstep); end
        %Correct "bad" places, if present
        idb=find(isnan(dp) | abs(dp)>1 | idf==1 | idf==NF);
        if ~isempty(idb)
            Ip(idnt(idb))=idf(idb);
            Fp(idnt(idb))=freq(idf(idb));
            Qp(idnt(idb))=a2(idb);
        end
        %Remove zeros and clear the indices
        TFR=TFR(2:end-1,:); clear idft idf idt idn dind idnt a1 a2 a3 dp;
        %Display
        if ~strcmpi(DispMode,'off')
            if ~strcmpi(DispMode,'notify')
                fprintf('(number of ridges: %d+-%d, from %d to %d)\n',round(mean(Np(tn1:tn2))),round(std(Np(tn1:tn2))),min(Np),max(Np));
            end
            idb=find(Np(tn1:tn2)==0); NB=length(idb);
            if NB>0, fprintf(2,sprintf('Warning: At %d times there are no peaks (using border points instead).\n',NB)); end
        end
        %If there are no peaks, assign border points
        idb=find(Np(tn1:tn2)==0); idb=tn1-1+idb; NB=length(idb);
        if NB>0, G4=abs(TFR([1;2;NF-1;NF],idb)); end
        for bn=1:NB
            tn=idb(bn); cn=1; cg=G4(:,bn);
            if cg(1)>cg(2) || cg(4)>cg(3)
                if cg(1)>cg(2), Ip(cn,tn)=1; Qp(cn,tn)=cg(1); Fp(cn,tn)=freq(1); cn=cn+1; end
                if cg(4)>cg(3), Ip(cn,tn)=NF; Qp(cn,tn)=cg(4); Fp(cn,tn)=freq(NF); cn=cn+1; end
            else
                Ip(1:2,tn)=[1;NF]; Qp(1:2,tn)=[cg(1);cg(4)]; Fp(1:2,tn)=[freq(1);freq(NF)]; cn=cn+2;
            end
            Np(tn)=cn-1;
        end
        clear idb NB G4;
    end
    if nargout>2, varargout{2}={Np,Ip,Fp,Qp}; end
    if strcmpi(NormMode,'on'), nfunc=tfrnormalize(abs(TFR(:,tn1:tn2)),freq); end
    ci=Ip; ci(isnan(ci))=NF+2; cm=ci-floor(ci); ci=floor(ci); nfunc=[nfunc(1);nfunc(:);nfunc(end);NaN;NaN];
    Rp=(1-cm).*nfunc(ci+1)+cm.*nfunc(ci+2); Wp=AmpFunc(Qp.*Rp); nfunc=nfunc(2:end-3); %apply the functional to amplitude peaks
    
elseif ~ischar(method) && length(method)>1 %frequency-based extraction
    if length(method)~=L
        error('The specified frequency profile ("Method" property) should be of the same length as signal.');
    end
    
    efreq=method; submethod=1; if max(abs(imag(efreq)))>0, submethod=2; efreq=imag(efreq); end
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        if submethod==1, fprintf('Extracting the ridge curve lying in the same TFR supports as the specified frequency profile.\n');
        else fprintf('Extracting the ridge curve lying nearest to the specified frequency profile.\n'); end
    end
    
    tn1=max([tn1,find(~isnan(efreq),1,'first')]); tn2=min([tn2,find(~isnan(efreq),1,'last')]);
    if fres==1, eind=1+floor(0.5+(efreq-freq(1))/fstep);
    else eind=1+floor(0.5+log(efreq/freq(1))/fstep); end
    eind(eind<1)=1; eind(eind>NF)=NF;
    
    %Extract the indices of the peaks
    for tn=tn1:tn2
        cind=eind(tn); cs=abs(TFR(:,tn));
        
        %Ridge point
        cpeak=cind;
        if cind>1 && cind<NF
            if cs(cind+1)==cs(cind-1) || submethod==2
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
    
    %Transform to frequencies
    pind=tfsupp(1,:); tfsupp(:,tn1:tn2)=freq(tfsupp(:,tn1:tn2));
    pamp(tn1:tn2)=abs(TFR(sub2ind(size(TFR),pind(tn1:tn2),tn1:tn2)));
    
    %Optional output arguments
    if nargout>1
        ec.efreq=efreq; ec.eind=eind;
        ec.pfreq=tfsupp(1,:); ec.pind=pind; ec.pamp=pamp; ec.idr=idr;
        varargout{1}=ec;
    end
    if nargout>2, varargout{2}=[]; end
    
    %Plotting (if needed)
    if ~isempty(strfind(DispMode,'plot'))
        scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,2*scrsz(3)/3,2*scrsz(4)/3]);
        ax=axes('Position',[0.1,0.1,0.8,0.8],'FontSize',16); hold all;
        title(ax(1),'Ridge curve \omega_p(t)/2\pi'); ylabel(ax(1),'Frequency (Hz)'); xlabel(ax(1),'Time (s)');
        plot(ax(1),(0:L-1)/fs,efreq,'--','Color',[0.5,0.5,0.5],'LineWidth',2,'DisplayName','Specified frequency profile');
        plot(ax(1),(0:L-1)/fs,tfsupp(1,:),'-k','LineWidth',2,'DisplayName','Extracted frequency profile');
        legend(ax(1),'show'); if fres==2, set(ax(1),'YScale','log'); end
    end
    if ~isempty(strfind(PlotMode,'on')), plotfinal(tfsupp,TFR,freq,fs,DispMode,PlotMode); end
    if nargout>2, varargout{2}=Skel; end
    
    return;
end
%//////////////////////////////////////////////////////////////////////////


%--------------------------- Global Maximum -------------------------------
if strcmpi(method,'max') || length(pars)==2
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        if sflag==0, fprintf('Extracting the curve by Global Maximum scheme.\n');
        else fprintf('Extracting positions of global maximums (needed to estimate initial parameters).\n'); end
    end
    
    if sflag==0
        if strcmpi(NormMode,'on')
            nfunc=tfrnormalize(abs(TFR(:,tn1:tn2)),freq);
            TFR=TFR.*(nfunc(:)*ones(1,L));
        end
        for tn=tn1:tn2, [pamp(tn),pind(tn)]=max(abs(TFR(:,tn))); end
        tfsupp(1,tn1:tn2)=freq(pind(tn1:tn2));
        if strcmpi(NormMode,'on')
            TFR=TFR./(nfunc(:)*ones(1,L));
            pamp(tn1:tn2)=pamp(tn1:tn2)./(nfunc(pind(tn1:tn2))');
        end
    else
        for tn=tn1:tn2, [pamp(tn),idr(tn)]=max(Wp(1:Np(tn),tn)); end
        lid=sub2ind(size(Fp),idr(tn1:tn2),tn1:tn2);
        tfsupp(1,tn1:tn2)=Fp(lid); pind(tn1:tn2)=round(Ip(lid)); pamp(tn1:tn2)=Qp(lid);
    end
    idz=tn1-1+find(pamp(tn1:tn2)==0 | isnan(pamp(tn1:tn2)));
    if ~isempty(idz)
        idnz=tn1:tn2; idnz=idnz(~ismember(idnz,idz));
        pind(idz)=interp1(idnz,pind(idnz),idz,'linear','extrap');
        pind(idz)=round(pind(idz));
        tfsupp(1,idz)=interp1(idnz,tfsupp(1,idnz),idz,'linear','extrap');
    end
    
    if nargout>1, ec.pfreq=tfsupp(1,:); ec.pind=pind; ec.pamp=pamp; ec.idr=idr; end
    if ~isempty(strfind(DispMode,'plot')) && strcmpi(method,'max')
        scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,2*scrsz(3)/3,2*scrsz(4)/3]);
        ax=axes('Position',[0.1,0.1,0.8,0.8],'FontSize',16); hold all;
        title(ax(1),'Ridge curve \omega_p(t)/2\pi'); ylabel(ax(1),'Frequency (Hz)'); xlabel(ax(1),'Time (s)');
        plot(ax(1),(0:L-1)/fs,tfsupp(1,:),'-k','LineWidth',2,'DisplayName','Global Maximum curve');
        legend(ax(1),'show'); if fres==2, set(ax(1),'YScale','log'); end
    end
end


%------------------------- Nearest neighbour ------------------------------
if strcmpi(method,'nearest')
    %Display, if needed
    [~,imax]=max(Wp(:)); [fimax,timax]=ind2sub([Mp,L],imax); idr(timax)=fimax;
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        fprintf('Extracting the curve by Nearest Neighbour scheme.\n');
        fprintf('The highest peak was found at time %0.3f s and frequency %0.3f Hz (indices %d and %d, respectively).\n',...
            (timax-1)/fs,Fp(fimax,timax),timax,Ip(fimax,timax));
        fprintf('Tracing the curve forward and backward from point of maximum.\n');
    end
    %Main part
    for tn=timax+1:tn2, [~,idr(tn)]=min(abs(Ip(1:Np(tn),tn)-idr(tn-1))); end
    for tn=timax-1:-1:tn1, [~,idr(tn)]=min(abs(Ip(1:Np(tn),tn)-idr(tn+1))); end
    lid=sub2ind(size(Fp),idr(tn1:tn2),tn1:tn2);
    tfsupp(1,tn1:tn2)=Fp(lid); pind(tn1:tn2)=round(Ip(lid)); pamp(tn1:tn2)=Qp(lid);
    %Assign the output structure and display, if needed
    if nargout>1, ec.pfreq=tfsupp(1,:); ec.pind=pind; ec.pamp=pamp; ec.idr=idr; end
    if ~isempty(strfind(DispMode,'plot'))
        scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,2*scrsz(3)/3,2*scrsz(4)/3]);
        ax=axes('Position',[0.1,0.1,0.8,0.8],'FontSize',16); hold all;
        title(ax(1),'Ridge curve \omega_p(t)/2\pi'); ylabel(ax(1),'Frequency (Hz)'); xlabel(ax(1),'Time (s)');
        plot(ax(1),(0:L-1)/fs,tfsupp(1,:),'-k','LineWidth',2,'DisplayName','Nearest Neighbour curve');
        plot(ax(1),(timax-1)/fs,Fp(fimax,timax),'ob','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r',...
            'DisplayName','Starting point (overall maximum of TFR amplitude)');
        legend(ax(1),'show'); if fres==2, set(ax(1),'YScale','log'); end
    end
end


%----------------------------- Method I -----------------------------------
if length(pars)==1
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        fprintf('Extracting the curve by I scheme.\n');
    end
    %Define the functionals
    if isempty(PenalFunc{1}), logw1=@(x)-pars(1)*abs(fs*x/DD); else logw1=@(x)PenalFunc{1}(fs*x,DD); end
    if isempty(PenalFunc{2}), logw2=[]; else logw2=@(x)PenalFunc{2}(x,DF); end
    %Main part
    if strcmpi(PathOpt,'on'), idr=pathopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode);
    else [idr,timax,fimax]=onestepopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode); end
    lid=sub2ind(size(Fp),idr(tn1:tn2),tn1:tn2);
    tfsupp(1,tn1:tn2)=Fp(lid); pind(tn1:tn2)=round(Ip(lid)); pamp(tn1:tn2)=Qp(lid);
    %Assign the output structure and display, if needed
    if nargout>1, ec.pfreq=tfsupp(1,:); ec.pind=pind; ec.pamp=pamp; ec.idr=idr; end
    if ~isempty(strfind(DispMode,'plot'))
        scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,2*scrsz(3)/3,2*scrsz(4)/3]);
        ax=axes('Position',[0.1,0.1,0.8,0.8],'FontSize',16); hold all;
        title(ax(1),'Ridge curve \omega_p(t)/2\pi'); ylabel(ax(1),'Frequency (Hz)'); xlabel(ax(1),'Time (s)');
        plot(ax(1),(0:L-1)/fs,tfsupp(1,:),'-k','LineWidth',2,'DisplayName','Extracted frequency profile');
        if ~strcmpi(PathOpt,'on')
            plot(ax(1),(timax-1)/fs,Fp(fimax,timax),'ob','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r',...
                'DisplayName','Starting point (overall maximum of local functional)');
        end
        legend(ax(1),'show'); if fres==2, set(ax(1),'YScale','log'); end
    end
end


%----------------------------- Method II ----------------------------------
if length(pars)==2 && method==2
    %Initialize the parameters
    pf=tfsupp(1,tn1:tn2); if fres==2, pf=log(pf); end, dpf=diff(pf);
    mv=[median(dpf),0,median(pf),0];
    ss1=sort(dpf); CL=length(ss1); mv(2)=ss1(round(0.75*CL))-ss1(round(0.25*CL));
    ss2=sort(pf); CL=length(ss2); mv(4)=ss2(round(0.75*CL))-ss2(round(0.25*CL));
    
    %Display, if needed
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        if fres==1
            fprintf(['Maximums frequencies (median+-range): ']);
            fprintf('%0.3f+-%0.3f Hz; frequency differences: %0.3f+-%0.3f Hz.\n',mv(3),mv(4),mv(1),mv(2));
        else
            fprintf(['Maximums frequencies (log-median*/range ratio): ']);
            fprintf('%0.3f*/%0.3f Hz; frequency ratios: %0.3f*/%0.3f.\n',exp(mv(3)),exp(mv(4)),exp(mv(1)),exp(mv(2)));
        end
        fprintf('Extracting the curve by II scheme: iteration discrepancy - ');
        if ~isempty(strfind(DispMode,'plot'))
            scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,2*scrsz(3)/3,2*scrsz(4)/3]);
            ax=zeros(3,1);
            ax(1)=axes('Position',[0.1,0.6,0.8,0.3],'FontSize',16); hold all;
            if fres==2, set(ax(1),'YScale','log'); end
            ax(2)=axes('Position',[0.1,0.1,0.35,0.35],'FontSize',16); hold all;
            ax(3)=axes('Position',[0.55,0.1,0.35,0.35],'FontSize',16); hold all;
            title(ax(1),'Ridge curve \omega_p(t)/2\pi'); ylabel(ax(1),'Frequency (Hz)'); xlabel(ax(1),'Time (s)');
            ylabel(ax(2),'Frequency (Hz)'); ylabel(ax(3),'Frequency (Hz)'); xlabel(ax(3),'Iteration number'); xlabel(ax(2),'Iteration number');
            title(ax(2),'${\rm m}[d\omega_p/dt]/2\pi$ (solid), ${\rm s}[d\omega_p/dt]/2\pi$ (dashed)','interpreter','Latex','FontSize',20);
            title(ax(3),'${\rm m}[\omega_p]/2\pi$ (solid), ${\rm s}[\omega_p]/2\pi$ (dashed)','interpreter','Latex','FontSize',20);
            line0=plot(ax(1),(0:L-1)/fs,tfsupp(1,:),':','Color',[0.5,0.5,0.5],'DisplayName','Global Maximum ridges');
            line1=plot(ax(2),0,fs*mv(1),'-sk','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k','DisplayName','m[d\omega_p/dt]/2\pi');
            line2=plot(ax(2),0,fs*mv(2),'--ok','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k','DisplayName','s[d\omega_p/dt]/2\pi');
            line3=plot(ax(3),0,mv(3),'-sk','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k','DisplayName','m[\omega_p]/2\pi');
            line4=plot(ax(3),0,mv(4),'--ok','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','k','DisplayName','s[\omega_p]/2\pi');
        end
    end
    
    %Iterate
    rdiff=NaN; itn=0; allpind=zeros(10,L); allpind(1,:)=pind;
    if nargout>1, ec.mv=mv; ec.rdiff=rdiff; end
    while rdiff~=0
        %Define the functionals
        smv=[mv(2),mv(4)]; %to avoid underflow
        if smv(1)<=0, smv(1)=10^(-32)/fs; end
        if smv(2)<=0, smv(2)=10^(-16); end
        if isempty(PenalFunc{1}), logw1=@(x)-pars(1)*abs((x-mv(1))/smv(1)); else logw1=@(x)PenalFunc{1}(x,mv(1),smv(1)); end
        if isempty(PenalFunc{2}), logw2=@(x)-pars(2)*abs((x-mv(3))/smv(2)); else logw2=@(x)PenalFunc{2}(x,mv(3),smv(2)); end
        %Calculate all
        pind0=pind;
        if strcmpi(PathOpt,'on'), idr=pathopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode);
        else [idr,timax,fimax]=onestepopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode); end
        lid=sub2ind(size(Fp),idr(tn1:tn2),tn1:tn2);
        tfsupp(1,tn1:tn2)=Fp(lid); pind(tn1:tn2)=round(Ip(lid)); pamp(tn1:tn2)=Qp(lid);
        rdiff=length(find(pind(tn1:tn2)-pind0(tn1:tn2)~=0))/(tn2-tn1+1);
        itn=itn+1;
        %Update the medians/ranges
        pf=tfsupp(1,tn1:tn2); if fres==2, pf=log(pf); end, dpf=diff(pf);
        mv=[median(dpf),0,median(pf),0];
        ss1=sort(dpf); CL=length(ss1); mv(2)=ss1(round(0.75*CL))-ss1(round(0.25*CL));
        ss2=sort(pf); CL=length(ss2); mv(4)=ss2(round(0.75*CL))-ss2(round(0.25*CL));
        %Update the information structure, if needed
        if nargout>1
            ec.pfreq=[ec.pfreq;tfsupp(1,:)]; ec.pind=[ec.pind;pind]; ec.pamp=[ec.pamp;pamp]; ec.idr=[ec.idr;idr];
            ec.mv=[ec.mv;mv]; ec.rdiff=[ec.rdiff,rdiff];
        end
        %Display, if needed
        if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
            fprintf('%0.2f%%; ',100*rdiff);
            if ~isempty(strfind(DispMode,'plot'))
                line0=plot(ax(1),(0:L-1)/fs,tfsupp(1,:),'DisplayName',sprintf('Iteration %d (discrepancy %0.2f%%)',itn,100*rdiff));
                set(line1,'XData',0:itn,'YData',[get(line1,'YData'),fs*mv(1)]);
                set(line2,'XData',0:itn,'YData',[get(line2,'YData'),fs*mv(2)]);
                set(line3,'XData',0:itn,'YData',[get(line3,'YData'),mv(3)]);
                set(line4,'XData',0:itn,'YData',[get(line4,'YData'),mv(4)]);
                if ~strcmpi(PathOpt,'on'),
                    mpt=plot(ax(1),(timax-1)/fs,Fp(fimax,timax),'ok','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','w',...
                        'DisplayName',['Starting point (iteration ',num2str(itn),')']);
                end
            end
        end
        %Stop if maximum number of iterations has been reached
        if itn>MaxIter && rdiff~=0
            if ~strcmpi(DispMode,'off')
                if ~strcmpi(DispMode,'notify'), fprintf('\n'); end
                fprintf('WARNING! Did not fully converge in %d iterations (current ''MaxIter''). Using the last estimate.',MaxIter);
            end
            break;
        end
        %Just in case, check for ``cycling'' (does not seem to occur in practice)
        allpind(itn+1,:)=pind; gg=Inf;
        if rdiff~=0 && itn>2
            for kn=2:itn-1, gg=min([gg,length(find(pind(tn1:tn2)-allpind(kn,tn1:tn2)~=0))]); end
        end
        if gg==0
            if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
                fprintf('converged to a cycle, terminating iteration.');
            end
            break;
        end
        
    end
    
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify'), fprintf('\n'); end
    if ~isempty(strfind(DispMode,'plot'))
        set(line0,'Color','k','LineWidth',2);
        if ~strcmpi(PathOpt,'on'), set(mpt,'Color',b,'MarkerFaceColor','k'); end
        if fres==2 %change plot if the resolution is logarithmic
            set(line1,'YData',exp(get(line1,'YData')),'DisplayName','exp(m[d\log\omega_p/dt])');
            set(line2,'YData',exp(get(line2,'YData'))-1,'DisplayName','exp(s[d\log\omega_p/dt])-1');
            set(line3,'YData',exp(get(line3,'YData')),'DisplayName','exp(m[\log\omega_p])/2\pi');
            set(line4,'YData',exp(get(line4,'YData'))-1,'DisplayName','exp(s[\log\omega_p])-1');
            set(ax(2:3),'YScale','log'); ylabel(ax(2),'Frequency Ratio'); ylabel(ax(3),'Frequency (Hz)');
            title(ax(2),'$e^{{\rm m}[d\log\omega_p/dt]}$ (solid), $e^{{\rm s}[d\log\omega_p/dt]}-1$ (dashed)','interpreter','Latex','FontSize',20);
            title(ax(3),'$e^{{\rm m}[\log\omega_p]}/2\pi$ (solid), $e^{{\rm s}[\log\omega_p]}-1$ (dashed)','interpreter','Latex','FontSize',20);
            iline=[get(line1,'YData'),get(line2,'YData')]; set(ax(2),'YLim',[0.75*min(iline),1.5*max(iline)]);
            iline=[get(line3,'YData'),get(line4,'YData')]; set(ax(3),'YLim',[0.75*min(iline),1.5*max(iline)]);
        end
    end
    
end
if nargout>1, varargout{1}=ec; end

%//////////////////////////////////////////////////////////////////////////
%Extract the time-frequency support around the ridge points
for tn=tn1:tn2
    cs=abs(TFR(:,tn)); cpeak=pind(tn);
    iup=[]; idown=[];
    if cpeak<NF-1, iup=cpeak+find(cs(cpeak+1:end-1)<=cs(cpeak:end-2) & cs(cpeak+1:end-1)<cs(cpeak+2:end),1,'first'); end
    if cpeak>2, idown=cpeak-find(cs(cpeak-1:-1:2)<=cs(cpeak:-1:3) & cs(cpeak-1:-1:2)<cs(cpeak-2:-1:1),1,'first'); end
    iup=min([iup,NF]); idown=max([idown,1]);
    tfsupp(2,tn)=idown; tfsupp(3,tn)=iup;
end
tfsupp(2:3,tn1:tn2)=freq(tfsupp(2:3,tn1:tn2));
%//////////////////////////////////////////////////////////////////////////

%Final display
if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
    fprintf('Curve extracted: ridge frequency %0.2f+-%0.2f Hz, lower support %0.2f+-%0.2f Hz, upper support %0.2f+-%0.2f Hz.\n',...
        mean(tfsupp(1,:)),std(tfsupp(1,:)),mean(tfsupp(2,:)),std(tfsupp(2,:)),mean(tfsupp(3,:)),std(tfsupp(3,:)));
end

%Plot (if needed)
if ~isempty(strfind(PlotMode,'on')), plotfinal(tfsupp,TFR,freq,fs,DispMode,PlotMode,nfunc); end

end

%==========================================================================
%========================= Support functions ==============================
%==========================================================================

%==================== Path optimization algorithm =========================
function idid=pathopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode)

[Mp,L]=size(Fp); NF=length(freq);
tn1=find(Np>0,1,'first'); tn2=find(Np>0,1,'last');
if min(freq)>0 && std(diff(freq))>std(diff(log(freq))), Fp=log(Fp); end

%Weighting functions
if ~isa(logw1,'function_handle')
    if isempty(logw1)
        logw1=zeros(2*NF+1,L);
    else
        logw1=[2*logw1(1)-logw1(2);logw1(:);2*logw1(end)-logw1(end-1)];
    end
end
if ~isa(logw2,'function_handle')
    if isempty(logw2)
        W2=zeros(Mp,L);
    else
        logw2=[2*logw2(1)-logw2(2);logw2(:);2*logw2(end)-logw2(end-1);NaN;NaN];
        ci=Ip; ci(isnan(ci))=NF+2; cm=ci-floor(ci); ci=floor(ci);
        W2=(1-cm).*logw2(ci+1)+cm.*logw2(ci+2); clear ci cm;
    end
else
    W2=logw2(Fp);
end
W2=Wp+W2;

%The algorithm by itself
q=zeros(Mp,L)*NaN; U=zeros(Mp,L)*NaN;
q(1:Np(tn1),tn1)=0; U(1:Np(tn1),tn1)=W2(1:Np(tn1),tn1);
if isa(logw1,'function_handle')
    for tn=tn1+1:tn2
        cf=Fp(1:Np(tn),tn)*ones(1,Np(tn-1))-ones(Np(tn),1)*Fp(1:Np(tn-1),tn-1)'; CW1=logw1(cf);
        [U(1:Np(tn),tn),q(1:Np(tn),tn)]=max(W2(1:Np(tn),tn)*ones(1,Np(tn-1))+CW1+ones(Np(tn),1)*U(1:Np(tn-1),tn-1)',[],2);
    end
else
    for tn=tn1+1:tn2
        ci=Ip(1:Np(tn),tn)*ones(1,Np(tn-1))-ones(Np(tn),1)*Ip(1:Np(tn-1),tn-1)';
        ci=ci+NF; cm=ci-floor(ci); ci=floor(ci);
        if Np(tn)>1, CW1=(1-cm).*logw1(ci+1)+cm.*logw1(ci+2);
        else CW1=(1-cm).*logw1(ci+1)'+cm.*logw1(ci+2)'; end
        [U(1:Np(tn),tn),q(1:Np(tn),tn)]=max(W2(1:Np(tn),tn)*ones(1,Np(tn-1))+CW1+ones(Np(tn),1)*U(1:Np(tn-1),tn-1)',[],2);
    end
end

%Recover the indices
idid=zeros(1,L)*NaN; [~,idid(tn2)]=max(U(:,tn2));
for tn=tn2-1:-1:tn1, idid(tn)=q(idid(tn+1),tn+1); end

%Plot if needed
%{
if ~isempty(strfind(DispMode,'plot+'))
    figure; axes('FontSize',16,'Box','on'); hold all;
    rlines=zeros(Np(tn2),1);
    for pn=1:Np(tn2)
        cridges=zeros(L,1)*NaN; cridges(tn2)=pn;
        for tn=tn2-1:-1:tn1
            cridges(tn)=q(cridges(tn+1),tn+1);
        end
        lind=sub2ind(size(Ip),cridges(tn1:tn2),(tn1:tn2)');
        cridges=[ones(tn1-1,1)*NaN;Ip(lind);ones(L-tn2,1)*NaN];
        crfreq=[ones(tn1-1,1);freq(cridges(tn1:tn2));ones(L-tn2,1)];
        rlines(pn)=plot((0:L-1)/fs,crfreq,'DisplayName',['U=',num2str(U(pn,tn2))]);
    end
    [~,imax]=max(U(:,tn2)); uistack(rlines(imax),'top');
    set(rlines(imax),'Color','k','LineWidth',2);
    title('Path to each of the last peaks maximizing the given functional');
    ylabel('Frequency (Hz)'); xlabel('Time (s)');
end
%}

end

%================== One-step optimization algorithm =======================
function [idid,varargout]=onestepopt(Np,Ip,Fp,Wp,logw1,logw2,freq,DispMode)

[Mp,L]=size(Fp); NF=length(freq);
tn1=find(Np>0,1,'first'); tn2=find(Np>0,1,'last');
if min(freq)>0 && std(diff(freq))>std(diff(log(freq))), Fp=log(Fp); end

%Weighting functions
if ~isa(logw1,'function_handle')
    if isempty(logw1)
        logw1=zeros(2*NF+1,L);
    else
        logw1=[2*logw1(1)-logw1(2);logw1(:);2*logw1(end)-logw1(end-1)];
    end
end
if ~isa(logw2,'function_handle')
    if isempty(logw2)
        W2=zeros(Mp,L);
    else
        logw2=[2*logw2(1)-logw2(2);logw2(:);2*logw2(end)-logw2(end-1);NaN;NaN];
        ci=Ip; ci(isnan(ci))=NF+2; cm=ci-floor(ci); ci=floor(ci);
        W2=(1-cm).*logw2(ci+1)+cm.*logw2(ci+2); clear ci cm;
    end
else
    W2=logw2(Fp);
end
W2=Wp+W2;

%The algorithm by itself
[~,imax]=max(W2(:)); [fimax,timax]=ind2sub([Mp,L],imax); %determine the starting point
idid=zeros(1,L)*NaN; idid(timax)=fimax;
if isa(logw1,'function_handle')
    for tn=timax+1:tn2
        cf=Fp(1:Np(tn),tn)-Fp(idid(tn-1),tn-1);
        CW1=logw1(cf); [~,idid(tn)]=max(W2(1:Np(tn),tn)+CW1);
    end
    for tn=timax-1:-1:tn1
        cf=-Fp(1:Np(tn),tn)+Fp(idid(tn+1),tn+1);
        CW1=logw1(cf); [~,idid(tn)]=max(W2(1:Np(tn),tn)+CW1);
    end
else
    for tn=timax+1:tn2
        ci=NF+Ip(1:Np(tn),tn)-Ip(idid(tn-1),tn-1); cm=ci-floor(ci); ci=floor(ci);
        CW1=(1-cm).*logw1(ci+1)+cm.*logw1(ci+2); [~,idid(tn)]=max(W2(1:Np(tn),tn)+CW1);
    end
    for tn=timax-1:-1:tn1
        ci=NF-Ip(1:Np(tn),tn)+Ip(idid(tn+1),tn+1); cm=ci-floor(ci); ci=floor(ci);
        CW1=(1-cm).*logw1(ci+1)+cm.*logw1(ci+2); [~,idid(tn)]=max(W2(1:Np(tn),tn)+CW1);
    end
end

if nargout>1, varargout{1}=timax; end
if nargout>2, varargout{2}=fimax; end

end

%========================= Plotting function ==============================
function plotfinal(tfsupp,TFR,freq,fs,DispMode,PlotMode,varargin)

[NF,L]=size(TFR); fres=1; if min(freq)>0 && std(diff(freq))>std(diff(log(freq))), fres=2; end
nfunc=ones(NF,1); if nargin>6 && ~isempty(varargin{1}), nfunc=varargin{1}; end
XX=(0:(L-1))/fs; YY=freq; ZZ=abs(TFR).*(nfunc(:)*ones(1,L)); scrsz=get(0,'ScreenSize');

MYL=round(scrsz(3)); MXL=round(scrsz(4)); %maximum number of points seen in plots
if isempty(strfind(lower(PlotMode),'wr')) && (size(ZZ,1)>MYL || size(ZZ,2)>MXL)
    if ~strcmpi(DispMode,'off') && ~strcmpi(DispMode,'notify')
        fprintf('Plotting: TFR contains more data points (%d x %d) than pixels in the plot, so for a\n',size(ZZ,1),size(ZZ,2));
        fprintf('          better performance its resampled version (%d x %d) will be displayed instead.\n',min([MYL,size(ZZ,1)]),min([MXL,size(ZZ,2)]));
    end
    XI=XX; YI=YY;
    if size(ZZ,2)>MXL, XI=linspace(XX(1),XX(end),MXL); end
    if fres==1
        if size(ZZ,1)>MYL, YI=linspace(YY(1),YY(end),MYL); end
        ZZ=aminterp(XX,YY,ZZ,XI,YI,'max');
    else
        if size(ZZ,1)>MYL, YI=exp(linspace(log(YY(1)),log(YY(end)),MYL)); end
        ZZ=aminterp(XX,log(YY),ZZ,XI,log(YI),'max');
    end
    XX=XI(:); YY=YI(:);
end

figure('Position',[scrsz(3)/6,scrsz(4)/4,2*scrsz(3)/3,2*scrsz(4)/3]);
ax0=axes('Position',[0.1,0.15,0.8,0.7],'NextPlot','Add','FontSize',18,...
    'XLim',[XX(1),XX(end)],'YLim',[YY(1),YY(end)],'Layer','top','Box','on');
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title({'TFR amplitude', '(black: extracted peaks; gray: their supports)'});
if fres==2, set(gca,'YScale','log'); end
pc=pcolor(XX,YY,ZZ); set(pc,'LineStyle','none');
plot((0:(L-1))/fs,tfsupp(1,:),'Color','k','LineWidth',2);
plot((0:(L-1))/fs,tfsupp(2,:),'Color',[0.5,0.5,0.5],'LineWidth',2);
plot((0:(L-1))/fs,tfsupp(3,:),'Color',[0.5,0.5,0.5],'LineWidth',2);
if ~isempty(find(nfunc~=1,1))
    set(ax0,'Position',[0.1,0.15,0.65,0.7]);
    title(ax0,{'Normalized TFR amplitude', '(black: extracted peaks; gray: their supports)'});
    axes('Position',[0.8,0.15,0.175,0.7],'FontSize',18,'XLim',[0,1.25],'YLim',[freq(1),freq(end)],'Box','on');
    hold all; plot(nfunc,freq,'-k','LineWidth',2); title({'Normalization','function'});
    set(gca,'XTick',[0,0.5,1],'YTickLabel',{}); if fres==2, set(gca,'YScale','log'); end
end

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
%--------------------------------------------------------------------------
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

%================= Normalization of the noise peaks =======================
function nfunc = tfrnormalize(TFR,freq)

%Calculate medians and percentiles
NF=length(freq); L=size(TFR,2); TFR=sort(TFR,2); mm=median(TFR,2);
pp=0.75; ss=TFR(:,round((0.5+pp/2)*L))-TFR(:,round((0.5-pp/2)*L));

%Calculate weightings
gg=ss./mm; gg(isnan(gg))=Inf; ii=find(isfinite(gg)); gg=gg-median(gg(ii));
zz=sort(abs(gg(ii))); zz=zz(round(0.25*length(zz))); rw=exp(-abs(gg)/zz/2);
rw=ones(NF,1);

%Fitting
Y=mm; ii=find(freq>0 & Y>0); CN=length(ii);
Y=rw(ii).*log(Y(ii)); X=log(freq(ii)); FM=(rw(ii)*ones(1,2)).*[ones(CN,1),X(:)];
b=pinv(FM)*Y(:);

%Construct the normalization function
nfunc=freq(freq>0).^b(2); nfunc=[nfunc(1)*ones(length(freq(freq<=0)),1);nfunc(:)];
dd=mm(:)-exp(b(1))*nfunc; pid=find(dd>0 & freq>0); [~,rid]=min(abs(cumsum(dd(pid))-0.5*sum(dd(pid)))); mff=freq(pid(rid));
nfunc=(mff^b(2))./nfunc; nfunc(nfunc>1 | isnan(nfunc))=1;
nfunc=1+(nfunc(:)-1).*(rw(:).^2);

end