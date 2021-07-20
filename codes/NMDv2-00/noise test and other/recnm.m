%================= Reconstruct the full Nonlinear Mode ====================
% Version 1.00 stable
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
% [NM,Optional:ah,ph,rhamp,rhphi,rhfreq]=recnm(hid,hamp,hphi,hfreq,Optional:ch,mode)
% - given the numbers of extracted harmonics [hid], their amplitudes [hamp],
%   phases [hphi] and frequencies [hfreq], returns the reconstructed
%   Nonlinear Mode [NM]; by optional input [ch] ( default =
%   ones(length(hid)),1) ) one can specify the expected absolute errors,
%   up to the constant multiplier, for each harmonic (e.g. based on the 
%   noise levels  in the corresponding frequency ranges), which are then
%   used in the refinement process to better recover the full [NM].
%   .......................................................................
%   If requested, returns also the ratios [ah] of the amplitudes of each
%   harmonic to the first one and the corresponding phase shifts [ph], so
%   that e.g. one can recover the waveform as
%   wff=@(phi)sum((ah(:)*ones(1,length(phi))).*cos(hid(:)*phi(:)'+ph(:)*ones(1,length(phi))),1)
%   (then it can be visualized by pp=linspace(-pi,pi,1000);
%   plot(pp,wff(pp));); also if requested returns the refined estimates
%   of the harmonic amplitudes [rhamp], phases [rhphi] and frequencies
%   [rhfreq]; the amplitude, phase and frequency of the full NM is the
%   amplitude, phase and frequency of its first harmonic rhamp(1,:),
%   rhphi(1,:) and rhfreq(1,:), respectively.
%
% mode:'mean'|'median'
% - specifies what to use, the mean or median values of the amplitudes and
%   phase-differences for determining the amplitude ratios and phase shifts
%   between harmonics.
%
%-------------------------------Examples-----------------------------------
%
% /Calculate first WFT of the reference signal [rsig] (see wft function):/
% [WFT,freq,wopt]=wft(rsig,fs);
% /Extract somehow the component time-frequency support from it, e.g./
% tfsupp=ecurve(WFT,freq);
% /Use it to extract all associated harmonics and their characteristics:/
% [hid,hamp,hphi,hfreq]=eharm(sig,tfsupp,TFR,freq,wopt);
% /Reconstruct the full Nonlinear Mode [NM]:/
% NM=recnm(hid,hamp,hphi,hfreq);
% /For WT the same procedure is used./
%
%--------------------------------------------------------------------------

function [NM,varargout] = recnm(hid,hamp,hphi,hfreq,varargin)

[NH,L]=size(hamp); for hn=1:NH, hphi(hn,:)=unwrap(hphi(hn,:)); end
ch=ones(NH,1); if nargin>4 && ~isempty(varargin{1}), ch=varargin{1}(:); end
method='mean'; if nargin>5 && ~isempty(varargin{2}), method=varargin{2}; end

for hn=1:NH, hphi(1,:)=unwrap(hphi(1,:)); end
ah=zeros(NH,1); ph=zeros(NH,NH);
rhamp=zeros(size(hamp)); rhphi=zeros(size(hphi)); rhfreq=zeros(size(hfreq));

%--------------------------------------------------------------------------
if strcmpi(method,'mean')
    for hn=1:NH
        ah(hn)=mean(hamp(hn,:));
        for kn=hn+1:NH
            ph(kn,hn)=angle(mean(exp(1i*(hid(hn)*hphi(kn,:)-hid(kn)*hphi(hn,:)))));
        end
    end
else
    for hn=1:NH
        ah(hn)=median(hamp(hn,:));
        for kn=hn+1:NH
            ph(kn,hn)=angle(median(exp(1i*(hid(hn)*hphi(kn,:)-hid(kn)*hphi(hn,:)))));
        end
    end
end
ah=ah/ah(1); ph(ph(:)>pi)=ph(ph(:)>pi)-2*pi; ph=ph-ph';
%--------------------------------------------------------------------------

for hn=1:NH, rhamp(1,:)=rhamp(1,:)+hamp(hn,:)/ch(hn); end, rhamp(1,:)=rhamp(1,:)/sum(ah./ch);
for hn=1:NH
    rhamp(hn,:)=ah(hn)*rhamp(1,:); cnr=0;
    for kn=1:NH
        cc=min([1,hid(kn)/hid(hn)])*(ah(kn)/ch(kn)); cnr=cnr+cc;
        rhphi(hn,:)=rhphi(hn,:)+cc*exp(1i*(hid(hn)*hphi(kn,:)-ph(kn,hn)-2*pi*round((hid(hn)*hphi(kn,:)-hid(kn)*hphi(hn,:)-ph(kn,hn))/2/pi))/hid(kn));
        %rhphi(hn,:)=rhphi(hn,:)+cc*(hid(hn)*hphi(kn,:)-ph(kn,hn)-2*pi*round((hid(hn)*hphi(kn,:)-hid(kn)*hphi(hn,:)-ph(kn,hn))/2/pi))/hid(kn);
        rhfreq(hn,:)=rhfreq(hn,:)+cc*hfreq(kn,:)*(hid(hn)/hid(kn));
    end
    rhphi(hn,:)=angle(rhphi(hn,:));
    %rhphi(hn,:)=rhphi(hn,:)/cnr;
    rhfreq(hn,:)=rhfreq(hn,:)/cnr;
end

%--------------------------------------------------------------------------

NM=sum(rhamp.*cos(rhphi),1);
if nargout>1, varargout{1}=ah; end
if nargout>2, varargout{2}=ph(:,1); end
if nargout>3, varargout{3}=rhamp; end
if nargout>4, varargout{4}=rhphi; end
if nargout>5, varargout{5}=rhfreq; end

end

