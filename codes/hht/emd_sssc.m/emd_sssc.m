function [imf,ort,fvs,iterNum] = emd_sssc(x,fs,varargin)
%
% ERPHM Code Log
% 2017-06-20 Created by Zhiliang Liu, Dandan Peng and Yaqiang Jin.
% 2019-01-19 Modified by Zhiliang Liu and Dandan Peng.
%
% This funcion performs the empirical mode decomposition(EMD) with a soft sifting stopping criterion (SSSC) on the input signal, and
% return intrinsic mode function (IMF). If you have any questions, please contact us via Zhiliang_Liu@uestc.edu.cn
%
% SYNTAX:
% imf = emd_sssc(x);
% [imf,ort,fvs,iterNum] = emd_sssc(x,fs);
% [imf,ort,fvs,iterNum] = emd_sssc(x,fs,options);
% [imf,ort,fvs,iterNum] = emd_sssc(x,fs,'option_name1',option_value1,....);
%
% INPUTS:
% [1] x: a row signal vecter.
% [2] option: a struct contains options' name(option_name) and
%         corresponding values(option_value).
% [2.1] 'display'
%       plot imf or not,1 plot,0 do not.
%       default: 0
% [2.2] 'max_iter'
%       max number of iterations in one imf sift procedure.
%       default: 30
% [2.3] 'max_imfs'
%   	max number of imf obtained in emd procedure.
%   	default: 10
% [2.4] 'ext_ratio'
%       end extension length to original data length ratio.
%       can be any positive real number from (0,1]
%   	default: 0.2
%
% OUTPUTS
% [1] imf
%       Intrinsic Mode Function (imf) Matrix of which each row is an imf and the last
%       row is the residual signal.
% [2] iterNum
%       iteration number for each imf.
% [3] fvs
%       funciton values.
% [4] ort
%       index of orthogonality.
%
% REFERENCE
% [1] Dandan Peng, Zhiliang Liu, Yaqiang Jin, Yong Qin. Improved EMD with a Soft
%     Sifting Stopping Criterion and Its Application to Fault Diagnosis of Rotating Machinery.
%     Journal of Mechanical Engineering. Accepted on Jan. 01, 2019.
% [2] Zhiliang Liu, Yaqiang Jin, Ming J. Zuo, and Zhipeng Feng. Time-frequency
%     representation based on robust local mean decomposition for multi-component
%     AM-FM signal analysis. Mechanical Systems and Signal Processing. 95: 468-487, 2017.
%
% EXAMPLE:
% clc; clear; close all;
% fs = 10000; % sampling frequency
% N = 30000; % data amount
% t = (1:N)/fs; % time vector
% x = (2+cos(2*pi*0.5*t)).*cos(2*pi*5*t+15*t.^2)+...
%      cos(2*pi*2*t);
% [imf,ort,fvs,iterNum] = emd_sssc(x,fs,'display',1);

[x, display, ssc, max_iter, max_imfs, smooth_mode,ext_ratio, x_energy, imf, iterNum, fvs]= initial(x,varargin{:});

% Initialize main loop
i = 0;
xs = x; % copy x to xs for sifting process, reserve original input as x.
nx = length(x);

while   i < max_imfs && ~stopemd(xs, x_energy) % outer loop for imf selection
    
    i = i+1;
    % initialize variables used in imf sifting loop
    s_j = zeros(max_iter,nx);
    
    % imf sifting iteration loop
    j = 0;
    stop_flag_sifting_process = 0;
    s = xs;
    while  j < max_iter && ~stop_flag_sifting_process %  inner loop for sifting process
        
        j = j+1;
        if j==1
            [m_j, n_extr] = emd_mean(s, smooth_mode, ext_ratio);
        end
        
        % force to stop iteration if number of extrema of s is smaller than 3.
        if n_extr < 3
            break;
        end
        s = s-m_j; % remove mean.
        s_j(j, :) = s;
        [m_j, n_extr] = emd_mean(s, smooth_mode,ext_ratio);
        [stop_flag_sifting_process,fvs(i,:)] = is_sifting_process_stop(m_j, s,j, fvs(i,:), ssc);
    end
    
    switch ssc
        case {'liu'}
            [~, opt0] = min(fvs(i,1:j)); % ***Critical Step***
            opt_IterNum = min(j, opt0); % in case iteration stop for n_extr<3
    end
    
    imf(i, :) =s_j(opt_IterNum, :); % gain intrinsic mode function
    xs = xs-imf(i, :); % remove imf just obtained from input signal
    iterNum(i) = opt_IterNum; % record the iteration number taken by each imf sifing
end

imf(i+1, :) = xs; % save residual in the last row of imf matrix.
imf(i+2:end,:) = [];
fvs(i+1:end,:) = [];
ort = io(x, imf);

% Output visualization
if display == 1
    emdplot(x,imf,fs);
end

end

% initialization
function [x, display, ssc, max_iter, max_imfs, smooth_mode, ext_ratio, x_energy, imf, iterNum, fvs] = initial(x,varargin)

% option fields(i.e. name)
optn_fields = {'display',  'ssc',  'max_iter',...
    'max_imfs', 'smooth_mode','ext_ratio', 'fix','fix_h'};

% set default options(def_opts)
def_optns.display = 0; % do not plot imf.
def_optns.ssc = 'liu'; % sifting stopping criterion
def_optns.max_iter = 30; % max iteration number in a IMF sifting process
def_optns.max_imfs = 10; % max number of imf
def_optns.smooth_mode = 'spline'; % cubic - moving average, spline - pchip
def_optns.ext_ratio = 0.2; % end extension length to original data

optns = def_optns; % opts stores the final options.

% get user input options(in_opts)
if nargin == 1 % use default options(see above).
    in_optns = def_optns;
elseif nargin == 2 && isstruct(varargin{1})
    % 1st argument is x, 2nd is options in a struct.
    in_optns = varargin{1};
elseif nargin > 2 % input options seperately.
    try
        in_optns = struct(varargin{:});
    catch
        error('wrong argmument syntax')
    end
else
    error('arguments error: maybe not enough or wrong syntax')
end

names = fieldnames(in_optns);% get input options' name and value

for k = names'
    if ~any(strcmpi(char(k), optn_fields))
        % find any wrong argument in syntax.
        error(['bad option field name: ',char(k)])
    end
    if ~isempty(eval(['in_optns.',char(k)]))
        % alter default option values with input, and empty input keep default.
        eval(['optns.',lower(char(k)),' = in_optns.',char(k)])
    end
end

display = optns.display;
ssc = optns.ssc;
max_iter = optns.max_iter;
max_imfs = optns.max_imfs;
smooth_mode = optns.smooth_mode;
ext_ratio = optns.ext_ratio;

% initialize x(input signal), x_energy, IMF, iterNum and fvs 
x = x(:)'; % make x a row vector.
nx = length(x);
x_energy = sum(x.^2); % energy = square summation.
imf = zeros(max_imfs,nx);
iterNum = zeros(1,max_imfs);
fvs = zeros(max_imfs,max_iter);

end

% Check whether there are enough (3) extrema to continue the decomposition
function stop = stopemd(xs, x_energy)
[indmin,indmax] = extr(xs);
peak = length(indmin) + length(indmax);
ratio = sum(xs.^2)/x_energy;
stop = peak < 3 | ratio < 0.001;
end

% Compute mean function of x in EMD
function [m, n_extr] = emd_mean(x, smooth_mode, ext_ratio)
% find extremum indices
[indmin, indmax, ~] = extr(x);

% total amount of extrema
n_extr = length(indmin)+length(indmax);

if n_extr < 3
    m = [];
    return
end

% extend original data to refrain end effect
[ext_indmin,ext_indmax,ext_x,cut_index] = extend(x, indmin, indmax, ext_ratio);

% compute local mean 
switch smooth_mode
    case 'spline'
        l = length(ext_x);
        ext_indmax(ext_indmax==1) = [];
        ext_indmax(ext_indmax==l) = [];
        ext_indmin(ext_indmin==1) = [];
        ext_indmin(ext_indmin==l) = [];
        max_ext_x = interp1([1,ext_indmax,l], ext_x([1,ext_indmax,l]),...
            1:l, 'spline'); % pchip
        min_ext_x = interp1([1,ext_indmin,l], ext_x([1,ext_indmin,1]),...
            1:l, 'spline');
        ext_m = 0.5*(max_ext_x+min_ext_x);
        m = ext_m(cut_index(1):cut_index(2));
    otherwise
        error(['mean computation method must be one of the following: ma,',...
            ' spline',' maspline'])
end
end

% sifting stopping criterion
function [stop_flag_sifting_process,fv_i] = is_sifting_process_stop(m, s, j, fv_i, ssc)
df = m;
[indmin, indmax, indzer] = extr(s);
lm = length(indmin);
lM = length(indmax);
nem = lm + lM;
nzm = length(indzer);
switch ssc
    
    case 'liu' % local optimal iteration.
        fv_i(j) =rms(df)+ abs(kurtosis(df)-3);
        if j >= 3 && abs(nzm-nem)<2
            if ((fv_i(j) >= fv_i(j-1)) && (fv_i(j-1) >= fv_i(j-2)))
                stop_flag_sifting_process = 1;
                return;
            end
        end
        
end
stop_flag_sifting_process = 0;

end

% Plot IMFs
function emdplot(x,imf,fs)
t = (1:size(imf,2))/fs;
imfn = size(imf,1);
figure
subplot(imfn+1,1,1);
plot(t,x);
title('Input Signal');
for imfi = 1:imfn
    subplot(imfn+1,1,imfi+1);
    plot(t,imf(imfi,:));
    if imfi < imfn
        title(['IMF',num2str(imfi)] );
    else
        title('Residual');
    end
    
end

end

% Compute the index of orthogonality
% ** Copied from emd toolbox by G.Rilling and P.Flandrin
% ** http://perso.ens-lyon.fr/patrick.flandrin/emd.html
function ort = io(x,imf)
% ort = IO(x,pfs) computes the index of orthogonality
%
% inputs : - x   : analyzed signal
%          - pfs  : production function

n = size(imf,1);
s = 0;

for i = 1:n
    for j =1:n
        if i~=j
            s = s + abs(sum(imf(i,:).*conj(imf(j,:)))/sum(x.^2));
        end
    end
end

ort = 0.5*s;
end

% Extracts the indices of extrema
% ** Copied from emd toolbox by G.Rilling and P.Flandrin
% ** http://perso.ens-lyon.fr/patrick.flandrin/emd.html
function [indmin, indmax, indzer] = extr(x)

m = length(x);

if nargout > 2
    x1 = x(1:m-1);
    x2 = x(2:m);
    indzer = find(x1.*x2<0);
    
    if any(x == 0)
        iz = find( x==0 );
        %     indz = [];
        if any(diff(iz)==1)
            zer = x == 0;
            dz = diff([0 zer 0]);
            debz = find(dz == 1);
            finz = find(dz == -1)-1;
            indz = round((debz+finz)/2);
        else
            indz = iz;
        end
        indzer = sort([indzer indz]);
    end
end

d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;

% when two or more successive points have the same value we consider only
% one extremum in the middle of the constant area (only works if the signal
% is uniformly sampled)

if any(d==0)
    
    imax = [];
    imin = [];
    
    bad = (d==0);
    dd = diff([0 bad 0]);
    debs = find(dd == 1);
    fins = find(dd == -1);
    if debs(1) == 1
        if length(debs) > 1
            debs = debs(2:end);
            fins = fins(2:end);
        else
            debs = [];
            fins = [];
        end
    end
    if ~isempty(debs)
        if fins(end) == m
            if length(debs) > 1
                debs = debs(1:(end-1));
                fins = fins(1:(end-1));
                
            else
                debs = [];
                fins = [];
            end
        end
    end
    lc = length(debs);
    if lc > 0
        for k = 1:lc
            if d(debs(k)-1) > 0
                if d(fins(k)) < 0
                    %           imax = [imax round((fins(k)+debs(k))/2)];
                end
            else
                if d(fins(k)) > 0
                    %           imin = [imin round((fins(k)+debs(k))/2)];
                end
            end
        end
    end
    
    if ~isempty(imax)
        indmax = sort([indmax imax]);
    end
    
    if ~isempty(imin)
        indmin = sort([indmin imin]);
    end
    
end
end

% Extend original data to refrain end effect
% ** Modified on emd by G.Rilling and P.Flandrin
% ** http://perso.ens-lyon.fr/patrick.flandrin/emd.html

function [ext_indmin, ext_indmax, ext_x, cut_index] = extend(x, indmin,...
    indmax, ext_ratio)
if ext_ratio == 0 % do not extend x
    ext_indmin = indmin;
    ext_indmax = indmax;
    ext_x = x;
    cut_index = [1,length(x)];
    return
end
nbsym = ceil(ext_ratio*length(indmax)); % number of extrema in extending end
xlen = length(x);
t = 1:xlen;

% boundary conditions for interpolations :

% left end extend
if indmax(1) < indmin(1) % first extremum is max
    if x(1) > x(indmin(1)) % first point > first min extremum
        lmax = fliplr(indmax(2:min(end,nbsym+1)));
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = indmax(1);
    else                  % first point < first min extremum
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
        lsym = 1;
    end
    
else                     % first extremum is maxmum
    
    if x(1) < x(indmax(1)) % first point < first maxmum
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = fliplr(indmin(2:min(end,nbsym+1)));
        lsym = indmin(1);
    else                   % first point > first minimum
        lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = 1;
    end
end

% right end
if indmax(end) < indmin(end) % last extremum is minimum
    if x(end) < x(indmax(end)) % last point < last maxmum
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = fliplr(indmin(max(end-nbsym,1):end-1));
        rsym = indmin(end);
    else                       % last point > last maxmum
        rmax = [xlen, fliplr(indmax(max(end-nbsym+2,1):end))];
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = xlen;
    end
else                         % last extremum is maxmum
    if x(end) > x(indmin(end)) % last point > last minimum
        rmax = fliplr(indmax(max(end-nbsym,1):end-1));
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = indmax(end);
    else                       % last point < last minimum
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = [xlen, fliplr(indmin(max(end-nbsym+2,1):end))];
        rsym = xlen;
    end
end

tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough
if tlmin(1) > t(1) || tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,nbsym)));
    else
        lmin = fliplr(indmin(1:min(end,nbsym)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    %     tlmin = 2*t(lsym)-t(lmin);
    %     tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(xlen) || trmax(end) < t(xlen)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
    else
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
    end
    if rsym == xlen
        error('bug')
    end
    rsym = xlen;
    %     trmin = 2*t(rsym)-t(rmin);
    %     trmax = 2*t(rsym)-t(rmax);
end

l_end = max(max(lmax, lmin));
r_end = min(min(rmax, rmin));

new_lmax = l_end+1-lmax;
new_lmin = l_end+1-lmin;
new_rmax = rsym-rmax;
new_rmin = rsym-rmin;
lx_length = l_end-lsym;
lx = fliplr(x(lsym+1:l_end));
rx = fliplr(x(r_end:rsym-1));

ext_x = [lx, x(lsym:rsym), rx];
ext_indmin = [new_lmin,indmin+lx_length-lsym+1,new_rmin+lx_length-lsym+1+...
    rsym];
ext_indmax = [new_lmax,indmax+lx_length-lsym+1,new_rmax+lx_length-lsym+1+...
    rsym];

% Index for cutting extension of x
cut_index = [lx_length-lsym+2, length(x)+lx_length-lsym+1];
end

