% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLC filtering on a single realization of the transmitted signal and a
% single realization of the additive noise. Figure 8.
%
%   Reference: 
%
%   [1] Ramirez, J., & Paredes, J. (2016). Recursive Weighted Myriad Based
%   Filters and their Optimizations. IEEE Transactions on Signal
%   Processing, 64(15), 4027-4039.
%
%   Author:
%   Juan Marcos Ramirez, M.S.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   September, 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%% Add Paths

addpath('RecursiveMyriadFilters/');
addpath('OthersNonlinearFilters/');

%%
%--------------------------------------------
% Training stage
%--------------------------------------------

symbols = sign(randn(200,1)).*randi(2,200,1);
alpha = 0.50;
dispersion = 0.01.^alpha;
kd = sqrt(alpha/(2-alpha))*(dispersion^(1/alpha));

so = [];
for ii = 1:length(symbols)
    so = [so symbols(ii)*ones(1,50)];
end

s1 = so + astable(1,10000,alpha,0,dispersion,0);
u = 0.001;

K1f = kd;
K2f = 25*K1f;
Kf  = K1f;

% Adaptive recursive weighted myriad filter (RWMy filter)
M = 8; N = 4;
g = (1/(M+N))*ones(1,M);
h = (1/(M+N))*ones(1,N);
[g2,h2,e2] = adaptiveRWMy(s1,so,g,h,K1f,K2f,u);

% Adaptive recursive hybrid myriad filter (RHMy filter)
[g4,h4] = adaptiveRHMy(s1,so,g,h,Kf,u);


symbols = sign(randn(20,1)).*randi(2,20,1);
so = [];
for ii = 1:length(symbols)
    so = [so symbols(ii)*ones(1,50)];
end
s1 = so + astable(1,1000,alpha,0,dispersion,0);

W = N + M;
y0 = zeros(size(s1));
y1 = zeros(size(s1));
z0 = zeros(size(s1));
z1 = zeros(size(s1));
z3 = zeros(size(s1));
z4 = zeros(size(s1));

for ii = W:length(s1);
    window = s1(ii - W + 1:ii);
    y0(ii) = median(window);
    y1(ii) = weightedMyriadFPSII(window,ones(1,W), K1f);
    z0(ii) = ITM(window,1);
    z1(ii) = ITTM(window,3);
    z3(ii) = MeridianFilter(window, 1);
    z4(ii) = MGCFilter(window, 0.756, 0.896);
end

% Recursive Weighted Myriad Filter on the chirp signal 
y2 = rwmyFilter(s1,g2,h2,K1f,K2f);
data(:,8) = y2';

% Recursive Hybrid Myriad Filter on the chirp signal 
y4 = rhmyFilter(s1,g4,h4,Kf);
data(:,10) = y4';

MAE(1,1) = mean(abs(so(102:899) - z0(105:902)));
MAE(1,2) = mean(abs(so(102:899) - z1(105:902)));
MAE(1,3) = mean(abs(so(102:899) - y0(107:904)));
MAE(1,4) = mean(abs(so(102:899) - y1(107:904)));
MAE(1,5) = mean(abs(so(102:899) - z3(107:904)));
MAE(1,6) = mean(abs(so(102:899) - z4(107:904)));
MAE(1,7) = mean(abs(so(102:899) - y2(103:900)));
MAE(1,8) = mean(abs(so(102:899) - y4(103:900)));

NMSE(1,1) = 10 * log10(mean((so(102:899) - z0(105:902)).^2)/mean(so(101:900).^2));
NMSE(1,2) = 10 * log10(mean((so(102:899) - z1(105:902)).^2)/mean(so(101:900).^2));
NMSE(1,3) = 10 * log10(mean((so(102:899) - y0(107:904)).^2)/mean(so(101:900).^2));
NMSE(1,4) = 10 * log10(mean((so(102:899) - y1(107:904)).^2)/mean(so(101:900).^2));
NMSE(1,5) = 10 * log10(mean((so(102:899) - z3(107:904)).^2)/mean(so(101:900).^2));
NMSE(1,6) = 10 * log10(mean((so(102:899) - z4(107:904)).^2)/mean(so(101:900).^2));
NMSE(1,7) = 10 * log10(mean((so(102:899) - y2(103:900)).^2)/mean(so(101:900).^2));
NMSE(1,8) = 10 * log10(mean((so(102:899) - y4(103:900)).^2)/mean(so(101:900).^2));


figure;
subplot(5,2,1);
plot(so);
title('Transmitted signal')
axis([1 1000 -3 3])

subplot(5,2,2);
plot(s1);
title('Noisy signal')
axis([1 1000 -3 3])

subplot(5,2,3);
plot(z0(105:902));
title(['ITM Filter. MAE = ' num2str(MAE(1,1))]);
axis([1 800 -3 3])

subplot(5,2,4);
plot(z1(105:902));
title(['ITTM Filter. MAE = ' num2str(MAE(1,2))]);
axis([1 800 -3 3])

subplot(5,2,5);
plot(y0(107:904));
title(['Median Filter. MAE = ' num2str(MAE(1,3))]);
axis([1 800 -3 3])

subplot(5,2,6);
plot(y1(107:904));
title(['Myriad Filter. MAE = ' num2str(MAE(1,4))]);
axis([1 800 -3 3])

subplot(5,2,7);
plot(z3(107:904));
title(['Meridian Filter. MAE = ' num2str(MAE(1,5))]);
axis([1 800 -3 3])

subplot(5,2,8);
plot(z4(107:904));
title(['M-GC Filter. MAE = ' num2str(MAE(1,6))]);
axis([1 800 -3 3])

subplot(5,2,9);
plot(y2(103:900));
title(['RWMy Filter. MAE = ' num2str(MAE(1,7))]);
axis([1 800 -3 3])

subplot(5,2,10);
plot(y4(103:900));
title(['RHMy Filter. MAE = ' num2str(MAE(1,8))]);
axis([1 800 -3 3])