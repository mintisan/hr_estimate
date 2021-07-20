% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NMSE[dB] and MAE yielded by various nonlinear filters in the context of
% powerline communications. Figure 9.
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

disp('---------------------------------------');
disp('This routine could take several minutes');
disp('---------------------------------------');

%% Add Paths

addpath('RecursiveMyriadFilters/');
addpath('OthersNonlinearFilters/');

%% Training Stage

alpha = 0.6:0.1:1.9;
dispersion = 0.01.^alpha;
trials = 10;
kd = sqrt(alpha./(2-alpha)).*(dispersion.^(1./alpha));


for jj = 1:length(alpha)
    disp(['Characteristic Exponent (alpha): ' num2str(alpha(jj))])
    tic;
    symbols = sign(randn(200,1)).*randi(2,200,1);
    so = [];
    for ii = 1:length(symbols)
        so = [so symbols(ii)*ones(1,50)];
    end
    
    s1 = so + astable(1,10000,alpha(jj),0,dispersion(jj),0);
    u = 0.001;
    
    if alpha(jj) ~= 2
        kd =  sqrt(alpha(jj)/(2-alpha(jj)))*(dispersion(jj)^(1/alpha(jj)));
    else
        kd =  sqrt(1.93/(2-1.93))*(dispersion(jj)^(1/1.925));
    end
    
    K1f = kd;
    K2f = kd*25;
    Kf  = kd;
    
    % Adaptive recursive weighted myriad filter (RWMy filter)
    M = 8; N = 4;
    g = (1/(M+N))*ones(1,M);
    h = (1/(M+N))*ones(1,N);
    [g2,h2,e2] = adaptiveRWMy(s1,so,g,h,K1f,K2f,u);
    
    
    % Adaptive recursive hybrid myriad filter (RHMy filter)
    [g4,h4] = adaptiveRHMy(s1,so,g,h,Kf,u);
    
    
    %% Test Stage
    disp(['Number of realizations: ' num2str(trials)])
    for kk = 1:trials
        disp(['Iteration:   ' num2str(kk)]);
        symbols = sign(randn(20,1)).*randi(2,20,1);
        so = [];
        for ii = 1:length(symbols)
            so = [so symbols(ii)*ones(1,50)];
        end
        s1 = so + astable(1,1000,alpha(jj),0,dispersion(jj),0);
        
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
            y1(ii) = weightedMyriadFPSII(window,ones(1,W),K1f);
            z0(ii) = ITM(window,1);
            z1(ii) = ITTM(window,3);
            z3(ii) = MeridianFilter(window, 1);
            z4(ii) = MGCFilter(window, 0.756, 0.896);
        end
        
        % Recursive Weighted Myriad Filter on the chirp signal
        y2 = rwmyFilter(s1,g2,h2,K1f,K2f);
        
        % Recursive Hybrid Myriad Filter on the chirp signal
        y4 = rhmyFilter(s1,g4,h4,Kf);
        
        MAE(kk,1) = mean(abs(so(102:899) - z0(105:902)));
        MAE(kk,2) = mean(abs(so(102:899) - z1(105:902)));
        MAE(kk,3) = mean(abs(so(102:899) - y0(107:904)));
        MAE(kk,4) = mean(abs(so(102:899) - y1(107:904)));
        MAE(kk,5) = mean(abs(so(102:899) - z3(107:904)));
        MAE(kk,6) = mean(abs(so(102:899) - z4(107:904)));
        MAE(kk,7) = mean(abs(so(102:899) - y2(103:900)));
        MAE(kk,8) = mean(abs(so(102:899) - y4(103:900)));
        
        NMSE(kk,1) = 10 * log10(mean((so(102:899) - z0(105:902)).^2)/mean(so(101:900).^2));
        NMSE(kk,2) = 10 * log10(mean((so(102:899) - z1(105:902)).^2)/mean(so(101:900).^2));
        NMSE(kk,3) = 10 * log10(mean((so(102:899) - y0(107:904)).^2)/mean(so(101:900).^2));
        NMSE(kk,4) = 10 * log10(mean((so(102:899) - y1(107:904)).^2)/mean(so(101:900).^2));
        NMSE(kk,5) = 10 * log10(mean((so(102:899) - z3(107:904)).^2)/mean(so(101:900).^2));
        NMSE(kk,6) = 10 * log10(mean((so(102:899) - z4(107:904)).^2)/mean(so(101:900).^2));
        NMSE(kk,7) = 10 * log10(mean((so(102:899) - y2(103:900)).^2)/mean(so(101:900).^2));
        NMSE(kk,8) = 10 * log10(mean((so(102:899) - y4(103:900)).^2)/mean(so(101:900).^2));
        
    end
    mae_mean(:,jj) = mean(MAE)';
    nmse_mean(:,jj) = mean(NMSE)';
    toc;
    disp('---------------------------------------');
end

%% Diplay Results

subplot(121);
plot(alpha,nmse_mean(3:8,:)','LineWidth',2);
ylabel('NMSE[dB]'), xlabel('Characteristic Exponent');
legend('Median','Myriad','Meridian','M-GC','RWMy','RHMy');

subplot(122);
plot(alpha,mae_mean(3:8,:)','LineWidth',2);
ylabel('MAE'), xlabel('Characteristic Exponent');
legend('Median','Myriad','Meridian','M-GC','RWMy','RHMy');