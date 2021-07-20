% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine that obtains the normalized computation time as the window size
% increases. Figure 4(f)
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

disp('------------------------------------');
disp('This routine could take some minutes');
disp('------------------------------------');

%% Add Paths

% Routine for obtaining the behavior of the adaptive algorithms
addpath('RecursiveMyriadFilters/');
addpath('OthersNonlinearFilters/');


%%
trials = 10;
t1 = zeros(5,trials);
variance = 0.1;

Mv = [15 30 45 60 75 90];
Nv1 = [10 20 30 40 50 60];
Nv2 = [5 10 15 20 25 30];

for ii = 1 : length(Mv)
    su = sign(randn(1,4047));
    so =  su(1:4000);
    s1 = so + laprnd(1,4000,0,sqrt(variance));
    KMAE = 1;
    u = 0.001;

    % Building the desired signal for adaptive filter algorithms
    win = Mv(ii);
    wfir = fir1(win-1,[.075 .125]);
    d = filter(wfir,1,so);
    ax = d(1:floor(Mv(ii)/2)-1);
    d=[d(floor(Mv(ii)/2):end) ax];

    % Adaptive recursive weighted median filter (RWM filter)
    disp(['--------------------------------']);
    disp(['Adaptive RWM filter Started']);
    tic;
    M = Nv1(ii); N = Nv2(ii);
    g = (1/(M+N))*ones(1,M);
    h = (1/(M+N))*ones(1,N);
    [gwm,hwm,ewm] = adaptiveWM(s1,d,g,h,0.2);
    t = toc;
    disp(['Elapsed time: ' num2str(t) ' seconds.']);
    disp(['Adaptive RWM filter Finished']);
    disp(['--------------------------------']);

    % Adaptive (non recursive) weighted myriad filter
    disp(['--------------------------------']);
    disp(['Adaptive Weighted Myriad Started']);
    tic;
    M = Mv(ii);
    w = (1/(M))*ones(1,M);
    [w1,e1] = adaptiveWMy(s1,d,w,1,u*10,1);
    t = toc;
    disp(['Elapsed time: ' num2str(t) ' seconds.']);
    disp(['Adaptive Weighted Myriad Finished']);
    disp(['--------------------------------']);

    % Adaptive scaled recursive weighted myriad filter (SRWMy filter)
    disp(['--------------------------------'])
    disp(['Adaptive SRWMyF Started']);
    tic;
    [go,ho,K1,K2] = adaptiveSRWMy(s1,d,g,h,1,1,u);
    t = toc;
    disp(['Elapsed time: ' num2str(t) ' seconds.']);
    disp(['Adaptive SRWMyF Finished']);
    disp(['--------------------------------']);

   
    % Adaptive scaled recursive hybrid myriad filter (SRHMy filter)
    disp(['--------------------------------'])
    disp(['Adaptive SRHMyF Started']);
    tic;
    [g6,h6,K16] = adaptiveSRHMy(s1,d,g,h,1,u);
    t = toc;
    disp(['Elapsed time: ' num2str(t) ' seconds.']);
    disp(['Adaptive SRHMyF Finished']);
    disp(['--------------------------------']);

    tic;
    for jj = 1:trials
        disp(['Window Size: ' num2str(Mv(ii)) '. Iteration: ' num2str(jj)]);
        fs = 2000;
        n = 0:1/fs:1.0-1/fs;
        s = chirp(n,0,1,400);
        s1 = s + sqrt(variance)*randn(1,length(s));

        % FIR filter parameters
        win1 = 96;
        wfir1 = fir1(win1-1,[.075 .125]);

        % Weighted Myriad Filter on the chirp signal
        tic;
        y1 = wmyFilter(s1,w1,1);
        t1(1,jj) = toc;
        
        % Scaled Recursive Weighted Myriad Filter on the chirp signal
        tic;
        y3 = srwmyFilter(s1,go,ho,K1,K2);
        t1(2,jj) = toc;
        
        % Scale Recursive Hybrid Myriad Filter on the chirp signal
        tic;
        y5 = srhmyFilter(s1,g6,h6,K16);
        t1(3,jj) = toc;
        
        % Dual Weighted Iterative Truncate Mean Filter (DWITM filter)
        tic;
        si      = [zeros(1,win1-1) s1];
        for i = 1:length(s1)
            samples_dwitm   = si(i:i+win1-1);
            weights_dwitm   = wfir1;
            y6(i)           = dwitm(samples_dwitm, weights_dwitm);
        end
        t1(4,jj) = toc;
        aux_y6  = y6(1:47); %Revisar fvtool...
        y6      = [y6(47+1:end) aux_y6];
        
        % Recursive Weighted Median Filter on the chirp signal
        tic;
        y7 = rwmFilter(s1,gwm,hwm);
        t1(5,jj) = toc;
        
        t1(:,jj) = t1(:,jj)*(1/t1(2,jj));
    end
    toc;
    t1_mean(ii,:) = mean(t1,2)
end

plot(Mv,t1_mean', 'LineWidth',2);
axis([12 120 0 4]);
ylabel('Normalized Computation time');
xlabel('Window Size');
legend('WMy', 'SRWMy', 'SRHMy', 'DWITM', 'RWM');
axis('tight');