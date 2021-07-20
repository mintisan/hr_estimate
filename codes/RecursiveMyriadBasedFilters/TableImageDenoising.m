% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Routine
% Routine for obtaining the performance of various nonlinear filters, under
% different metrics, in the context of image denoising. Table II.
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

addpath('RecursiveMyriadImages/');
%% Loading Images

Ii = imread('Data/girl512.bmp');
n_row = size(Ii, 1);
n_col = size(Ii, 2);

% Salt and pepper noise
In = imnoise(Ii, 'salt & pepper', 0.1);

Iid = double(Ii);
Ind = double(In);

%% Adaptive filter algorithms

patch_size = 64;
Ii_block = Iid(end - patch_size +1:end, end - patch_size + 1:end);
In_block = Ind(end - patch_size +1:end, end - patch_size + 1:end);
window_size = 5;
Weights = ones(window_size^2,1)/(window_size^2);

% Adaptive Weighted Median Filter
[Weights, ewm] = AdaptiveWeightedMedianImages(In_block, Ii_block, Weights, window_size, 0.0001);

gi = ones(1, window_size^2)/(window_size^2 + floor(window_size^2/2));
hi = ones(1, floor(window_size^2/2))/(window_size^2 + floor(window_size^2/2));
K1_rwmy = 5;
K2_rwmy = 15;
u0 = 0.0001;

% Adaptive Recursive Weighted Myriad Filter
[g0, h0]            = AdaptiveRWMyImages(Ii_block, In_block, gi, hi, K1_rwmy, K2_rwmy, u0, window_size);

K1i = 0.05;
u0 = 0.000001;
% Adaptive Recursive Hybrid Myriad Filter
[g2, h2, K12, eh2] = AdaptiveRHMyImages(Ii_block, In_block, gi, hi, K1i, u0, window_size);

% Adaptive Recursive Weighted Median Filter
[grwm, hrwm, erwm] = AdaptiveRecursiveWeightedMedianImages(Ii_block, In_block, gi, hi, 0.001, window_size);


trials = 10;
%% Filtering noisy image
disp(['Number of realizations: ' num2str(trials)])
for kk = 1:trials
    disp(['Iteration:   ' num2str(kk)]);
    
    In = imnoise(Ii, 'salt & pepper', 0.2);
    Iid = double(Ii);
    Ind = double(In);
    
    Ind1 = [Ind(:,1:floor(window_size/2)) Ind Ind(:,end - floor(window_size/2) +1:end)];
    Ind1 = [Ind1(1:floor(window_size/2),:); Ind1; Ind1(end - floor(window_size/2) +1:end,:)];
    
    
    tic;
    Is0  = WITMFilterImage(Ind1, Weights, window_size);
    t(1,kk) = toc;
    disp(['WITM Filter. Elapsed time:   ' num2str(t(1,kk),2) ' seconds']);
    
    tic;
    Is1  = RWMyImages(Ind1, g0, h0, K1_rwmy, K2_rwmy, window_size);
    t(2,kk) = toc;
    disp(['RWMy Filter. Elapsed time:   ' num2str(t(2,kk),2) ' seconds']);
    
    tic;
    Is2  = RHMyImages(Ind1, g2, h2, K12, window_size);
    t(3,kk) = toc;
    disp(['RHMy Filter. Elapsed time:   ' num2str(t(3,kk),2) ' seconds']);

    tic;
    Is3 = RecursiveWeightedMedianImages(Ind1, grwm, hrwm, window_size);
    t(4,kk) = toc;
    disp(['RWM Filter. Elapsed time:   ' num2str(t(4,kk),2) ' seconds']);

    tic;
    Is4  = ITMFilterImage(Ind1, window_size);
    t(5,kk) = toc;
    disp(['ITM Filter. Elapsed time:   ' num2str(t(5,kk),2) ' seconds']);
    
    tic;
    Is5  = ITTMFilterImage(Ind1, window_size);
    t(6,kk) = toc;
    disp(['ITTM Filter. Elapsed time:   ' num2str(t(6,kk),2) ' seconds']);
 
    tic;
    Is6  = WeightedMedianFilterImage(Ind1, Weights, window_size);
    t(7,kk) = toc;
    disp(['WM Filter. Elapsed time:   ' num2str(t(7,kk),2) ' seconds']);

    tic;
    Is7  = WeightedMeanFilterImage(Ind1, Weights, window_size);
    t(8,kk) = toc;
    disp(['WMean Filter. Elapsed time:   ' num2str(t(8,kk),2) ' seconds']);
    disp('---------------------------------------');
    
    Ir0 = Is0(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir1 = Is1(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir2 = Is2(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir3 = Is3(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir4 = Is4(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir5 = Is5(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir6 = Is6(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    Ir7 = Is7(ceil(window_size/2):end - floor(window_size/2),ceil(window_size/2):end - floor(window_size/2));
    
    % Computing errors
    MAE(1,kk) = mean(abs(Iid(:) - Ind(:)));
    MAE(2,kk) = mean(abs(Iid(:) - Ir0(:)));
    MAE(3,kk) = mean(abs(Iid(:) - Ir1(:)));
    MAE(4,kk) = mean(abs(Iid(:) - Ir2(:)));
    MAE(5,kk) = mean(abs(Iid(:) - Ir3(:)));
    MAE(6,kk) = mean(abs(Iid(:) - Ir4(:)));
    MAE(7,kk) = mean(abs(Iid(:) - Ir5(:)));
    MAE(8,kk) = mean(abs(Iid(:) - Ir6(:)));
    MAE(9,kk) = mean(abs(Iid(:) - Ir7(:)));
    
    MSE(1,kk) = mean((Iid(:) - Ind(:)).^2);
    MSE(2,kk) = mean((Iid(:) - Ir0(:)).^2);
    MSE(3,kk) = mean((Iid(:) - Ir1(:)).^2);
    MSE(4,kk) = mean((Iid(:) - Ir2(:)).^2);
    MSE(5,kk) = mean((Iid(:) - Ir3(:)).^2);
    MSE(6,kk) = mean((Iid(:) - Ir4(:)).^2);
    MSE(7,kk) = mean((Iid(:) - Ir5(:)).^2);
    MSE(8,kk) = mean((Iid(:) - Ir6(:)).^2);
    MSE(9,kk) = mean((Iid(:) - Ir7(:)).^2);

    PSNR(1,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ind(:)).^2)));
    PSNR(2,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir0(:)).^2)));
    PSNR(3,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir1(:)).^2)));
    PSNR(4,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir2(:)).^2)));
    PSNR(5,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir3(:)).^2)));
    PSNR(6,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir4(:)).^2)));
    PSNR(7,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir5(:)).^2)));
    PSNR(8,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir6(:)).^2)));
    PSNR(9,kk) = 10 * log10((255^2)/(mean((Iid(:) - Ir7(:)).^2)));
end

t_mean = mean(t,2);
MAE_mean = mean(MAE,2);
MSE_mean = mean(MSE,2);
PSNR_mean = mean(PSNR,2);

disp('Performance metrics');
disp('--------------------------------------------------');
disp(['Filter       MAE         MSE         PSNR'])
disp('--------------------------------------------------');
disp(['WMean        ' num2str(MAE_mean(9),4) '       ' num2str(MSE_mean(9),4) '       ' num2str(PSNR_mean(9),4)] );
disp(['WM           ' num2str(MAE_mean(8),4) '       ' num2str(MSE_mean(8),4) '       ' num2str(PSNR_mean(8),4)]);
disp(['ITM          ' num2str(MAE_mean(6),4) '       ' num2str(MSE_mean(6),4) '       ' num2str(PSNR_mean(6),4)] );
disp(['ITTM         ' num2str(MAE_mean(7),4) '       ' num2str(MSE_mean(7),4) '       ' num2str(PSNR_mean(7),4)]);
disp(['WITM         ' num2str(MAE_mean(2),4) '       ' num2str(MSE_mean(2),4) '       ' num2str(PSNR_mean(2),4)]);
disp(['RWM          ' num2str(MAE_mean(5),4) '       ' num2str(MSE_mean(5),4) '       ' num2str(PSNR_mean(5),4)]);
disp('--------------------------------------------------');
disp(['RWMy         ' num2str(MAE_mean(3),4) '       ' num2str(MSE_mean(3),4) '       ' num2str(PSNR_mean(3),4)]);
disp(['RHMy         ' num2str(MAE_mean(4),4) '       ' num2str(MSE_mean(4),4) '       ' num2str(PSNR_mean(4),4)]);
disp('--------------------------------------------------');