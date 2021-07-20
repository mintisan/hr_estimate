function [g,h,K1,K2,e,gTrack,hTrack,K1Track,K2Track] = adaptiveSRWMy(trainingSignal,desiredSignal,g,h,K1,K2,initialMu)
% adaptiveSRWMy obtains the filter parameters (filter coefficients and 
% linearity parameters) of a scaled recursive weighted myriad (SRWMy) filter
% using a LMA based steepest descent adaptive algorithm [1]. Furthermore,
% the expression for computing each instantaneous gradient is obtained
% under the equation error formulation
%
% [g,h,K1,K2,e] = adaptiveSRWMy(trainingSignal,desiredSignal, g, h, K1, K2,
%                               initialMu)
%
%   Inputs:
%   trainingSignal  = input signal vector
%   desiredSignal   = desired signal vector
%   g               = initial feedforward weight vector
%   h               = initial feedback weight vector
%   K1              = feedforward linear parameter
%   K2              = feedback linear parameter
%   initialMu       = initial step size
%
%   Outputs:
%   g               = final feedforward weight vector
%   h               = final feedback weight vector
%   K1              = final feedforward linearity parameter (K1)
%   K2              = final feedback linearity parameter (K2)
%   e               = Iterative error update
%
%   Please refer to [1] if you use this software.
%
%   Reference:
%
%   [1] Ramirez, J., & Paredes, J. (2016). Recursive Weighted Myriad Based
%   Filters and their Optimizations. IEEE Transactions on Signal
%   Processing, 64(15), 4027-4039.
%
%   Authors:
%   Juan Marcos Ramirez, M.S.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   September, 2016
%
%   Copyright 2016 Juan Marcos Ramirez Rondon.  [juanmarcos26-at-gmail.com]


% Initial Settings
N       = length(trainingSignal);
M1      = length(g);
M2      = length(h);
kk1     = K1^2;
kk2     = K2^2;
e       = zeros(N - M1 + 1, 1);
gTrack = zeros(N - M1 + 1,M1);
hTrack = zeros(N - M1 + 1,M2);
K1Track = zeros(1,N - M1 + 1);
K2Track = zeros(1,N - M1 + 1);

m       = 0;

% Iterative update of the filter coefficients
for i = M2: N - M2
    % Variable step-size
    us = initialMu *exp(-m/1000);
    m = m + 1;

    % Computing the output error
    windowg  = trainingSignal(i - M2 + 1:i + M2);
    windowd  = desiredSignal(i - M2 + 1: i);
    d        = desiredSignal(i + 1);
    Y   = recursiveWeightedMyriadFPSII(windowg,windowd,g,h,K1,K2);
    e(m) = (sum(abs(g))+sum(abs(h)))*Y - d;

    % Instantaneous gradient
    x = sign(g).*windowg;
    z = sign(h).*windowd;
    ag = abs(g);
    ah = abs(h);
    
    num1 = K1^2*(sign(g).*(x - Y))./((K1^2 + ag.*((x - Y).^2)).^2);
    num2 = K2^2*(sign(h).*(z - Y))./((K2^2 + ah.*((z - Y).^2)).^2);
    num3 = sum((ag.*(Y - x))./((K1^2 + ag.*((x - Y).^2)).^2));
    num4 = sum((ah.*(Y - z))./((K2^2 + ah.*((z - Y).^2)).^2));

    den = sum(ag.*((K1^2 - ag.*((x - Y).^2))./(K1^2 + ag.*((x - Y).^2)).^2)) +...
        sum(ah.*((K2^2 - ah.*((z - Y).^2))./(K2^2 + ah.*((z - Y).^2)).^2));
    
    % Steepest descent update
    g = g - us*sign(e(m)).*(sign(g).*Y + (sum(abs(g))+sum(abs(h)))*num1/den);
    h = h - us*sign(e(m)).*(sign(h).*Y + (sum(abs(g))+sum(abs(h)))*num2/den);
    kk1 = abs(kk1 - us*sign(e(m)).*((sum(abs(g))+sum(abs(h)))*num3/den));
    kk2 = abs(kk2 - us*sign(e(m)).*((sum(abs(g))+sum(abs(h)))*num4/den));

    K1 = sqrt(kk1);
    K2 = sqrt(kk2);
    
    gTrack(m,:) = g;
    hTrack(m,:) = h;
    K1Track(m)   = K1;
    K2Track(m)   = K2;
end

return