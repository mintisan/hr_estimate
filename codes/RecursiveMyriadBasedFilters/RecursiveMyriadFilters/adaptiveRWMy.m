function [gFinal, hFinal, ierror,gTrack,hTrack] = adaptiveRWMy(trainingSignal,desiredSignal, g, h, K1, K2, initialMu)
% adaptiveRWMy obtains the filter coefficients of a recursive weighted
% myriad (RWMy) filter using a LMA based steepest descent adaptive
% algorithm [1]. Furthermore, the expression for computing each
% instantaneous gradient is obtained under the equation error formulation
% framework.
%
%   [gFinal, hFinal, ierror] = adaptiveRWMy(trainingSignal,...
%                               desiredSignal, g, h, K1, K2, initialMu)
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
%   gFinal          = final feedforward weight vector
%   hFinal          = final feedback weight vector
%   ierror          = Iterative error update
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
%   August, 2016
%
%   Copyright 2016 Juan Marcos Ramirez Rondon.  [juanmarcos26-at-gmail.com]


% Initial Settings
N       = length(trainingSignal);
M1      = length(g);
M2      = length(h);
ierror  = zeros(N - M1 + 1, 1);
gTrack  = zeros(N - M1 + 1,M1);
hTrack  = zeros(N - M1 + 1,M2);
m = 1;

% Iterative update of the filter coefficients
for n = M2 : N - M2
    % Variable step-size
    u = initialMu*exp(-(m-1)/1000);
    
    % Computing the output error
    windowg = trainingSignal(n - M2 + 1:n + M2);
    windowd = desiredSignal(n - M2 + 1:n);
    Y = recursiveWeightedMyriadFPSII(windowg,windowd,g,h,K1,K2);
    ierror(m) = Y - desiredSignal(n+1);
    
    % Instantaneous gradient
    x = sign(g).*windowg;
    z = sign(h).*windowd;
    ag = abs(g);
    ah = abs(h);
    
    
    num1 = K1^2*(sign(g).*(x - Y))./((K1^2 + ag.*((x - Y).^2)).^2);
    num2 = K2^2*(sign(h).*(z - Y))./((K2^2 + ah.*((z - Y).^2)).^2);
    den = sum(ag.*((K1^2 - ag.*((x - Y).^2))./(K1^2 + ag.*((x - Y).^2)).^2)) +...
        sum(ah.*((K2^2 - ah.*((z - Y).^2))./(K2^2 + ah.*((z - Y).^2)).^2));
    
    % Iterative update
    g = g - u*sign(ierror(m)).*(num1/den);
    h = h - u*sign(ierror(m)).*(num2/den);
    
    gTrack(m,:) = g;
    hTrack(m,:) = h;
   
    m = m + 1;
end

gFinal = g;
hFinal = h;
return;