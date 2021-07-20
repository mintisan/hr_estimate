function [g,h,ierror,gTrack,hTrack] = adaptiveRHMy(trainingSignal,desiredSignal,g,h,K1,initialMu)
% adaptiveRHMy obtains the filter coefficients of a recursive hybrid
% myriad (RHMy) filter using a LMA based steepest descent adaptive
% algorithm [1]. Furthermore, the expression for each instantaneous
% gradient is obtained under the equation error formulation framework.
%
% [g, h, ierror] = adaptiveRHMy(trainingSignal, desiredSignal, g, h,...
%                   K1, initialMu)
%
%   Inputs:
%   trainingSignal  = input signal vector
%   desiredSignal   = desired signal vector
%   g               = initial feedforward weight vector
%   h               = initial feedback weight vector
%   K1              = feedforward linear parameter
%   initialMu       = initial step size
%
%   Outputs:
%   g               = final feedforward weight vector
%   h               = final feedback weight vector
%   ierror          = Iterative error update
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


N       = length(trainingSignal);
M1      = length(g);
M2      = length(h);
ierror  = zeros(N - M1 + 1, 1);
gTrack  = zeros(N - M1 + 1,M1);
hTrack  = zeros(N - M1 + 1,M2);
m = 1;

for n = M2 : N - M2
    u = initialMu*exp(-(m-1)/1000);
    
    windowg = trainingSignal(n - M2 + 1:n + M2);
    windowd = desiredSignal(n - M2 + 1:n);
    
    x = sign(g).*windowg;
    ag = abs(g);
    ah = abs(h);
    
    Y = recursiveHybridMyriadFPSII(windowg,windowd,g,h,K1);
    
    ierror(m) = Y - desiredSignal(n+1);
    
    num1 = K1^2*(sign(g).*(x - Y))./((K1^2 + ag.*((x - Y).^2)).^2);
    num2 = windowd - sign(h)*Y;
    den = sum(ag.*((K1^2 - ag.*((x - Y).^2))./(K1^2 + ag.*((x - Y).^2)).^2)) +...
        sum(ah);
    
    g = g - u*sign(ierror(m)).*(num1/den);
    h = h - u*sign(ierror(m)).*(num2/den);
    
    gTrack(m,:) = g;
    hTrack(m,:) = h;

    m = m + 1;
end

return;