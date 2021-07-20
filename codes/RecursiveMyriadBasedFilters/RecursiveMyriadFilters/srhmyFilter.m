function outputSignal = srhmyFilter(inputSignal, g, h, K1)
% srhmyFilter filters the inputSignal vector using a recursive hybrid
% myriad (RHMy) filter which is defined by the following parameters:
%
%   g               = feedforward weight vector
%   h               = feedback weight vector
%   K1              = feedforward linear parameter
%
%   outputSignal = srhmyFilter(inputSignal, g, h, K1)
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

N1  = length(inputSignal);
M2  = length(h);
outputSignal   = zeros(N1,1);

for i = M2 : N1 - M2
    windowg = inputSignal(i - M2 + 1:i + M2);
    windowd = outputSignal(i - M2 + 1:i)';
    outputSignal(i+1)  = (sum(abs(g))+sum(abs(h)))*...
        recursiveHybridMyriadFPSII(windowg,windowd,g,h,K1);
end