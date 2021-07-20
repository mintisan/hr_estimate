function outputSignal = rwmyFilter(inputSignal,g,h,K1,K2)
%   rwmyFilter filters the inputSignal vector using a recursive weighted
%   myriad (RWMy) filter which is defined by the following parameters:
%
%   g               = feedforward weight vector
%   h               = feedback weight vector
%   K1              = feedforward linear parameter
%   K2              = feedback linear parameter
%
%   The outputSignal vector contains the components of the filtered signal.
%
%   outputSignal = rwmyFilter(inputSignal,g,h,K1,K2)
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

% Initial settings
N1              = length(inputSignal);
outputSignal    = zeros(1,N1);
M2              = length(h);

% Obtaining the filtered signal
for i = M2 : N1 - M2
    windowg             = inputSignal(i - M2 + 1:i + M2);
    windowd             = outputSignal(i - M2 + 1:i);
    outputSignal(i + 1) = recursiveWeightedMyriadFPSII(windowg,windowd,g,h,K1,K2);
end