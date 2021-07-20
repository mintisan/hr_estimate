function outputSignal = wmyFilter(inputSignal, weights, K)
%   wmyFilter filters the inputSignal vector using a weighted myriad
%   filtering structure defined by both the vector of weights and the
%   linearity parameter (K).
%
%   outputSignal = wmyFilter(inputSignal, weights, K)
%
%   Reference:
%   [1] Kalluri, S., Arce, G. Robust frequency-selective filtering
%   using weighted myriad filters admitting real-valued weights. Signal
%   Processing, IEEE Transactions on, 49(11), 2721-2733. 2001.
%
%   Authors:
%   Antonny Alarcon, B.S.E.E., B.S.M.E.
%   Juan Marcos Ramirez, M.S.B.E.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   August, 2016

% Initial settings
N1              = length(inputSignal);  % Input signal length
M               = length(weights);      % Filter length
outputSignal    = zeros(N1, 1);

% Filtering the input signal
for i = 1 : N1 - M + 1
    window = inputSignal(i:i + M -1);    
    outputSignal(i) = weightedMyriadFPSII(window, weights, K);                   %Causal structure
end