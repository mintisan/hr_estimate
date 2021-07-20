function outputSignal = rwmFilter(inputSignal,g,h)
%   rwmFilter filters the inputSignal vector using a recursive weighted
%   median filter defined by both the vector of feedforward weights (g) and
%   the vector of feedback weights (h). The vector outputSignal contains
%   the samples of the filtered signal.
%
%   outputSignal = rwmFilter(inputSignal,g,h)
%
%   Reference:
%
%   [1] Arce, G. R., & Paredes, J. L. (2000). Recursive weighted median
%   filters admitting negative weights and their optimization. IEEE
%   Transactions on Signal Processing, 48(3), 768-779.
%
%   Authors:
%   Juan Marcos Ramirez, M.S.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   August, 2016


N1  = length(inputSignal);
outputSignal   = zeros(1,N1);
M2  = length(h);

for i = M2 : N1 - M2
    windowg = inputSignal(i - M2 + 1:i + M2);
    windowd = outputSignal(i - M2 + 1:i);
    
    x = sign(g).*windowg;
    z = sign(h).*windowd;
    ag = abs(g);
    ah = abs(h);
    xv = [x z];
    wv = [ag ah];
    
    outputSignal(i+1) = wmedian(xv,wv,mean(xv));
end

return;