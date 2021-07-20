function [g, h, ierror] = adaptiveRWM(trainingSignal, desiredSignal, g, h, initialMu)
% adaptiveRWM obtains the filter coefficients of a recursive weighted
% median filter using the fast LMA based steepest descent adaptive
% algorithm [1].
%
%   [g, h, ierror] = adaptiveRWM(trainingSignal, desiredSignal, g, h,...
%                   initialMu)
%
%   Inputs:
%   trainingSignal  = input signal vector
%   desiredSignal   = desired signal vector
%   g               = initial feedforward weight vector
%   h               = initial feedback weight vector
%   initialMu       = initial step size
%
%   Outputs:
%   g               = final feedforward weight vector
%   h               = final feedback weight vector
%   ierror          = Iterative error update
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


N  = length(trainingSignal);
M2  = length(h);
m = 1;

for n = M2 : N - M2
    % Variable step-size
    u = initialMu*exp(-(m-1)/1000);
    
    windowg = trainingSignal(n - M2 + 1:n + M2);
    windowd = desiredSignal(n - M2 + 1:n);

    x = sign(g).*windowg;
    z = sign(h).*windowd;
    ag = abs(g);
    ah = abs(h);
    xv = [x z];
    wv = [ag ah];
    
    Y = wmedian(xv,wv,mean(xv));     
    
    ierror = desiredSignal(n+1) - Y;
    
    g = g + u * ierror * sign(g) .* sign(x - Y);
    h = h + u * ierror * sign(h) .* sign(z - Y);
    
    m = m + 1;
end

return;