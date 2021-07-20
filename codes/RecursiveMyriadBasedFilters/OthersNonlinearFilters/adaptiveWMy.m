function [weights, K, ierror] = adaptiveWMy(trainingSignal, desiredSignal, weights, K, mu, algorithm)
% adaptiveWMy obtains the filter coefficients of a weighted myriad
% filter using a steepest descent based adaptive algorithm.
%
% [weights, K, ierror] = adaptiveWMy(trainingSignal, desiredSignal,...
%                           weights, K, mu, algorithm)
%
%   Inputs:
%   trainingSignal  = input signal
%   desiredSignal   = desired signal
%   weights         = initial weight vector
%   K               = linearity parameter
%   mu              = initial step size parameter
%   algorithm       = this option select the adaptive algorithm
%                 'LMS' -> Least mean square based algorithm [1]
%                 'LMA' -> Least mean absolute based algorithm [2]
%
%   Outputs:
%   weights         = output weight vector
%   K               = linearity parameter
%   ierror          = Iterative error update
%
%   References:
%   [1] Kalluri, S., Arce, G. Robust frequency-selective filtering
%   using weighted myriad filters admitting real-valued weights. Signal
%   Processing, IEEE Transactions on, 49(11), 2721-2733. 2001.
%
%   [2] Delgado, M., Reyes, H., Bautista, L., Ramirez, Juan. Filtrado
%   Robusto de Imagenes con Estructuras Optimizadas basadas en Myriad.
%   Tercera Conferencia Nacional de Computacion, Informatica y Sistemas,
%   CoNCISa 2015, Valencia, Venezuela. 2015
%
%   Authors:
%   Juan Marcos Ramirez, M.S.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   August, 2016

% Initial settings
N   = length(trainingSignal);
M   = length(weights);
m   = 0;
ierror   = zeros(N - M + 1, 1);

% Loop of the adative algorithm
for i = 1: N - M + 1
    % Computing the current step size parameter
    muVariable = mu *exp(-m/1000);
    m = m + 1;
    
    % Computing the output error
    window          = trainingSignal(i:i + M - 1);
    desiredOutput   = desiredSignal(i);
    Y               = weightedMyriadFPSII(window, weights,K); % Filter Output
    ierror(m)       = Y - desiredOutput;
    
    % Instataneous gradient
    x   = sign(weights).*window;
    aW  = abs(weights);   
    num = K^2*(sign(weights).*(Y-x))./((K^2 + aW.*((x - Y).^2)).^2);
    den = sum(aW.*((K^2 - aW.*((Y-x).^2))./(K^2 + aW.*((x - Y).^2)).^2));
    
    % Iterative update
    if (algorithm == 'LMS')
        weights   = weights     +  muVariable * (ierror(m)) * num/den;
    end
    
    if (algorithm == 'LMA')
        weights   = weights     +  muVariable * sign(ierror(m)) * num/den;
    end
end