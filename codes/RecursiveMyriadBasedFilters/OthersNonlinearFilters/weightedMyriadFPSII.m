function [weightedMyriad] = weightedMyriadFPSII(samples, weights, K)
% weightedMyriadFPSII computes the weighted myriad of the sample vector
% (samples) weighted by the weight vector (weights), with a linearity
% parameter given by K. The weight vector and the sample vector must have
% the same length. This function uses a fast version of the fixed-point
% search algorithm for computing the weighted myriad output [1].
%
%   [weightedMyriad] = weightedMyriadFPSII(samples,weights,K)
%
%   Reference: 
%
%   [1] Kalluri, S., & Arce, G. R. (2000). Fast algorithms for weighted
%   myriad computation by fixed-point search. IEEE Transactions on Signal
%   Processing, 48(1), 159-171.
%
%   Authors:
%   Antonny Alarcon, B.S.E.E., B.S.M.E.
%   Juan Marcos Ramirez, M.S.B.E.
%   Universidad de Los Andes, Merida, Venezuela
%   email: juanra@ula.ve, juanmarcos26@gmail.com
%
%   Date:
%   August, 2016

signedSamples       = sign(weights).*samples;
magnitudeWeights    = abs(weights);
N                   = length(signedSamples); 
k2                  = K^2./magnitudeWeights;

tempCostFunction = zeros(N,N);
for n = 1:N
    tempCostFunction(n,:) = K^2 + ...
        magnitudeWeights.*(signedSamples - samples(n)).^2;
end

costFunction        = sum(log(tempCostFunction));
[~,I]               = min(costFunction);

den                 = k2 + (signedSamples-samples(I)).^2; 
previousEstimate    = samples(I);
currentEstimate     = sum(signedSamples./den)/sum(1./den);

while abs(currentEstimate - previousEstimate) > 0.000001
    den = k2 + (signedSamples-currentEstimate).^2;
    previousEstimate = currentEstimate;
    currentEstimate = sum(signedSamples./den)/sum(1./den);
end

weightedMyriad = currentEstimate;
return;