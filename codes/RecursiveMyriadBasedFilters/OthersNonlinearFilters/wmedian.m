function [median, index] = wmedian(samples,weights,varargin)
%WMEDIAN Traditional Weighted Median Operator with real weights function
%
%   X = WMEDIAN(SAMPLES) returns the weighted median of SAMPLES with 
%   unitary weights, or the traditional median without interpolation.
%
%   X = WMEDIAN(SAMPLES,WEIGHTS) returns the weighted median of SAMPLES
%   weighted by WEIGHTS.
%
%   [X,I] = WMEDIAN(SAMPLES,WEIGHTS) returns the weighted median of SAMPLES
%   weighted by WEIGHTS in X and the median index I.
%
%   If SAMPLES and WEIGHTS are matrices use IND2SUB function to recover a
%   subscript from the INDEX value.
%
%   Reference:
%       A General Weighted Median Filter Structure Admitting Negative 
%       Weights
%       http://dx.doi.org/10.1109/78.735296

if nargin == 1
    weights = ones(numel(samples), 1);
end

if length(samples) ~= numel(samples)
    weights = weights(:);
    samples = samples(:);
end

samples = samples.*(sign(weights) + (weights == 0));
weights = abs(weights);

[~, sorted_index] = sort(samples,'descend');
threshold = sum(weights)/2;

partialSum = 0;
i = 1;

while partialSum < threshold
    partialSum = partialSum + weights(sorted_index(i));
    i = i + 1;
end

index = sorted_index(i - 1);
median = samples(index);