function [beta] = recursiveHybridMyriadFPSII(X,Y,g,h,K)
% recursiveHybridMyriadFPSII computes the recursive hybrid myriad
% (beta) on both a subset of input samples X and a subset of previous
% outputs [1]. For obtaining the operator output the function uses a fixed
% point search approach. Input parameters are the following:
%
%   g       = feedforward weight vector (must be have the same length of X) 
%   h       = feedback weight vector (must be have the same length of Y)
%   K      = feedforward linearity parameter
%
%   [beta] = recursiveHybridMyriadFPSII(X,Y,g,h,K)
%
%   Please refer to [1] if you use this software. 
%
%   References:
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
%   Licensed under the Simplified BSD License [see BSD-2.txt]


Beta = [X mean(h.*Y)]; 

X = sign(g).*X;
Y = sign(h).*Y;
g = abs(g);
h = abs(h);

N = length(X);
T = N+1;
k1 = K^2./g; 

Costo = zeros(1,T); 

for n = 1:T  
    Costo(n) = sum(log(K^2 + g.*(X - Beta(n)).^2)) + sum(h.*(Y - Beta(n)).^2);
end

[~,I] = min(Costo); 

den = k1 + (X-Beta(I)).^2; 
band = Beta(I); 
bd = (sum(X./den) + sum(h.*Y))/(sum(1./den) + sum(h));

while abs(bd - band) > 0.000001 
    den = k1 + (X-bd).^2; 
    band = bd;
    bd = (sum(X./den) + sum(h.*Y))/(sum(1./den) + sum(h));
end

beta = bd; 
return;