function [beta] = recursiveWeightedMyriadFPSII(X,Y,g,h,K1,K2)
% recursiveWeightedMyriadFPSII computes the recursive weighted myriad
% (beta) on both a subset of input samples X and a subset of previous
% outputs [1]. For obtaining the operator output the function uses a fixed
% point search approach. Input parameters are the following:
%
%   g       = feedforward weight vector (must be have the same length of X) 
%   h       = feedback weight vector (must be have the same length of Y)
%   K1      = feedforward linearity parameter
%   K2      = feedback linearity parameter
%
%   [beta] = recursiveWeightedMyriadFPSII(X,Y,g,h,K1,K2)
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

Beta = [X Y]; %Vector de puntos de minimizacion de la funcion costo recursivo
              %Beta(1:N) = X; %Elementos no recursivos
              %Beta((N+1):T) = Y; %Elementos recursivos
              
X = sign(g).*X;
Y = sign(h).*Y;
g = abs(g);
h = abs(h);

N = length(X);%Ancho de la muestra o del vector de observacion
M = length(Y); T = N+M;
k12 = K1^2./g; %Simplificacion matem�tica
k22 = K2^2./h;

Costo = zeros(1,T); %Inicializacion de la funcion costo

for n = 1:T  %Evaluacion de la funcion costo en las muestras de observacion y los elementos recursivos
    Costo(n) = prod(K1^2 + g.*(X - Beta(n)).^2)*prod(K2^2 + h.*(Y - Beta(n)).^2);
end

[~,I] = min(Costo); %Minimo de la funcion costo recursivo

%Proceso de minimizaci�n. Estimaci�n...
den1 = k12 + (X-Beta(I)).^2; %Denominador de elementos no recursivos
den2 = k22 + (Y-Beta(I)).^2; %Denominador de elementos recursivos
band = Beta(I); %Primera estimaci�n
bd = (sum(X./den1) + sum(Y./den2))/(sum(1./den1) + sum(1./den2));

while abs(bd - band) > 0.000001 %Iteraciones, Proceso de minimizaci�n. Estimaci�n...
    den1 = k12 + (X-bd).^2; %Actualizaci�n de los denominadores.
    den2 = k22 + (Y-bd).^2;
    band = bd;
    bd = (sum(X./den1) + sum(Y./den2))/(sum(1./den1) + sum(1./den2));
end

beta = bd; %Valor de Beta estimado. Salida de la estimaci�n recursiva.
return;