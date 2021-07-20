function [beta] = RWMyFfpsII(X,Y,g,h,K1,K2)

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


xt = repmat(X,T,1);
yt = repmat(Y,T,1);
beta1 = repmat(Beta',1,length(g));
beta2 = repmat(Beta',1,length(h));

cst1 = K1^2 * ones(T,length(g)) + repmat(g,T,1) .* (xt - beta1).^2;
cst2 = K2^2 * ones(T,length(h)) + repmat(h,T,1) .* (yt - beta2).^2;
Costo = prod([cst1 cst2],2);

[Q,I] = min(Costo); %Minimo de la funcion costo recursivo

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