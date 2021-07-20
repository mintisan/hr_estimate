function [beta] = HyFfpsII(X,Y,g,h,K)

Beta = [X mean(h.*Y)]; %Vector de puntos de minimizaci�n de la funci�n costo recursivo
              %Beta(1:N) = X; %Elementos no recursivos
              %Beta((N+1):T) = Y; %Elementos recursivos

X = sign(g).*X;
Y = sign(h).*Y;
g = abs(g);
h = abs(h);

N = length(X);%Ancho de la muestra o del vector de observaci�n
M = length(Y); T = N+1;
k1 = K^2./g; %Simplificaci�n matem�tica

xt = repmat(X,T,1);
yt = repmat(Y,T,1);
beta1 = repmat(Beta',1,length(g));
beta2 = repmat(Beta',1,length(h));

cst1 = K^2 * ones(T,length(g)) + repmat(g,T,1) .* (xt - beta1).^2;
cst2 = sum(repmat(h,T,1) .* (yt - beta2).^2, 2);
Costo = prod([cst1 exp(cst2)],2);

[Q,I] = min(Costo); %Minimo de la funci�n costo recursivo

%Proceso de minimizaci�n. Estimaci�n...
den = k1 + (X-Beta(I)).^2; %Denominador de elementos no recursivos
band = Beta(I); %Primera estimaci�n
bd = (sum(X./den) + sum(h.*Y))/(sum(1./den) + sum(h));

while abs(bd - band) > 0.000001 %Iteraciones, Proceso de minimizaci�n. Estimaci�n...
    den = k1 + (X-bd).^2; %Actualizaci�n del denominador.
    band = bd;
    bd = (sum(X./den) + sum(h.*Y))/(sum(1./den) + sum(h));
end

beta = bd; %Valor de Beta estimado. Salida de la estimaci�n recursiva.

return;