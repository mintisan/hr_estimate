function [beta] = WMyFfpsII(X,W,K)

B = X;
X = sign(W).*X;
W = abs(W);   
N = length(X); %Ancho de la muestra o del vector de observacion
k2 = K^2./W; %Simplificacion matem�tica

Costo = zeros(1,N); %Inicializacion de la funcion costo
cst = zeros(N,N);
for n = 1:N %Evaluacion de la funcion costo en las muestras de observacion
    cst(n,:) = K^2 + W.*(X - B(n)).^2;
end

Costo = sum(log(cst));
[Q,I] = min(Costo); %Minimo de la funcion costo

den = k2 + (X-B(I)).^2; %Proceso de minimizaci�n. Estimaci�n...
band = B(I);
bd = sum(X./den)/sum(1./den);

while abs(bd - band) > 0.000001        
    den = k2 + (X-bd).^2;
    band = bd;
    bd = sum(X./den)/sum(1./den);
end

beta = bd; %Valor de Beta estimado. Salida del filtro

return;