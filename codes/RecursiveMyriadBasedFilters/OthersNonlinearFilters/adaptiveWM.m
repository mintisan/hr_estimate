function [g,h,e] = adaptiveWM(X,d,g,h,u0)

N  = length(X);
M2  = length(h);
m = 1;

for n = M2 : N - M2
    u = u0*exp(-(m-1)/1000);
    
    windowg = X(n - M2 + 1:n + M2);
    windowd = d(n - M2 + 1:n);

    x = sign(g).*windowg;
    z = sign(h).*windowd;
    ag = abs(g);
    ah = abs(h);
    xv = [x z];
    wv = [ag ah];
    
    Y = wmedian(xv,wv);     
    
    e = d(n+1) - Y;
    
    g = g + u * e * sign(g) .* sign(x - Y);
    h = h + u * e * sign(h) .* sign(z - Y);
    
    m = m + 1;
end

return;