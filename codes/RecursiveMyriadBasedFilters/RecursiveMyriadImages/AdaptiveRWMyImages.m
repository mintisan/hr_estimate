function [g, h, K1, K2, e] = AdaptiveRWMyImages(Io, In, g, h, K1, K2, u0, window_size)

[N1,N2] = size(Io);
mid_pixel = ceil(window_size^2/2);
e = zeros((N2 - window_size + 1)*(N1 - window_size + 1),1);
n = 1;
for j = 1:N2 - window_size + 1
    for i  = 1:N1 - window_size + 1
        Iiw = Io(i:i + window_size - 1,j:j + window_size - 1);
        Iix = Iiw(:);
        windowd = Iix(1:mid_pixel-1)';
        
        Inw = In(i:i + window_size - 1,j:j + window_size - 1);
        windowg = Inw(:)';
        
        u = u0*exp(-(n-1)/1000);
        x = sign(g).*windowg;
        z = sign(h).*windowd;
        ag = abs(g);
        ah = abs(h);
        d  = Iix(mid_pixel);
    
        Y(n) = RWMyFfpsII(windowg,windowd,g,h,K1,K2);
    
        e(n) = Y(n) - d;
    
        num1 = K1^2*(sign(g).*(x - Y(n)))./((K1^2 + ag.*((x - Y(n)).^2)).^2);
        num2 = K2^2*(sign(h).*(z - Y(n)))./((K2^2 + ah.*((z - Y(n)).^2)).^2);
 
        den = sum(ag.*((K1^2 - ag.*((x - Y(n)).^2))./(K1^2 + ag.*((x - Y(n)).^2)).^2)) +...
        sum(ah.*((K2^2 - ah.*((z - Y(n)).^2))./(K2^2 + ah.*((z - Y(n)).^2)).^2));
    
        g = g - u*sign(e(n)).*(num1/den);
        h = h - u*sign(e(n)).*(num2/den);

        n = n + 1;        
    end
end