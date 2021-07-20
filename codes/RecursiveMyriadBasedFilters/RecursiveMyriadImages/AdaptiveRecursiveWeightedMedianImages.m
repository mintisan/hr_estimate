function [g, h, e] = AdaptiveRecursiveWeightedMedianImages(Io, In, g, h, u0, window_size)

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
        xv = [x z];
        wv = [ag ah];
        
        d  = Iix(mid_pixel);
    
        Y = wmedian(xv,wv);
        
        e(n) = d - Y;
        g = g + u * e(n) * sign(g) .* sign(x - Y);
        h = h + u * e(n) * sign(h) .* sign(z - Y);
        n = n + 1;        
    end
end