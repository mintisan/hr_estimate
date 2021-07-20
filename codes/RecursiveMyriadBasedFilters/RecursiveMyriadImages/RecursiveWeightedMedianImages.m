function [Iout] = RecursiveWeightedMedianImages(Iin, g, h, window_size)

[n_row, n_col]  = size(Iin);
Iout            = zeros(n_row,n_col);
mid_pixel       = ceil(window_size^2/2);
wv              = [abs(g) abs(h)];

for j = 1: n_col - window_size + 1
    for i  = 1: n_row - window_size + 1   
        Isw = Iout(i:i + window_size - 1,j:j + window_size - 1);
        Isx = Isw(:);
        windowh = Isx(1:mid_pixel-1)';
        
        Inw = Iin(i:i + window_size - 1,j:j + window_size - 1);
        windowg = Inw(:)';
        
        xv = [sign(g).*windowg sign(h).*windowh];        
        
        Iout(i + floor(window_size/2),j + floor(window_size/2)) = wmedian(xv, wv);
    end
end

end

