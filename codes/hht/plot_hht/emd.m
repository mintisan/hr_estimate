% https://www.cnblogs.com/hdu-zsk/p/4799470.html
% https://zhuanlan.zhihu.com/p/150706759
%主函数
function imf = emd(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x)
% Func : findpeaks
    x = transpose(x(:));
    imf = [];
    while ~ismonotonic(x)
        x1 = x;
        sd = Inf;
        while (sd > 0.1) || ~isimf(x1)
            s1 = getspline(x1);
            s2 = -getspline(-x1);
            x2 = x1-(s1+s2)/2;
       
            sd = sum((x1-x2).^2)/sum(x1.^2);
            x1 = x2;
        end
    
        imf{end+1} = x1;
        x = x-x1;
    end
    imf{end+1} = x;
end

% FUNCTIONS

%非主函数，被调用<br>%判断x是否单调，返回0代表不是单调，返回1代表是单调
function u = ismonotonic(x)
    u1 = length(findpeaks(x))*length(findpeaks(-x));
    if u1 > 0
        u = 0;
    else
        u = 1; 
    end
end

%非主函数，被调用。判断当前x是不是真IMF
function u = isimf(x)
    N  = length(x);
    u1 = sum(x(1:N-1).*x(2:N) < 0);%求x与y=0轴交点的个数
    u2 = length(findpeaks(x))+length(findpeaks(-x));%求极值点个数
    if abs(u1-u2) > 1
        u = 0;
    else
        u = 1;
    end
end

%非主函数，被调用，作用是获得x的包络线
function s = getspline(x)
    N = length(x);
    p = findpeaks(x);
    s = spline([0 p N+1],[0 x(p) 0],1:N);
end
