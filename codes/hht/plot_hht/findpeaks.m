%非主函数，被调用
function n = findpeaks(x)%用于寻找极值点，该函数只会求极大值
%   Find peaks.
%   n = findpeaks(x)
    n = find(diff(diff(x)>0)<0);%一阶导数大于0二阶导数小于0的点
    u = find(x(n+1)>x(n));
    n(u) = n(u) + 1;
end