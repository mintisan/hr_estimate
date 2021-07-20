function y = MGCFilter(x, p, sigma)

N = length(x);

for ii = 1:N
    C(ii) = sum(log(sigma^p + abs(x - x(ii)).^p));
end

[minvaule, index] = min(C);
y = x(index);
