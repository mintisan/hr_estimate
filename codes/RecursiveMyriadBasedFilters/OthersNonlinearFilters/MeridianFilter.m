function y = MeridianFilter(x, delta)

N = length(x);

for ii = 1:N
    C(ii) = sum(log(delta + abs(x - x(ii))));
end

[minvaule, index] = min(C);
y = x(index);
