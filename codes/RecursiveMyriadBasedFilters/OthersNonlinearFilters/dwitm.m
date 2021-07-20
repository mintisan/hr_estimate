function Y = dwitm(X, W)

X0 = X(:);
W0 = W(:);

IndN    = find(W0 < 0);
XN      = X0(IndN);
WN      = W0(IndN);

IndP    = find(W0 >= 0);
XP      = X0(IndP);
WP      = W0(IndP);

Y = witm(XP, WP)*sum(abs(WP))/sum(abs(W0)) - witm(XN, WN)*sum(abs(WN))/sum(abs(W0));