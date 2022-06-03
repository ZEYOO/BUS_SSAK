function y = G1(X)
p = X(:,1);
x = X(:,2);
mu = 3;
sig = 0.3;
y = log(p)-(-(x-mu).^2/(2*sig^2));
end