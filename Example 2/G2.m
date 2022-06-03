function y = G2(X)
n = 2;
p = X(:,1);
x = X(:,2:n+1);
ce = 1e-4;
sig_l = 0.2;
a = -2*(1 + sig_l^2);
b = sqrt(2*pi*sqrt(1 + sig_l^2));
c = log(ce^(1/n)*b);
mu_l = sqrt(a*c);
x_hat = (x - mu_l)/sig_l;
y = log(p) + sum(x_hat.^2/2,2);
end