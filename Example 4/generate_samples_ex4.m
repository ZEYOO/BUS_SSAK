function x_s = generate_samples_ex2(Num)
% First random variable
p = rand(Num ,1);

% R = evrnd(6.5e4,6.5e3,Num,1);
m = 1.52; coeff = 0.1; std = m*coeff; v = std^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
Cs = lognrnd(mu,sigma,Num,1);

m = 2.5e-4; coeff = 0.1; std = m*coeff; v = std^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
D1 = lognrnd(mu,sigma,Num,1);

m = 1.5e-4; coeff = 0.1; std = m*coeff; v = std^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
D2 = lognrnd(mu,sigma,Num,1);

m = 0.004; coeff = 0.1; std = m*coeff; v = std^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
k = lognrnd(mu,sigma,Num,1);

x_s = [p, Cs, D1, D2, k];
end