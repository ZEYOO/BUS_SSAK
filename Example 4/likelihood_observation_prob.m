function L = likelihood_observation_prob(C)
Num = length(C);
m = 1.823; coeff = 0.05; std = m*coeff; v = std^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
chloride_inf = lognrnd(mu,sigma,Num,1);
L = normpdf(C-chloride_inf,0,0.01);
end