function L = likelihood_observation_1(C)
% First Likelihood function
chloride_inf = 1.823;
% L = normpdf(C-chloride_inf,0,0.05);
L = log(1/(0.05*sqrt(2*pi)))-0.5*(C-chloride_inf).^2/0.05^2;
end