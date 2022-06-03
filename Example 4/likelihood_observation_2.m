function L = likelihood_observation_2(C)
% First Likelihood function
chloride_inf = 0.0057;
% L = normpdf(C-chloride_inf,0,0.0005);
L = log(1/(0.0005*sqrt(2*pi)))-0.5*(C-chloride_inf).^2/0.0005^2;
end