function y = G4(x_s)
p = x_s(:,1); Cs = x_s(:,2); D1 = x_s(:,3); D2 = x_s(:,4); k = x_s(:,5);
% The first observation
x1 = 0.005; t1 = 5; Co = 0;
[C, dis_carb] = chloride_carbonation(x1,t1,D1,D2,k,Cs,Co);
L1 = likelihood_observation_1(C);
L2 = likelihood_observation_2(dis_carb);
L_prob = likelihood_observation_prob(C);
% y = p -L_prob;
% y = p -L1.*L2;
% y = log(p) -log(L1)-log(L2);
y = log(p) -L1-L2;
end