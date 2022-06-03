function y = G3(x)
 p = x(:,1);
[f1,f2] = two_DOF(x(:,2:3));

% First Likelihood function
f1_p = 3.13;
sd_ep = 1/16;
J1 = (f1.^2/f1_p^2-1).^2;
L1 = -J1/(2*sd_ep^2);

f2_p = 9.83;
sd_ep = 1/16;
J2 = (f2.^2/f2_p^2-1).^2;
L2 = -J2/(2*sd_ep^2);

y = log(p)-L1-L2;
end