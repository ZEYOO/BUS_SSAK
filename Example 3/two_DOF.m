function [f1,f2] = two_DOF(x)
kn = 29.7e6;
m1 = 16.531e3;
m2 = 16.131e3;
k1 = kn*x(:,1);
k2 = kn*x(:,2);
a = m1*m2;
b = -((k1+k2)*m2+m1*k2);
c = k1.*k2;
omeg_sqr_2 = (b-sqrt(b.^2-4*a.*c))./(-2*a);
omeg_sqr_1 = (b+sqrt(b.^2-4*a.*c))./(-2*a);
f2 = sqrt(omeg_sqr_2)/(2*pi);
f1 = sqrt(omeg_sqr_1)/(2*pi);
end