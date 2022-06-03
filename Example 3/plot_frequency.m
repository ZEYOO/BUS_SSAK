[X1,X2] = meshgrid(0:0.1:3,0:0.1:3);
kn = 29.7e6;
m1 = 16.531e3;
m2 = 16.131e3;
k1 = kn*X1;
k2 = kn*X2;

a = m1*m2;
b = -((k1+k2)*m2+m1*k2);
c = k1.*k2;
omeg_sqr_2 = (b-sqrt(b.^2-4*a.*c))./(-2*a);
omeg_sqr_1 = (b+sqrt(b.^2-4*a.*c))./(-2*a);

f2 = sqrt(omeg_sqr_2)/(2*pi);
f1 = sqrt(omeg_sqr_1)/(2*pi);

f1_p = 3.13;
f2_p = 9.83;
sd_ep = 1/16;
J = (f1.^2/f1_p^2-1).^2+(f2.^2/f2_p^2-1).^2;
% J = (f1.^2/f1_p^2-1).^2;
L = exp(-J/(2*sd_ep^2));


figure(1);
surf(X1,X2,f1);

figure(2);
surf(X1,X2,f2);

figure(3);
surf(X1,X2,L);




