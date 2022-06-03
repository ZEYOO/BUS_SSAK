
[x,y] = meshgrid(0:0.01:5,0:0.01:5);
mode = exp(x-y.^2);
variance = (exp(y.^2)-1).*exp(2*x+y.^2);
point  = [x(abs(mode-0.8)<2e-2&abs(variance-1)<1e-2) y(abs(mode-0.8)<2e-2&abs(variance-1)<1e-2)];
% figure(1);
% surf(x,y,mode);

% x = 0:0.01:5;
% y = 0:0.01:5;

