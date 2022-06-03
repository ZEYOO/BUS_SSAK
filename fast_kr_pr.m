function [y_kr, mse_kr] = fast_kr_pr(x, dmodel)
thr = 1e4;
m = size(x,1);
p = rem(m,thr);
n = (m - p)/thr;
[y_kr, mse_kr] = predictor(x(1:thr,:),dmodel);
for i = 2:n
[y,mse] = predictor(x(thr*(i-1)+1:thr*i,:),dmodel);
y_kr = [y_kr;y;];
mse_kr = [mse_kr;mse;];
end
end

