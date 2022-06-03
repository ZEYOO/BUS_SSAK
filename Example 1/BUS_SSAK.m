
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw the samples for MCS
Mu = 0; Sigma = 1; Num = 5e3; p_i = 0.1; pd=makedist('Normal');
% Use Latin Hypercube Sampling to generate random variable samples
x_s = [rand(Num ,1) lhsnorm(0,1,Num)];

%% Find the value t0
Num_initial = 10;
% Use Latin Hypercube Sampling to generate random variable samples
x_i = [rand(Num_initial,1) lhsnorm(0,1,Num_initial)];

x_tr = x_i;y_kr = zeros(length(x_s ),1);mse_kr = zeros(length(x_s ),1);

%% Use traditional kriging model 
% use the kernel function to define the correlation matrix 
lob = [0.001 0.001];upb = [1e3 1e3];theta = [1 1];
m = 1; ratio = 1;
T0 = [];
while(ratio>=1e-2) 
% while(length(x_tr)<60) 
if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G1(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = 0.1:0.1:30;
matrix_y_kr = repmat(y_kr,1,length(lsf_line));
matrix_mse_kr = repmat(mse_kr,1,length(lsf_line));
matrix_lsf_line = repmat(lsf_line,Num,1);
U_pc = -(matrix_y_kr-matrix_lsf_line)./sqrt(matrix_mse_kr);
I_pc = cdf(pd,U_pc);
var_pc = I_pc.*(1-I_pc);

pf_pc = sum(I_pc,1)/Num;
pf_sd = sqrt(sum(var_pc,1))/Num;

pf_up = pf_pc + 1.96*pf_sd;
pf_low = pf_pc - 1.96*pf_sd;

pf_abs = abs(pf_pc-p_i);
pf_min = min(pf_abs);
t0 = lsf_line(pf_abs==pf_min);
T0 = [T0; t0;];
% Use the 'U' learning function  
U = abs(y_kr-t0)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Ratio is: '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end
% 
% figure;
% plot(lsf_line,pf_pc,'k','LineWidth',1,'MarkerSize',5);
% hold on;
% plot(lsf_line,pf_up,'r','LineWidth',1,'MarkerSize',5);
% hold on;
% plot(lsf_line,pf_low,'r','LineWidth',1,'MarkerSize',5);
% 
% % grid on;
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{\hat{t}_{1}}$$',...
%     'FontName','times','FontSize',15)
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{\hat{P}_{o}}$$',...
%     'FontName','times','FontSize',15)
% set(gca,'fontsize',15);
% 
% legend('$$\it{\hat{P}_{o}}$$','$$\it{95\%\ CI\ of\hat{P}_{o}}$$');
% set(get(gca,'legend'),'Interpreter','latex',...
%     'FontName','times','FontSize',15);
% legend boxoff;
% % axis([1.2 1.4  0.05 0.11]);
% axis([1.2 1.3  0.085 0.11]);


%% Use the same sample to do monte carlo simulation
y_mcs = G1(x_s);
len_x = length(x_s);
I = zeros(len_x,1);
I(y_mcs-t0<0) = 1;
I_p = sum(I);
P_f_MCS = I_p/len_x;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['The true Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t0'' is:  '  num2str(t0)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw the samples from the first subset with MCS
x_s = [];
while(length(x_s)<Num)
x_sub = [rand(Num ,1) lhsnorm(0,1,Num)];
[y_kr, ~] = predictor(x_sub, dmodel);
x_sub = x_sub(y_kr<=t0,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);



%% Find the value t1
x_tr = x_i;
T1 = [];
% use the kernel function to define the correlation matrix 
lob = [0.001 0.001];upb = [1e3 1e3];theta = [1 1]; m = 1; ratio = 1;

while(ratio>=1e-2) 
% while(length(x_tr)<60) 
if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G1(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = 0.01:0.01:5;
matrix_y_kr = repmat(y_kr,1,length(lsf_line));
matrix_mse_kr = repmat(mse_kr,1,length(lsf_line));
matrix_lsf_line = repmat(lsf_line,Num,1);
U_pc = -(matrix_y_kr-matrix_lsf_line)./sqrt(matrix_mse_kr);
I_pc = cdf(pd,U_pc);
var_pc = I_pc.*(1-I_pc);


pf_pc = sum(I_pc,1)/Num;
pf_sd = sqrt(sum(var_pc,1))/Num;

pf_up = pf_pc + 1.96*pf_sd;
pf_low = pf_pc - 1.96*pf_sd;

pf_abs = abs(pf_pc-p_i);
pf_min = min(pf_abs);
t1 = lsf_line(pf_abs==pf_min);
T1 = [T1; t1;];
% Use the 'U' learning function  
U = abs(y_kr-t1)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Ratio is: '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end

% figure(2);
% plot(lsf_line,pf_pc,'r');
% hold on;
% plot(lsf_line,pf_up,'b');
% hold on;
% plot(lsf_line,pf_low,'b');

y_mcs = G1(x_s);
len_x = length(x_s);
I = zeros(len_x,1);
I(y_mcs-t1<0) = 1;
I_p = sum(I);
P_f_MCS = I_p/len_x;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['The true Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t1'' is:  '  num2str(t1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw the samples from the second subset with MCS
x_s = [];
while(length(x_s)<Num)
x_sub = [rand(Num ,1) lhsnorm(0,1,Num)];
[y_kr, ~] = predictor(x_sub, dmodel);
x_sub = x_sub(y_kr<=t1,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);


%% Find the value t2
x_tr = x_i; 
T2 = [];
% use the kernel function to define the correlation matrix 
lob = [0.001 0.001];upb = [1e3 1e3];theta = [1 1]; m = 1; ratio = 1;
while(ratio>=1e-2) 
% while(length(x_tr)<60) 
if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G1(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = -5:0.01:5;
matrix_y_kr = repmat(y_kr,1,length(lsf_line));
matrix_mse_kr = repmat(mse_kr,1,length(lsf_line));
matrix_lsf_line = repmat(lsf_line,Num,1);
U_pc = -(matrix_y_kr-matrix_lsf_line)./sqrt(matrix_mse_kr);
I_pc = cdf(pd,U_pc);
var_pc = I_pc.*(1-I_pc);


pf_pc = sum(I_pc,1)/Num;
pf_sd = sqrt(sum(var_pc,1))/Num;

pf_up = pf_pc + 1.96*pf_sd;
pf_low = pf_pc - 1.96*pf_sd;

pf_abs = abs(pf_pc-p_i);
pf_min = min(pf_abs);
t2 = lsf_line(pf_abs==pf_min);
T2 = [T2; t2;];
% Use the 'U' learning function  
U = abs(y_kr-t2)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Ratio is: '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end

% figure;
% plot(lsf_line,pf_pc,'r');
% hold on;
% plot(lsf_line,pf_up,'b');
% hold on;
% plot(lsf_line,pf_low,'b');

y_mcs = G1(x_s);
len_x = length(x_s);
I = zeros(len_x,1);
I(y_mcs-t2<0) = 1;
I_p = sum(I);
P_f_MCS = I_p/len_x;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['The true Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t2'' is:  '  num2str(t2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 4 %%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate the failure probability in last subset simulation
x_tr = x_i;
% use the kernel function to define the correlation matrix 
lob = [0.001 0.001];upb = [1e3 1e3];theta = [1 1]; m = 1; ratio = 1;

% while(ratio>=1e-3) 
while(length(x_tr)<40) 
if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G1(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Use the 'U' learning function  
U = abs(y_kr)./sqrt(mse_kr);

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Now ratio is:  '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);

%% Estimate the failure probability
U_pc = -y_kr./sqrt(mse_kr);
I_pc = cdf(pd,U_pc);
pf_pc = sum(I_pc)/Num;
var_pc = I_pc.*(1-I_pc);
pf_sd = sqrt(sum(var_pc))/Num;
ratio = pf_sd/pf_pc;
p_f = pf_pc*1e-2;
disp(['Probability of failure by MCS:  '  num2str(p_f)]);
 
 m = m + 1;
end

%% Estimate the failure probability
y_mcs = G1(x_s);
I = zeros(Num,1);I(y_mcs<=0) = 1;I_p = sum(I); Pf_mcs = I_p/Num;
Pf_mcs = Pf_mcs*1e-2;
disp(['Probability of failure by Subset:  '  num2str(Pf_mcs)]);


x_f = x_s(y_kr<=0,:);
e_mu_bus = mean(x_f(:,2)); 
e_std_bus = std(x_f(:,2)); 
disp(['e_mu_bus is:  '  num2str(e_mu_bus)]);
disp(['e_std_bus is:  '  num2str(e_std_bus)]);
disp(['the ratio of e_mu_bus is:  '  num2str(e_mu_bus/2.75)]);
disp(['the ratio of e_std_bus is:  '  num2str(e_std_bus/0.287)]);


