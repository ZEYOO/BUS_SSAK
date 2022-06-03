
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw the samples for MCS
Num = 5e3; p_i = 0.1; pd=makedist('Normal'); n = 2; interval = 0.075;
mu1 = 0.5;sigma1 = 0.5; mu2 = 0.16;sigma2 = 0.63;
x1 = lognrnd(mu1,sigma1,Num,1); x2 = lognrnd(mu2,sigma2,Num,1);
p = rand(Num,1); x_s = [p x1 x2];

% figure(1);
% plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',8);

%% Find the value t0
Num_initial = 10;
% Use Latin Hypercube Sampling to generate random variable samples
x_i = x_s(1:Num_initial,:);
x_tr = x_i;

figure(1);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',6);
axis([0 3 0 1.5]);
grid on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);

figure(2);
plot(x_tr(:,2),x_tr(:,3),'*k','MarkerSize',8);
axis([0 3 0 1.5]);
grid on;
legend('$$\it{x_{in}}$$');
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);

x_tr_old = x_tr;
size_old = length(x_tr_old);

%% Use traditional kriging model 
% use the kernel function to define the correlation matrix 
lob = 1e-3*ones(1,n+1);upb = 10*ones(1,n+1);theta = ones(1,n+1);
m = 1;ratio = 1;
T0 = [];

while(ratio>=1e-2) 
    
% while(length(x_tr)<30) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G3(x_tr);
%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = 20:0.1:60;
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
disp(['Now ratio is:  '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end

%% Use the same sample to do monte carlo simulation
y_mcs = G3(x_s); y_sort = sort(y_mcs); t0_true = y_sort(0.1*length(y_sort));
len_x = length(x_s);
I = zeros(len_x,1);
I(y_mcs-t0<0) = 1;
I_p = sum(I);
P_f_MCS = I_p/len_x;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['The true Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t0'' is:  '  num2str(t0)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%
% Draw samples in the first subset
y_ss = predictor(x_s, dmodel); x_seed = x_s(y_ss<=t0,:); Num_seed = length(x_seed); 
x_seed = x_seed(1:Num_seed,:); mu_matrix = 0; Sigma_matrix = 1; x_mcmc = []; num_mcmc = 0;

while(num_mcmc<3*Num)    
    
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*lognpdf(x_seed(:,2),mu1,sigma1).*lognpdf(x_seed(:,3),mu2,sigma2);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,n+1); 

% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*lognpdf(x_proposal(:,2),mu1,sigma1).*lognpdf(x_proposal(:,3),mu2,sigma2);

y_proposal = predictor(x_proposal, dmodel); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=t0)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_s = x_mcmc(end-Num+1:end,:);

figure(3);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',6);
axis([0 3 0 1.5]);
grid on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);

figure(4);
plot(x_tr(size_old+1:end,2),x_tr(size_old+1:end,3),'*k','MarkerSize',8);
hold on;
plot(x_tr_old(:,2),x_tr_old(:,3),'ok','MarkerSize',8);
axis([0 3 0 1.5]);
grid on;
legend('$$\it{New\ x_{tr}}$$','$$\it{Old\ x_{tr}}$$');
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);


x_tr_old = x_tr;
size_old = length(x_tr_old);

%% Find the value t1
% x_tr = x_i; 
T1 = [];
% use the kernel function to define the correlation matrix 
lob = 1e-3*ones(1,n+1);upb = 10*ones(1,n+1);theta = ones(1,n+1); m = 1; ratio = 1;
while(ratio>=1e-2) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G3(x_tr);
%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = -5:0.1:10;
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
disp(['Now ratio is:  '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end

y_mcs = G3(x_s); y_sort = sort(y_mcs); t1_true = y_sort(0.1*length(y_sort));
len_x = length(x_s);
I = zeros(len_x,1);
I(y_mcs-t1<0) = 1;
I_p = sum(I);
P_f_MCS = I_p/len_x;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['The true Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t1'' is:  '  num2str(t1)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%
% Draw the samples from the second subset with MCS
y_ss = predictor(x_s, dmodel); x_seed = x_s(y_ss<=t1,:); Num_seed = length(x_seed); 
x_seed = x_seed(1:Num_seed,:); mu_matrix = 0; Sigma_matrix = 1; x_mcmc = []; num_mcmc = 0;

while(num_mcmc<3*Num) 
    
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*lognpdf(x_seed(:,2),mu1,sigma1).*lognpdf(x_seed(:,3),mu2,sigma2);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,n+1); 

% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*lognpdf(x_proposal(:,2),mu1,sigma1).*lognpdf(x_proposal(:,3),mu2,sigma2);
y_proposal = predictor(x_proposal, dmodel); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=t1)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_s = x_mcmc(end-Num+1:end,:);

figure(5);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',6);
axis([0 3 0 1.5]);
grid on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);

figure(6);
plot(x_tr(size_old+1:end,2),x_tr(size_old+1:end,3),'*k','MarkerSize',8);
hold on;
plot(x_tr_old(:,2),x_tr_old(:,3),'ok','MarkerSize',8);
axis([0 3 0 1.5]);
grid on;
legend('$$\it{New\ x_{tr}}$$','$$\it{Old\ x_{tr}}$$');
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);


x_tr_old = x_tr;
size_old = length(x_tr_old);


%% Find the value t2
% x_tr = x_i; 
T2 = [];
% use the kernel function to define the correlation matrix 
lob = 1e-3*ones(1,n+1);upb = 10*ones(1,n+1);theta = ones(1,n+1); m = 1; ratio = 1;
while(ratio>=1e-2) 
if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G3(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = -2:0.01:3;
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
disp(['Now ratio is:  '  num2str(ratio)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
m = m + 1;
end

y_mcs = G3(x_s); y_sort = sort(y_mcs); t2_true = y_sort(0.1*length(y_sort));
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
% x_tr = x_i;
% use the kernel function to define the correlation matrix 
lob = 1e-3*ones(1,n+1);upb = 10*ones(1,n+1);theta = 1e-2*ones(1,n+1); m = 1; ratio = 1;
while(ratio>=1e-2) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G3(x_tr);

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
 m = m + 1;
end



interval = 0.5;
%% Collect the seeds in the final subset
% Draw the samples from the second subset with MCS
y_ss = predictor(x_s, dmodel); x_seed = x_s(y_ss<=0,:); Num_seed = length(x_seed); 
x_seed = x_seed(1:Num_seed,:); mu_matrix = 0; Sigma_matrix = 1; x_mcmc = []; num_mcmc = 0;

while(num_mcmc<3*Num) 
    
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*lognpdf(x_seed(:,2),mu1,sigma1).*lognpdf(x_seed(:,3),mu2,sigma2);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,n+1); 

% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*lognpdf(x_proposal(:,2),mu1,sigma1).*lognpdf(x_proposal(:,3),mu2,sigma2);
y_proposal = predictor(x_proposal, dmodel); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=0)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_s = x_mcmc(end-Num+1:end,:);

figure(7);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',6);
axis([0 3 0 1.5]);
grid on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);

figure(8);
plot(x_tr(size_old+1:end,2),x_tr(size_old+1:end,3),'*k','MarkerSize',8);
hold on;
plot(x_tr_old(:,2),x_tr_old(:,3),'ok','MarkerSize',8);
axis([0 3 0 1.5]);
grid on;
legend('$$\it{New\ x_{tr}}$$','$$\it{Old\ x_{tr}}$$');
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',12);
x_tr_old = x_tr;
size_old = length(x_tr_old);


%% Report the results
% Estimate the mean and standard deviation
x_posterior_est = [x_s(:,2) x_s(:,3)];

e_mu_L = mean(x_posterior_est(x_posterior_est(:,1)<1,1)); 
e_std_L = std(x_posterior_est(x_posterior_est(:,1)<1,1)); 

disp(['e_mu_L is:  '  num2str(e_mu_L )]);
disp(['e_std_L is:  '  num2str(e_std_L)]);

e_mu_R = mean(x_posterior_est(x_posterior_est(:,1)>1,1)); 
e_std_R = std(x_posterior_est(x_posterior_est(:,1)>1,1)); 

disp(['e_mu_R is:  '  num2str(e_mu_R )]);
disp(['e_std_R is:  '  num2str(e_std_R)]);
% 
% figure(5);
% plot(x_posterior_est(:,1),x_posterior_est(:,2),'xk');
% hold on;
% axis([0 3 0 1.5]);
% 
% legend('$$\it{Accepted\ samples}$$');
% 
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
%     'FontName','times','FontSize',15);
% 
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
%     'FontName','times','FontSize',15);
% 
% set(gca,'fontsize',15);
% 
% set(get(gca,'legend'),'Interpreter','latex',...
%     'FontName','times','FontSize',15);
% legend boxoff;

