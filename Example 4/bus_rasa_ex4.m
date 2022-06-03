%% Draw the samples for MCS
Num = 5e3; p_i = 0.1; pd=makedist('Normal'); n_dim = 5;
x_s = generate_samples_ex4(Num);

%% Find the value t0
Num_initial = 30;
% Use Latin Hypercube Sampling to generate random variable samples
x_i = x_s(1:Num_initial,:); 
x_tr = x_i;

%% Use traditional kriging model 
% use the kernel function to define the correlation matrix 
lob = 1e-3*ones(1,n_dim); upb = 10*ones(1,n_dim); theta = 1e-1*ones(1,n_dim);
m = 1; ratio = 1; T1 = [];
% indicator_t = []; cov_t = 1; num_det = 10;
while(ratio>1e-2) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G4(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta);
% [dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = 15:0.01:25;
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
% % second stopping criterion
% indicator_t = [indicator_t; t1];
% if m>=num_det
%     std_t = std(indicator_t(end-num_det+1:end));
%     mean_t = mean(indicator_t(end-num_det+1:end));
%     cov_t = std_t/mean_t;
% end

% Use the 'U' learning function  
U = abs(y_kr-t1)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Now ratio is:  '  num2str(ratio)]);
disp(['The value of ''t1'' is:  '  num2str(t1)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);
 
 m = m + 1;
end

y_mcs = G4(x_s);
true_t1 = true_find_t(p_i,y_mcs);

% figure(1);
% n_call = (length(x_tr)-length(T1)+1):length(x_tr); true_T1 = true_t1*ones(2*size(x_tr,1),1);
% semilogy(n_call,T1,'-k','LineWidth',1.5,'MarkerSize',5);
% hold on;
% semilogy(true_T1,'--r','LineWidth',1.5,'MarkerSize',4);
% grid on;
% hold on;
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{N_{call}}$$',...
%     'FontName','times','FontSize',15)
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{\hat{t}_{1}}$$',...
%     'FontName','times','FontSize',15)
% set(gca,'fontsize',15);
% 
% legend('$$\it{\hat{t}_{1}}$$','$$\it{t_{1}}$$');
% set(get(gca,'legend'),'Interpreter','latex',...
%     'FontName','times','FontSize',15);
% axis([30 250 0.016 0.0175]);
% legend boxoff;

%% Draw the samples from the first subset
x_s = [];
while(length(x_s)<Num)
x_sub = generate_samples_ex4(1e3);
[y_kr, ~] = predictor(x_sub, dmodel);
x_sub = x_sub(y_kr<=t1,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RASA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the value t2
% num_pre = round(size(x_tr,1)*0.8);
% x_tr = x_i; 
T2 = [];
% use the kernel function to define the correlation matrix 
% lob = 1e-4*ones(1,n_dim); upb = 10*ones(1,n_dim); theta = 1e-4*ones(1,n_dim);
m = 1; ratio = 1;
% indicator_t = []; cov_t = 1; 
while(ratio>1e-2) 
% while(ratio>1e-2&&cov_t>1e-20) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G4(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta);
% [dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
[y_kr, mse_kr] = predictor(x_s, dmodel);

% Estimate the 't' value
lsf_line = 0:0.01:5;
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
% % second stopping criterion
% indicator_t = [indicator_t; t2];
% if m>=num_det
%     std_t = std(indicator_t(end-num_det+1:end));
%     mean_t = mean(indicator_t(end-num_det+1:end));
%     cov_t = std_t/mean_t;
% end

% Use the 'U' learning function  
U = abs(y_kr-t2)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Now ratio is:  '  num2str(ratio)]);
disp(['The value of ''t2'' is:  '  num2str(t2)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);

 m = m + 1;
end

y_mcs = G4(x_s);
true_t2 = true_find_t(p_i,y_mcs);


%% Draw the samples from the second subset
x_s = [];
while(length(x_s)<Num)
x_sub = generate_samples_ex4(1e3);
[y_kr, ~] = predictor(x_sub, dmodel);
x_sub = x_sub(y_kr<=t2,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RASA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    3    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find the value t2
% num_pre = round(size(x_tr,1)*0.8);
%% Find the value t3
% x_tr = x_i;
% use the kernel function to define the correlation matrix 
% lob = 1e-3*ones(1,n_dim); upb = 10*ones(1,n_dim); theta = 1e-3*ones(1,n_dim);
m = 1; ratio = 1; T3 = [];
% indicator_t = []; cov_t = 1;
while(ratio>1e-2) 
% while(ratio>1e-2&&cov_t>1e-20) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G4(x_tr);

%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta);
% [dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
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
t3 = lsf_line(pf_abs==pf_min);
T3 = [T3; t3;];
% % second stopping criterion
% indicator_t = [indicator_t; t3];
% if m>=num_det
%     std_t = std(indicator_t(end-num_det+1:end));
%     mean_t = mean(indicator_t(end-num_det+1:end));
%     cov_t = std_t/mean_t;
% end

% Use the 'U' learning function  
U = abs(y_kr-t3)./sqrt(mse_kr);

%% Error control
stopping_pf = pf_pc(pf_abs==pf_min);
stopping_sd = pf_sd(pf_abs==pf_min);
ratio = stopping_sd/stopping_pf;

% Find the new points
disp(['Now m is:  '  num2str(m)]);
disp(['Now ratio is:  '  num2str(ratio)]);
disp(['The value of ''t3'' is:  '  num2str(t3)]);
P_new = find(U==min(U));
x_new = x_s(P_new,:);

 m = m + 1;
end

y_mcs = G4(x_s);
true_t3 = true_find_t(p_i,y_mcs);

% figure(3);
% n_call = (length(x_tr)-length(T3)+1):length(x_tr); true_T3 = true_t3*ones(2*size(x_tr,1),1);
% plot(n_call,T3,'-k','LineWidth',1.5,'MarkerSize',5);
% hold on;
% plot(true_T3,'--r','LineWidth',1.5,'MarkerSize',4);
% grid on;
% hold on;
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{N_{call}}$$',...
%     'FontName','times','FontSize',15)
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{\hat{t}_{3}}$$',...
%     'FontName','times','FontSize',15)
% set(gca,'fontsize',15);
% 
% legend('$$\it{\hat{t}_{3}}$$','$$\it{t_{3}}$$');
% set(get(gca,'legend'),'Interpreter','latex',...
%     'FontName','times','FontSize',15);
% axis([20 450 -0.005 -0.002]);
% legend boxoff;

%% Estimate the probability of failure in the last subset
% x_tr = x_i;
% use the kernel function to define the correlation matrix 
% lob = 1e-6*ones(1,n_dim); upb = 10*ones(1,n_dim); theta = 1e-6*ones(1,n_dim);
m = 1; ratio = 1;
% indicator_pf = []; std_pf = 1; 
while(ratio>1e-2) 
% while(ratio>1e-2&&cov_t>1e-20) 

if m>1
   x_tr = [x_tr;x_new;];
end

y_tr = G4(x_tr);
%% Dace tool box
[dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta);
% [dmodel, perf] = dacefit(x_tr, y_tr, @regpoly0, @corrgauss, theta,lob,upb);
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



% % second stopping criterion
% indicator_pf = [indicator_pf; p_f];
% if m>=num_det
%     std_p_f = std(indicator_pf(end-num_det+1:end));
% end

 m = m + 1;
end

%% Estimate the failure probability
U_pc = -y_kr./sqrt(mse_kr);
I_pc = cdf(pd,U_pc);
pf_pc = sum(I_pc)/Num;
p_f = pf_pc*1e-2;
disp(['Probability of failure by MCS:  '  num2str(p_f)]);



