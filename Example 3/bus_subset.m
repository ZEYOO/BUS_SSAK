%% BUS + Subset
Num = 5e3; po = 0.1; Num_seed = Num*po; interval = 0.075;n = 2;

mu1 = 0.5;sigma1 = 0.5;
mu2 = 0.16;sigma2 = 0.63;
x1 = lognrnd(mu1,sigma1,Num,1); x2 = lognrnd(mu2,sigma2,Num,1);
p = rand(Num,1);
x_ss = [p x1 x2];

t = 10;
m = 1;

x_posterior_est = [x_ss(:,2) x_ss(:,3)];
figure(m);
plot(x_ss(:,1),x_ss(:,2),'xk');
hold on;
axis([0 3 0 1.5]);

legend('$$\it{Accepted\ samples}$$');

set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);

set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);

set(gca,'fontsize',15);

set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
legend boxoff;



while(t>0)
%% Find the t
y_ss = G3(x_ss); y_sort = sort(y_ss); 
t = y_sort(Num_seed);
x_seed = x_ss(y_ss<=t,:); x_seed = x_seed(1:Num_seed,:); x_mcmc = [];
num_mcmc = 0;

disp(['Now t',num2str(m),' is:  '  num2str(t)]);

% Draw samples in the first subset
while(num_mcmc<3*Num)   
 
% estimate the pdf here
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*lognpdf(x_seed(:,2),mu1,sigma1).*lognpdf(x_seed(:,3),mu2,sigma2);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,n+1); 

% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*lognpdf(x_proposal(:,2),mu1,sigma1).*lognpdf(x_proposal(:,3),mu2,sigma2);
y_proposal = G3(x_proposal); 

I_proposal = zeros(Num_seed,1);
I_proposal(y_proposal<=t)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_ss = x_mcmc(end-Num+1:end,:);
m = m +1;

x_posterior_est = [x_ss(:,2) x_ss(:,3)];
figure(m);
plot(x_posterior_est(:,1),x_posterior_est(:,2),'xk');
hold on;
axis([0 3 0 1.5]);

legend('$$\it{Accepted\ samples}$$');

set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);

set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);

set(gca,'fontsize',15);

set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
legend boxoff;



end

interval = 0.5;
%% Collect the seeds in the final subset
y_ss = G3(x_ss); x_seed = x_ss(y_ss<=0,:); x_mcmc = [];
num_mcmc = 0; 

while(num_mcmc<3*Num)
Num_seed_final = length(x_seed);

% estimate the pdf here
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*lognpdf(x_seed(:,2),mu1,sigma1).*lognpdf(x_seed(:,3),mu2,sigma2);

x_proposal = x_seed + unifrnd(-interval,interval,Num_seed_final,n+1); 
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*lognpdf(x_proposal(:,2),mu1,sigma1).*lognpdf(x_proposal(:,3),mu2,sigma2);
y_proposal = G3(x_proposal); 


I_proposal = zeros(Num_seed_final,1);
I_proposal(y_proposal<=0)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed_final,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end
x_ss = x_mcmc(end-Num+1:end,:);

%% Report the results
% Estimate the mean and standard deviation
x_posterior_est = [x_ss(:,2) x_ss(:,3)];

e_mu_L = mean(x_posterior_est(x_posterior_est(:,1)<1,1)); 
e_std_L = std(x_posterior_est(x_posterior_est(:,1)<1,1)); 

disp(['e_mu_L is:  '  num2str(e_mu_L )]);
disp(['e_std_L is:  '  num2str(e_std_L)]);

e_mu_R = mean(x_posterior_est(x_posterior_est(:,1)>1,1)); 
e_std_R = std(x_posterior_est(x_posterior_est(:,1)>1,1)); 

disp(['e_mu_R is:  '  num2str(e_mu_R )]);
disp(['e_std_R is:  '  num2str(e_std_R)]);

figure(5);
plot(x_posterior_est(:,1),x_posterior_est(:,2),'xk');
hold on;
axis([0 3 0 1.5]);

legend('$$\it{Accepted\ samples}$$');

set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{X_{1}}$$',...
    'FontName','times','FontSize',15);

set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{X_{2}}$$',...
    'FontName','times','FontSize',15);

set(gca,'fontsize',15);

set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
legend boxoff;

