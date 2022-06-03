%% BUS + Subset
Num = 5e3; po = 0.1; Num_seed = Num*po; interval = 1;n = 2;
x_s = [rand(Num ,1) reshape(lhsnorm(0,1,Num*n),Num,n)];

t = 10;
m = 1; fig_num = 1;
while(t>0)
%% Find the t
y_s = G2(x_s); y_sort = sort(y_s); t = y_sort(Num_seed);
x_seed = x_s(y_s<=t,:); x_seed = x_seed(1:Num_seed,:);  mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = [];
num_mcmc = 0;
disp(['Now t',num2str(m),' is:  '  num2str(t)]);

% Draw samples in the first subset
while(num_mcmc<3*Num)   
% estimate the pdf here
pdf_x_seed = unifpdf(x_seed(:,1),0,1);
for i = 2:n+1
pdf_x_seed = pdf_x_seed.*mvnpdf(x_seed(:,i),mu_matrix,Sigma_matrix);
end

x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,n+1); 
% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1);
for i = 2:n+1
pdf_x_proposal = pdf_x_proposal.*mvnpdf(x_proposal(:,i),mu_matrix,Sigma_matrix);
end

y_proposal = G2(x_proposal); I_proposal = zeros(Num_seed,1);I_proposal(y_proposal<=t)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc); x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_s = x_mcmc(end-Num+1:end,:);
m = m +1;
figure(fig_num);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',8);
fig_num = fig_num + 1;
end

%% Collect the seeds in the final subset
y_s = G2(x_s); x_seed = x_s(y_s<=0,:);  mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = []; num_mcmc = 0; 
while(num_mcmc<3*Num)
    
Num_seed_final = length(x_seed);
% estimate the pdf here
pdf_x_seed = unifpdf(x_seed(:,1),0,1);
for i = 2:n+1
pdf_x_seed = pdf_x_seed.*mvnpdf(x_seed(:,i),mu_matrix,Sigma_matrix);
end

x_proposal = x_seed + unifrnd(-interval,interval,Num_seed_final,n+1); 
% estimate the pdf here
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1);
for i = 2:n+1
pdf_x_proposal = pdf_x_proposal.*mvnpdf(x_proposal(:,i),mu_matrix,Sigma_matrix);
end

y_proposal = G2(x_proposal); I_proposal = zeros(Num_seed_final,1);
I_proposal(y_proposal<=0)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed_final,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_s = x_mcmc(end-Num+1:end,:); 
figure(fig_num);
plot(x_s(:,2),x_s(:,3),'.k','MarkerSize',8);
%% Report the results
e_mu_bus = mean(x_s(:,2)); 
e_std_bus = std(x_s(:,2)); 
disp(['e_mu_bus is:  '  num2str(e_mu_bus)]);
disp(['e_std_bus is:  '  num2str(e_std_bus)]);

ce = 1e-4;
sig_l = 0.2;
a = -2*(1 + sig_l^2);
b = sqrt(2*pi*sqrt(1 + sig_l^2));
c = log(ce^(1/n)*b);
mu_l = sqrt(a*c); 
mu_prime = mu_l/(sig_l^2*(1+1/(sig_l^2)));
sig_prime = sqrt(1/(1+1/(sig_l^2)));
disp(['mu_prime is:  '  num2str(mu_prime)]);
disp(['sig_prime is:  '  num2str(sig_prime)]);
disp(['the ratio of e_mu_bus is:  '  num2str(e_mu_bus/mu_prime)]);
disp(['the ratio of e_std_bus is:  '  num2str(e_std_bus/sig_prime)]);


% figure(4);
% plot(x_s(:,1),x_s(:,2),'.k','MarkerSize',8);
% hold on;
% set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{x_{1}}$$',...
%     'FontName','times','FontSize',15);
% set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{x_{2}}$$',...
%     'FontName','times','FontSize',15);
% set(gca,'fontsize',15);
% set(get(gca,'legend'),'Interpreter','latex',...
%     'FontName','times','FontSize',12);
% axis([1.5 3.5 1.5 3.5]);

