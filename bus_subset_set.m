%% Implement MCMC based on true performance
Num = 1e4; po = 0.1; Num_seed = Num*po; interval = 1.2;
x_ss = [rand(Num ,1) lhsnorm(0,1,Num)];


% Find the t1
y_ss = lsf_1(x_ss); y_sort = sort(y_ss); 
t1 = y_sort(Num_seed);
x_seed = x_ss(y_ss<=t1,:); x_seed = x_seed(1:Num_seed,:);  mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = [];
num_mcmc = 0;

disp(['Now t1 is:  '  num2str(t1)]);

% Draw samples in the first subset
while(num_mcmc<3*Num)    
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*mvnpdf(x_seed(:,2),mu_matrix,Sigma_matrix);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,2); 
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*mvnpdf(x_proposal(:,2),mu_matrix,Sigma_matrix);
y_proposal = lsf_1(x_proposal); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=t1)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_ss = x_mcmc(end-Num+1:end,:);

% Find the t2
y_ss = lsf_1(x_ss); y_sort = sort(y_ss); 
t2 = y_sort(Num_seed);
x_seed = x_ss(y_ss<=t2,:); x_seed = x_seed(1:Num_seed,:); mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = [];
num_mcmc = 0;

disp(['Now t2 is:  '  num2str(t2)]);

% Draw samples in the second subset
while(num_mcmc<3*Num)    
pdf_x_seed = unifpdf(x_seed(:,1),0,1).*mvnpdf(x_seed(:,2),mu_matrix,Sigma_matrix);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,2); 
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*mvnpdf(x_proposal(:,2),mu_matrix,Sigma_matrix);
y_proposal = lsf_1(x_proposal); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=t2)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end

x_ss = x_mcmc(end-Num+1:end,:);

% Find the t3
y_ss = lsf_1(x_ss); y_sort = sort(y_ss); 
t3 = y_sort(Num_seed);
x_seed = x_ss(y_ss<=t3,:); x_seed = x_seed(1:Num_seed,:); mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = [];
num_mcmc = 0;

disp(['Now t3 is:  '  num2str(t3)]);

% Draw samples in the third subset
while(num_mcmc<3*Num)    
pdf_x_seed = mvnpdf(x_seed(:,2),mu_matrix,Sigma_matrix);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed,2); 
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*mvnpdf(x_proposal(:,2),mu_matrix,Sigma_matrix);
y_proposal = lsf_1(x_proposal); I_proposal = zeros(Num_seed,1); I_proposal(y_proposal<=t3)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end
x_ss = x_mcmc(end-Num+1:end,:);

% Collect the seeds in the final subset
y_ss = lsf_1(x_ss); x_seed = x_ss(y_ss<=0,:);  mu_matrix = 0; Sigma_matrix = 1;  x_mcmc = [];
num_mcmc = 0; 

while(num_mcmc<3*Num)
Num_seed_final = length(y_ss);
pdf_x_seed = mvnpdf(x_seed(:,2),mu_matrix,Sigma_matrix);
x_proposal = x_seed + unifrnd(-interval,interval,Num_seed_final,2); 
pdf_x_proposal = unifpdf(x_proposal(:,1),0,1).*mvnpdf(x_proposal(:,2),mu_matrix,Sigma_matrix);
y_proposal = lsf_1(x_proposal); I_proposal = zeros(Num_seed_final,1); I_proposal(y_proposal<=0)=1;
r_proposal = I_proposal.*pdf_x_proposal./pdf_x_seed; r_uni = unifrnd(0,1,Num_seed_final,1);
x_mcmc = [x_mcmc;x_proposal(r_proposal>r_uni,:)];
[num_mcmc,~]= size(x_mcmc);
x_seed(r_proposal>r_uni,:)=x_proposal(r_proposal>r_uni,:);
end
x_ss = x_mcmc(end-Num+1:end,:);

e_mu_bus = mean(x_ss(:,2)); 
e_std_bus = std(x_ss(:,2)); 

disp(['e_mu_bus is:  '  num2str(e_mu_bus)]);
disp(['e_std_bus is:  '  num2str(e_std_bus)]);
disp(['the ratio of e_mu_bus is:  '  num2str(e_mu_bus/2.75)]);
disp(['the ratio of e_std_bus is:  '  num2str(e_std_bus/0.287)]);

