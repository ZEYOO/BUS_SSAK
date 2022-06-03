%% Draw the samples for MCS
Mu = 0; Sigma = 1; Num = 5e4; p_i = 0.1;
% Use Latin Hypercube Sampling to generate random variable samples
x_s = [rand(Num ,1) lhsnorm(0,1,Num)];

%% Find the value t1
y_mcs = G1(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t1 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t1<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t1'' is:  '  num2str(t1)]);

%% Draw the samples from the first subset
Num_i = 1e6; x_s = [rand(Num_i,1) lhsnorm(0,1,Num_i)];
y_mcs = G1(x_s); x_s = x_s(y_mcs<=t1,:); x_s = x_s(1:Num,:);


%% Find the value t2
y_mcs = G1(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t2 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t2<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t2'' is:  '  num2str(t2)]);

%% Draw the samples from the second subset
Num_i = 1e7; x_s = [rand(Num_i,1) lhsnorm(0,1,Num_i)];
y_mcs = G1(x_s); x_s = x_s(y_mcs<=t2,:); x_s = x_s(1:Num,:);


%% Find the value t3
y_mcs = G1(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t3 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t3<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t3'' is:  '  num2str(t3)]);

%% Estimate the failure probability
y_mcs = G1(x_s);
I = zeros(Num,1);
I(y_mcs<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
p_f = P_f_MCS*1e-2;
disp(['Probability of failure by MCS:  '  num2str(p_f)]);
