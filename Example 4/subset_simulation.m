%% Draw the samples for MCS
Num = 5e3; p_i = 0.1;
x_s = generate_samples_ex4(Num);

%% Find the value t1
y_mcs = G4(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t1 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t1<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t1'' is:  '  num2str(t1)]);

%% Draw the samples from the first subset
x_s = [];
while(length(x_s)<Num)
x_sub = generate_samples_ex4(1e5);
y_mcs = G4(x_sub); 
x_sub = x_sub(y_mcs<=t1,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);

%% Find the value t2
y_mcs = G4(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t2 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t2<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t2'' is:  '  num2str(t2)]);

%% Draw the samples from the second subset
x_s = [];
while(length(x_s)<Num)
x_sub = generate_samples_ex4(1e5);
y_mcs = G4(x_sub); 
x_sub = x_sub(y_mcs<=t2,:);
x_s = [x_s;x_sub;];
disp(['Number of points in subset: '  num2str(length(x_s))]);
end
x_s = x_s(1:Num,:);

%% Find the value t3
y_mcs = G4(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
t3 = y_sort(t_index);
I = zeros(Num,1);
I(y_mcs-t3<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
disp(['Number of design points:  '  num2str(length(x_s))]);
disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
disp(['The value of ''t3'' is:  '  num2str(t3)]);

% %% Draw the samples from the second subset
% x_s = [];
% while(length(x_s)<Num)
% x_sub = generate_samples_ex2(1e5);
% y_mcs = G2(x_sub); 
% x_sub = x_sub(y_mcs<=t3,:);
% x_s = [x_s;x_sub;];
% disp(['Number of points in subset: '  num2str(length(x_s))]);
% end
% x_s = x_s(1:Num,:);


% %% Find the value t4
% y_mcs = G2(x_s); y_sort = sort(y_mcs); t_index = Num*p_i;
% t4 = y_sort(t_index);
% I = zeros(Num,1);
% I(y_mcs-t4<=0) = 1;
% I_p = sum(I);
% P_f_MCS = I_p/Num;
% disp(['Number of design points:  '  num2str(length(x_s))]);
% disp(['Probability of failure by MCS:  '  num2str(P_f_MCS)]);
% disp(['The value of ''t4'' is:  '  num2str(t4)]);
% 

%% Estimate the failure probability
I = zeros(Num,1);
I(y_mcs<=0) = 1;
I_p = sum(I);
P_f_MCS = I_p/Num;
p_f = P_f_MCS*1e-2;
disp(['Probability of failure by MCS:  '  num2str(p_f)]);

