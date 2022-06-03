Num_initial = 10;
T = [T0;T1;T2;T3;]; 


true_T0 = t0_true*ones(length(T0),1); true_T1 = t1_true*ones(length(T1),1); true_T2 = t2_true*ones(length(T2),1);
true_T3 = t3_true*ones(length(T3),1); 

T_true = [true_T0;true_T1;true_T2;true_T3;]; 

figure(13);
n_call = Num_initial+1:Num_initial+length(T);
plot(n_call,T,'--k','LineWidth',1.5);
hold on;
plot(n_call,T_true,'-k','LineWidth',0.75);
hold on;



grid on;
hold on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{N_{call}}$$',...
    'FontName','times','FontSize',15)
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{\hat{t}_{i}}$$',...
    'FontName','times','FontSize',15)
set(gca,'fontsize',15);

legend('$$\it{\hat{t}_{i}}$$','$$\it{t_{i}}$$');
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);

legend boxoff;
axis([10 35 -10 1.1*max(T)]);