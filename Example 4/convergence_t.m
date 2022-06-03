Num_initial = 12;
T = [T1;T2;T3]; true_T1 = t1*ones(length(T1),1); true_T2 = t2*ones(length(T2),1); true_T3 = t3*ones(length(T3),1);
T_true = [true_T1;true_T2;true_T3]; 

figure(1);
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
axis([0.8*Num_initial 1.1*length(T) 2*min(T) 1.1*max(T)]);
