Num_initial = 10;
T = [T0;T1;T2]; true_T0 = t0*ones(length(T0),1); true_T1 = t1*ones(length(T1),1); true_T2 = t2*ones(length(T2),1);
T_true = [true_T0;true_T1;true_T2]; 

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