% Co = 0;
% Cs = 1.52; D1 = 2.5e-4; D2 = 1.5e-4; k = 0.004;
% [x,t] = meshgrid(0:0.001:0.05,0:0.1:50);
% [C, dis_carb] = chloride_carbonation(x,t,D1,D2,k,Cs,Co);
% surf(x,t,C);

Co = 0;
Cs = 1.52; D1 = 2.5e-4; D2 = 1.5e-4; k = 0.004; t = 10;
x = 0:0.001:0.06;
[C, dis_carb] = chloride_carbonation(x,t,D1,D2,k,Cs,Co);

x_l_b = [-0.001   -1];
x_l_t = [-0.001    5];
x_r_b = [dis_carb -1];
x_r_t = [dis_carb  5];
x_esr = [x_l_b;x_l_t;x_r_t;x_r_b;x_l_b;];

figure(1);
carbonation_fill = fill(x_esr(:,1),x_esr(:,2),[192 192 192]/255,'EdgeColor',[1 1 1]);
hold on;
% surf(x,t,C);
plot(x(2:end),C(2:end),'-k');

grid off;
hold on;
set(get(gca,'XLabel'),'Interpreter','latex','String','$$\it{x(m)}$$',...
    'FontName','times','FontSize',15)
set(get(gca,'YLabel'),'Interpreter','latex','String','$$\it{{C}_{cl}(mol/m^3)}$$',...
    'FontName','times','FontSize',15)
set(gca,'fontsize',15);

legend('$$\it{Carbonated\ region}$$','$$\it{{C}_{cl}}$$');
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
axis([0.001 0.06 0.95*min(C) 1.1*Cs]);
legend boxoff;



