figure(1);
x = 1:0.01:4; y = 0:0.01:1;
[X,Y] = meshgrid(x,y); 
mu = 3;
sig = 0.3;
Z = log(Y)-(-(X-mu).^2/(2*sig^2));

v = [0,0];
[C,h] = contour(X,Y,Z,v,'ShowText','off','Fill','off','color','k');

figure(2);
x_plot_1 = C(1,2:end); y_plot_1 = C(2,2:end); 
hold on; fill_1 = fill(x_plot_1,y_plot_1,[0.9206 0.9216 0.9204]);
axis([1 4 0 1]);

% Plot the subset samples
[C,h] = contour(X,Y,Z,v,'ShowText','off','Fill','off','color','k');
legend('$${g(x,p)\leq0}$$');
% set(get(get(fill_1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(gca,'XLabel'),'Interpreter','latex','String','$${x}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$${p}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
legend boxoff;
axis([1 4 0 1]);
box on;
% box on;