figure(1);
x = 0:0.01:4; y = 0:0.01:1;
[X,Y] = meshgrid(x,y); 
X_one = reshape(X,size(X,1)*size(X,2),1);
Y_one = reshape(Y,size(Y,1)*size(Y,2),1);
X_S = [Y_one X_one];

mu = 3;
sig = 0.3;
Z = log(Y)-(-(X-mu).^2/(2*sig^2));
[Z_one, mse_kr] = predictor(X_S, dmodel);
Z_kr = reshape(Z_one,size(X,1),size(X,2),1);

v = [15.3,15.3];
[C,h] = contour(X,Y,Z,v,'ShowText','off','Fill','off','color','k'); 
[K,h] = contour(X,Y,Z_kr,v,'--','ShowText','off','Fill','off','color','r');

figure(2);
plot_1 = contour(X,Y,Z,v,'ShowText','off','Fill','off','color','k'); hold on;
x_plot_1 = [C(1,2:end) 10 10]; y_plot_1 = [C(2,2:end) 10 -10]; 
hold on; 
% plot_2 = contour(X,Y,Z,v,'ShowText','off','Fill','off','color','k');
fill_1 = fill(x_plot_1,y_plot_1,[0.9206 0.9216 0.9204]);
axis([0 2 0 1]);

% Plot the estimated limit state


plot_2 = contour(X,Y,Z_kr,v,'--','ShowText','off','Fill','off','color','r'); hold on;
scatter(x_tr(:,2),x_tr(:,1),20,'filled','MarkerFaceColor',[118 128 105]/255); hold on;

legend('$${h(x,p)=t_1}$$','$${\Omega_{1}}$$','$${\hat{h}(x,p)=t_1}$$','$${x_{tr}}$$');
% set(get(get(fill_1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(gca,'XLabel'),'Interpreter','latex','String','$${x}$$',...
    'FontName','times','FontSize',15);
set(get(gca,'YLabel'),'Interpreter','latex','String','$${p}$$',...
    'FontName','times','FontSize',15);
set(gca,'fontsize',15);
set(get(gca,'legend'),'Interpreter','latex',...
    'FontName','times','FontSize',15);
% legend boxoff;
axis([0 2 0 1]);
box on;
% box on;
