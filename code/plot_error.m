function y = plot_error(error_bvfim,error_rhg,error_cg,error_Tregion,funname,x0str,y0str)
oldfolder = cd;
epsfolder = strcat(oldfolder,'/eps');
cd(epsfolder);

flag = 1;
% if nargin==5 
%     flag = 0;
% elseif nargin==6
%     flag = 1;
% else
%    error('plot_error: wrong input number');
% end
figure(1)
semilogy(1:length(error_bvfim.F), error_bvfim.F, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
semilogy(1:length(error_rhg.F), error_rhg.F, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
% hold on
% semilogy(1:length(error_trhg.F), error_trhg.F, '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
hold on
semilogy(1:length(error_cg.F), error_cg.F, '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
semilogy(1:length(error_Tregion.F), error_Tregion.F, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM','RHG','CG','BTRIPM','Location','South');
set(gca,'FontSize',12); 
ylabel('$|F-F^*|$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iterations','FontSize',15,'LineWidth',2);
if flag==0
    epsname = strcat(funname,'error_F.eps');
elseif flag == 1
    epsname = strcat(x0str,y0str,funname,'error_F.eps');
end
print(figure(1), '-depsc',epsname);

figure(2)
semilogy(1:length(error_bvfim.f), error_bvfim.f, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
semilogy(1:length(error_rhg.f), error_rhg.f, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
% hold on
% semilogy(1:length(error_trhg.f), error_trhg.f, '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
hold on
semilogy(1:length(error_cg.f), error_cg.f, '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
semilogy(1:length(error_Tregion.f), error_Tregion.f, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM','RHG','CG','BTRIPM','Location','South');
ylabel('$|f-f^*|$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iterations','FontSize',15,'LineWidth',2);
if flag==0
    epsname = strcat(funname,'error_fval.eps');
elseif flag == 1
    epsname = strcat(x0str,y0str,funname,'error_fval.eps');
end
print(figure(2), '-depsc',epsname);


figure(3)
semilogy(1:length(error_bvfim.x), error_bvfim.x, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
semilogy(1:length(error_rhg.x), error_rhg.x, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
% hold on
% semilogy(1:length(error_trhg.x), error_trhg.x, '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
hold on
semilogy(1:length(error_cg.x), error_cg.x, '-', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
semilogy(1:length(error_Tregion.x), error_Tregion.x, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM','RHG','CG','BTRIPM','Location','South');
ylabel('$|x-x^*|$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iteration','FontSize',15,'LineWidth',2);
if flag==0
    epsname = strcat(funname,'error_x.eps');
elseif flag == 1
    epsname = strcat(x0str,y0str,funname,'error_x.eps');
end
print(figure(1), '-depsc',epsname);

figure(4)
semilogy(1:length(error_bvfim.y), error_bvfim.y, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
semilogy(1:length(error_rhg.y), error_rhg.y, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
% hold on
% semilogy(1:length(error_trhg.y), error_trhg.y, '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
hold on
semilogy(1:length(error_cg.y), error_cg.y, '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
semilogy(1:length(error_Tregion.y), error_Tregion.y, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM','RHG','CG','BTRIPM','Location','South');
ylabel('$|y-y^*|$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iteration','FontSize',15,'LineWidth',2);
if flag==0
    epsname = strcat(funname,'error_y.eps');
elseif flag == 1
    epsname = strcat(x0str,y0str,funname,'error_y.eps');
end
print(figure(4), '-depsc',epsname);
cd(oldfolder);
end