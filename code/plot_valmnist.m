function y = plot_valmnist(funname,val_bvfim,val_Tregion,val_bvfimzj,val_Tregionzj,F0,f0)
oldfolder = cd;
epsfolder = strcat(oldfolder,'/eps');
cd(epsfolder);

figure(1)
plot(0:length(val_bvfim.F), [F0;val_bvfim.F], '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
plot(0:length(val_bvfimzj.F),[F0;val_bvfimzj.F], ':', 'Color',[0.8 0.1 0.2], 'LineWidth',3);
hold on
plot(0:length(val_Tregion.F),[F0;val_Tregion.F], '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
hold on
plot(0:length(val_Tregionzj.F),[F0;val_Tregionzj.F], '-.', 'Color',[0.8 0.4 0.2], 'LineWidth',3);
 legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM1','BVFIM2','BTRIPM1','BTRIPM2','Location','best'); 
ylabel('$F(\mathbf{x},\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iteration','FontSize',15,'LineWidth',2);
epsname = strcat(funname,'val_F.eps');
print(figure(1), '-depsc',epsname);

figure(2)
plot(0:length(val_bvfim.f), [f0;val_bvfim.f], '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
plot(0:length(val_bvfimzj.f),[F0;val_bvfimzj.f], '', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
hold on
plot(0:length(val_Tregion.f), [f0;val_Tregion.f], '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
hold on
plot(0:length(val_Tregionzj.f),[F0;val_Tregionzj.f], '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',3);
 legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM1','BVFIM2','BTRIPM1','BTRIPM2','Location','best'); 
ylabel('$f(\mathbf{x},\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
xlabel('Iterations','FontSize',15,'LineWidth',2);
epsname = strcat(funname,'fval.eps');
print(figure(2), '-depsc',epsname);

cd(oldfolder);
end