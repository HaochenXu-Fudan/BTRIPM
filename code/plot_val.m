function y = plot_val(funname,val_bvfim,val_Tregion,val_rhg,val_trhg,val_cg,F0,f0)
oldfolder = cd;
epsfolder = strcat(oldfolder,'/eps');
cd(epsfolder);

figure(1)
plot(0:length(val_bvfim.F), [F0;val_bvfim.F], '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
plot(0:length(val_Tregion.F),[F0;val_Tregion.F], '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
hold on
plot(0:length(val_cg.F),[F0;val_cg.F], '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
plot(0:length(val_rhg.F),[F0;val_rhg.F], ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
hold on
plot(0:length(val_trhg.F),[F0;val_trhg.F], '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
legend('BVFIM','BTRIPM','CG','RHG','TRHG','Location','best'); 
% end
set(gca,'FontSize',12); 
if strcmp(funname,'wifi')
    ylabel('$F(x,\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
else
    ylabel('$F(\mathbf{x},\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
end
xlabel('Iteration','FontSize',15,'LineWidth',2);
epsname = strcat(funname,'val_F.eps');
print(figure(1), '-depsc',epsname);

figure(2)
plot(0:length(val_bvfim.f), [f0;val_bvfim.f], '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
plot(0:length(val_Tregion.f), [f0;val_Tregion.f], '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
hold on
plot(0:length(val_cg.f),[F0;val_cg.f], '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',2);
hold on
plot(0:length(val_rhg.f),[F0;val_rhg.f], ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
hold on
plot(0:length(val_trhg.f),[F0;val_trhg.f], '-', 'Color',[0.1 0.6 0.6], 'LineWidth',2);
set(gca,'FontSize',12); 
legend('BVFIM','BTRIPM','CG','RHG','TRHG','Location','best'); 
if strcmp(funname,'wifi')
    ylabel('$f(x,\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
else
    ylabel('$f(\mathbf{x},\mathbf{y})$','interpreter','latex','FontSize',15,'LineWidth',2);
end
xlabel('Iteration','FontSize',15,'LineWidth',2);
epsname = strcat(funname,'fval.eps');
print(figure(2), '-depsc',epsname);

cd(oldfolder);
end