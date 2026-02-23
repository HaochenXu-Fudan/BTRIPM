function y = plot_accuracy(funname,accuracy,FI)

oldfolder = cd;
epsfolder = strcat(oldfolder,'\eps');
cd(epsfolder);

figure(3)
plot(1:length(accuracy.bvfim), accuracy.bvfim, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
hold on
plot(1:length(accuracy.bvfimzj), accuracy.bvfimzj, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
hold on
plot(1:length(accuracy.Tregion),accuracy.Tregion, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
hold on
plot(1:length(accuracy.Tregionzj), accuracy.Tregionzj, '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',3);
legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
legend('BVFIM1','BVFIM2','BTRIPM1','BTRIPM2','Location','best'); 
set(gca,'FontSize',12); 
ylabel('Accuracy','FontSize',15,'LineWidth',2);
xlabel('Iteration','FontSize',15,'LineWidth',2);
epsname = strcat(funname,'accuracy.eps');
print(figure(1), '-depsc',epsname);
if FI.plot_yes
    figure(4)
    plot(1:length(FI.bvfim), FI.bvfim, '--', 'Color',[0.99 0.01 0.01], 'LineWidth',2);
    hold on
    plot(1:length(FI.bvfimzj), FI.bvfimzj, ':', 'Color',[0.4 0.2 0.7], 'LineWidth',3);
    hold on
    plot(1:length(FI.Tregion),FI.Tregion, '-', 'Color',[0.2 0.8 0.4], 'LineWidth',2);
    hold on
    plot(1:length(FI.Tregionzj), FI.Tregionzj, '-.', 'Color',[0.8 0.1 0.3], 'LineWidth',3);
    legend('FontName','Times New Roman','FontSize',15,'LineWidth',1);
    legend('BVFIM1','BVFIM2','BTRIPM1','BTRIPM2','Location','best'); 
    set(gca,'FontSize',12); 
    ylabel('F1 score','FontSize',15,'LineWidth',2);
    xlabel('Iteration','FontSize',15,'LineWidth',2);
    epsname = strcat(funname,'FI.eps');
    print(figure(2), '-depsc',epsname);
end
cd(oldfolder);
end