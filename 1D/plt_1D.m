load dump_ana.bin;
load dump_LF.bin;
load dump_SF.bin;

figure; hold on;

plot(plt_ana(:,1),plt_ana(:,2),'k-','Linewidth',2);
plot(plt_lf(:,1),plt_lf(:,2),'k--','LineWidth',2);
plot(plt_SF(:,1),plt_SF(:,2),'k.-', 'LineWidth',1,'MarkerSize',20);
xlabel('normalized distance x/L','FontSize',16);
ylabel('water saturation','FontSize', 16);
ylim([0.1 0.85])
xlim([0.0 0.99])
box on; grid on;
lgd = legend('Analytic','FD','MFSF');
legend('boxoff');
set(lgd,'FontSize',16);

%% compute error norm

xs = plt_ana(end-1,1);
id_lf_2 = find(plt_lf(:,1) > xs);
id_lf_1 = find(plt_lf(:,1) <= xs);

rmse_lf = sqrt((sum(abs(interp1(plt_ana(1:end-2,1),plt_ana(1:end-2,2),plt_lf(id_lf_1,1), 'linear') - plt_lf(id_lf_1,2)).^2) + 
sum(abs(plt_lf(id_lf_2,2)).^2))/length(plt_lf))

rmse_SF = sqrt(sum((plt_SF(:,2)-plt_ana(:,2)).^2) / length(plt_SF))

text(0.1,0.4,['RMSE_{FD}   = ', num2str(rmse_lf)],'FontSize',18)
text(0.1,0.35,['RMSE_{MFSF} = ', num2str(rmse_SF)],'FontSize',18)