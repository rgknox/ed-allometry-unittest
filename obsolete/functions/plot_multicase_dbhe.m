function plot_multicase_dbhe(dbh,dbhe,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(dbh',dbhe','LineWidth',2);
xlabel('DBH [cm]');
ylabel('Effective DBH [cm]');
title(sprintf('Allometry Visual Tests\nIntercomparison of Effective Diameter and Diameter Curves'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');