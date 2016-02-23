function plot_multicase_bag(dbh,bag,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(dbh',bag'./1000,'LineWidth',2);
xlabel('DBH [cm]');
ylabel('Aboveground Biomass [Mg]');
title(sprintf('Allometry Visual Tests\nIntercomparison of AGB-Diameter Curves'));
grid on;
box on;

legend(traitp.tag,'Location','NorthWest');