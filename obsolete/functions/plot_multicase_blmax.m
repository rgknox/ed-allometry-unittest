function plot_multicase_blmax(dbh,blmax,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);


[nsets,ntimes] = size(blmax);

usesets = 1:nsets;

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,9,4.5],'Color','w');
subplot(1,2,1);
phan=plot(dbh(usesets,:)',blmax(usesets,:)','LineWidth',2);
xlim([0,20]);
xlabel('DBH [cm]');
ylabel('Maximum Leaf Biomass [kg]');
title(sprintf('Allometry Visual Tests\nIntercomparison of Leaf Biomas and Diameter Curves'));
grid on;
box on;
legend(traitp.tag,'Location','NorthWest');

subplot(1,2,2);
phan=loglog(dbh(usesets,:)',blmax(usesets,:)','LineWidth',2);
xlim([0,max(max(dbh(usesets,:)))]);
xlabel('DBH [cm]');
ylabel('Maximum Leaf Biomass [kg]');
title(sprintf('Allometry Visual Tests\nIntercomparison of Leaf Biomas and Diameter Curves'));
grid on;
box on;
