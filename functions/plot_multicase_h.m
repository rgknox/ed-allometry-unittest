function plot_multicase_h(dbh,h,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(dbh',h','LineWidth',2);
xlabel('DBH [cm]');
ylabel('Height [m]');
title(sprintf('Allometry Visual Tests\nIntercomparison of Diameter Height Curves'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');













end