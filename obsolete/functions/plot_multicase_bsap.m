function plot_multicase_bsap(dbh,bsap,traitp,fign)

[ncases,~] = size(bsap);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(dbh',bsap','LineWidth',2);
xlabel('dbh [cm]');
ylabel('bsap [cm]');
title(sprintf('Allometry Visual Tests\nSapwood Biomas [kgC] per diameter'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');