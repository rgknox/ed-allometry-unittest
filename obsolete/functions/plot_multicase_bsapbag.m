function plot_multicase_bsapbag(dbh,bsapbag,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(dbh',bsapbag','LineWidth',2);
xlabel('DBH [cm]');
ylabel('Fraction of AGB in Sapwood');
title(sprintf('Allometry Visual Tests\nIntercomparison of Sapwood Fraction'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');