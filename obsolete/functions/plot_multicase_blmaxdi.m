function plot_multicase_blmaxdi(blmaxi,blmaxd,traitp,fign)

[ncases,~] = size(blmaxi);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(blmaxi',blmaxd','LineWidth',2);
xlabel('blmax (direct) [kgC]');
ylabel('blmax (integrated) [kgC]');
title(sprintf('Allometry Visual Tests\nIntercomparison of Integrated and Direct Leaf Biomass'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');