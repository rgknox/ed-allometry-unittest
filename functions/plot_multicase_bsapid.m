function plot_multicase_bsapid(bsapd,bsapi,traitp,fign)

[ncases,~] = size(bsapd);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,6,4.5],'Color','w');
phan=plot(bsapd,bsapi,'LineWidth',2);
xlabel('bsap (direct method) [kgC]');
ylabel('bsap (integrated method) [kgC');
title(sprintf('Sapwood Biomas [kgC] Integrator Check'));
grid on;
box on;

legend(traitp.tag,'Location','SouthEast');