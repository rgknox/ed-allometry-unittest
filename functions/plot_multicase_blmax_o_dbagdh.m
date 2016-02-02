function plot_multicase_blmax_o_dbagdh(dbh,blmax_o_dhdbag,traitp,fign)

[ncases,~] = size(dbh);

cmap = cbrewer('qual','Set2',ncases);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);


[nsets,ntimes] = size(dbh);


blmax_o_dhdbag(blmax_o_dhdbag>1e7 | blmax_o_dhdbag<1e-7)=NaN;

usesets = 1:nsets;

figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,7,6],'Color','w');
phan=plot(dbh(usesets,:)',blmax_o_dhdbag(usesets,:)','LineWidth',2);
xlabel('DBH [cm]');
ylabel('bl / dAGB/dH [m]');
title(sprintf('Leaf Biomas per Relative Carbon Gain Rate per Height Rate '));
grid on;
box on;
legend(traitp.tag,'Location','NorthWest');
