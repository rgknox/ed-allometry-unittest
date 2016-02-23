function plot_singlecase_cfractions(dbh,blmax,bfrmax,bcr,bsap,bdead,traitp,ic,fign)


[ndbh] = numel(blmax);


cmap = cbrewer('qual','Set2',4);
set(0,'defaultAxesColorOrder',cmap);
set(0,'defaultAxesFontSize',13);
set(0,'defaultAxesLineWidth',1);


btot = blmax+bfrmax+bdead+bsap;

afrac = zeros(ndbh,4);
atot  = zeros(ndbh,4);

afrac(:,4) = bdead./btot;
afrac(:,3) = bsap./btot;
afrac(:,1) = blmax./btot;
afrac(:,2) = bfrmax./btot;
atot(:,4) = bdead;
atot(:,3) = bsap;
atot(:,1) = blmax;
atot(:,2) = bfrmax;
atags = {'blmax','bfrmax','bsap','bdead'};


figure(fign);
clf;
set(fign,'Units','inches','Position',[fign,fign,7.5,4],'Color','w');
subplot(1,2,1);
ahan=area(dbh,afrac);
for ihan=1:length(ahan)
    set(ahan(ihan),'FaceColor',cmap(ihan,:));
end
set(gca,'Xscale','log');
xlabel('dbh [cm] (logscale)');
ylabel('[-]');
ylim([0,1]);
title(sprintf('Fraction of Total Carbon Allocated\nCase: %s',traitp.tag{ic}));
grid on;
box on;
legend(atags,'Location','SouthEast');

subplot(1,2,2);
syhan = semilogy(dbh,atot,'LineWidth',2);
for ihan=1:length(syhan)
    set(syhan(ihan),'Color',cmap(ihan,:));
    if(strcmp(atags{ihan},'bfrmax'))
        display('brow');
        set(syhan(ihan),'LineStyle','--');
    end
end
title(sprintf('Total Carbon Allocated\nCase: %s',traitp.tag{ic}));
grid on;
box on;
xlabel('dbh [cm]');
ylabel('[kgC]');
