%==========================================================================
% drive_allomtests.m
%==========================================================================

clear all;
close all;

addpath('functions/');
addpath('tools/cbrewer/');
addpath('tools/xml2struct/');

xmlfile = 'allom_params.xml';

%addpath('~/local/Matlab/cbrewer/');

% Initialize the allometry library and the pointers
[f_h,f_bag,f_blmax,f_h2d, ...
    f_bsap,f_bcr,f_bfrmax,f_bdead] = allom_lib_v4;

% Some testing constants
ndbh = 2000;
maxdbh = 150.0;

% Load Parameter Values from XML
% =========================================================================
[pftcon] = read_xml_params(xmlfile);
n_pfts = numel(pftcon.tag);

% =========================================================================
% Generate a vector of diameters that starts at the smallest known diameter
% and extends to 150cm

% =========================================================================
% Initialize Output Arrays

blmaxi  = zeros(n_pfts,ndbh);
blmaxd  = zeros(n_pfts,ndbh);

bfrmax = zeros(n_pfts,ndbh);
hi     = zeros(n_pfts,ndbh);
bagi   = zeros(n_pfts,ndbh);
bagd   = zeros(n_pfts,ndbh);
dbh    = zeros(n_pfts,ndbh);
bcr    = zeros(n_pfts,ndbh);
bsapi  = zeros(n_pfts,ndbh);
bsapd  = zeros(n_pfts,ndbh);
bdead  = zeros(n_pfts,ndbh);
dbhe   = zeros(n_pfts,ndbh);

blmax_o_dbagdh = zeros(n_pfts,ndbh);


for ip=1:n_pfts
    [pftcon.dbh_min(ip),~]  = f_h2d(pftcon.h_min(ip),pftcon,ip);
    [pftcon.dbh_maxh(ip),~] = f_h2d(pftcon.h_max(ip),pftcon,ip);
    dbh(ip,:)               = linspace(pftcon.dbh_min(ip),maxdbh,ndbh);
end

% =========================================================================
% Define Tests
%
% 1) all relationships are monotonip increasing or flat
% 2) real numbers
% 3) no zeros in any non-derivatives
% 4) minimal error between integrated and explipit functions
%                                  | dbh_min -> dbh_hmax
% 5) Sanity Checks
% 6) 0.1% bag  < blmax < 10% bag   | dbh_min -> dbh_hmax
% 7) 0.1% bag  < bfrmax < 10% bag  | dbh_min -> dbh_hmax
% 8) 3% bag < bcr < 90% bag        | dbh_min -> dbh_hmax
% 9) 0.1% bag  < bsap < 60% bag    | dbh_min -> dbh_hmax
% 10) Visual tests, plot all functional relationships x-y
%     all cases are displayed in each x/y graph for:
%     d2h, h2d, d2bag, d2blmax, d2sap, d2bcr, d2bfrmax, d2bdead
%
% =========================================================================

% Calculate Results from Cases
% =========================================================================
for ip=1:n_pfts
    
    % Initialize Both Integrated and Absolute Quantities
    % ---------------------------------------------------------------------
    d0 = dbh(ip,1);
    [h0,dhdd0]           = f_h(d0,pftcon,ip);
    % Height Integrated
    hi(ip,1)          = h0;
    % AGB
    [bagi(ip,1),~]   = f_bag(d0,h0,pftcon,ip);
    % AGB
    [bagd(ip,1),dbagd0]    = f_bag(d0,h0,pftcon,ip);
    % blmax (integrated)
    [blmaxi(ip,1),~]  = f_blmax(d0,h0,pftcon,ip);
    % blmax (integrated)
    [blmaxd(ip,1),~]  = f_blmax(d0,0,pftcon,ip);
    % bfrmax
    [bfrmax(ip,1),~] = f_bfrmax(d0,blmaxi(ip,1),0,pftcon,ip);
    % bcr
    [bcr(ip,1),~]    = f_bcr(d0,bagi(ip,1),0,pftcon,ip);
    % bsap (integrated)
    [bsapi(ip,1),~]   = f_bsap(d0,hi(ip,1),blmaxi(ip,1),0,0,pftcon,ip);
    % bsap (direct)
    [bsapd(ip,1),~]   = f_bsap(d0,hi(ip,1),blmaxi(ip,1),0,0,pftcon,ip);
    % bdead
    [bdead(ip,1),~]  = f_bdead(bagi(ip,1),bcr(ip,1),blmaxi(ip,1),bsapi(ip,1),0,0,0,0);
    % Reverse (effective diameter)
    dbhe(ip,1)       = dbh(ip,1);
    
    
    blmax_o_dbagdh(ip,1)  = blmaxi(ip,1)./(dbagd0/dhdd0);
    
    % Walk through a range of diameters
    for id=2:ndbh
        dp  = dbh(ip,id-1);
        dc  = dbh(ip,id);
        dd = dc-dp;
        
        % Height
        [~,dhdd] = f_h(dp,pftcon,ip);
        hi(ip,id) = hi(ip,id-1) + dhdd*dd;
        
        % Effective diameter( hd2(h) )
        %if(ip==3)
        [dbhe(ip,id),~] = f_h2d(hi(ip,id),pftcon,ip);
        %end
        
        % AGB (integrated)
        [~,dbagdd]  = f_bag(dp,hi(ip,id-1),pftcon,ip);
        bagi(ip,id) = bagi(ip,id-1) + dbagdd*dd;
        
        % AGB (direct)
        [bagd(ip,id),dbagdd] = f_bag(dc,hi(ip,id),pftcon,ip);
        
        % blmax (integrated)
        [~,dblmaxdd]  = f_blmax(dp,hi(ip,id-1),pftcon,ip);
        blmaxi(ip,id) = blmaxi(ip,id-1) + dblmaxdd*dd;
        
        % blmax (direct)
        [blmaxd(ip,id),~] = f_blmax(dc,hi(ip,id),pftcon,ip);
        
        % bfrmax
        [~,dbfrmaxdd] = f_bfrmax(dp,blmaxi(ip,id-1),dblmaxdd,pftcon,ip);
        bfrmax(ip,id) = bfrmax(ip,id-1) + dbfrmaxdd*dd;
        
        % bcr
        [~,dbcrdd]            = f_bcr(dp,bagi(ip,id-1),dbagdd,pftcon,ip);
        bcr(ip,id) = bcr(ip,id-1) + dbcrdd*dd;
        
        % bsap (integrated)
        [~,dbsapdd]  = f_bsap(dp,hi(ip,id-1),blmaxi(ip,id-1),dblmaxdd,dhdd,pftcon,ip);
        bsapi(ip,id) = bsapi(ip,id-1) + dbsapdd*dd;
        
        % bsap (direct)
        [bsapd(ip,id),~] = f_bsap(dp,hi(ip,id),blmaxi(ip,id),0,0,pftcon,ip);
        
        if(hi(ip,id)>=pftcon.h_max(ip))
            blmax_o_dbagdh(ip,id) = NaN;
        else
            blmax_o_dbagdh(ip,id)  = blmaxi(ip,id)./(dbagdd./dhdd);
        end
        
        %display([bagi(ip,id)+bcr(ip,id)-blmaxi(ip,id)-bsapi(ip,id),bagi(ip,id),bcr(ip,id),blmaxi(ip,id),bsapi(ip,id)])
        %pause;
        
        % bdead
        [~,dbdeaddd]            = f_bdead(bagi(ip,id-1),bcr(ip,id-1), ...
                                blmaxi(ip,id-1),bsapi(ip,id-1), ...
                                dbagdd,dbcrdd,dblmaxdd,dbsapdd);
        bdead(ip,id) = bdead(ip,id-1) + dbdeaddd*dd;
        
    end
    
    
end



% Test 1: Make direct plots of each
fig=1; plot_multicase_h(dbh,hi,pftcon,fig);
fig=fig+1; plot_multicase_bag(dbh,bagi,pftcon,fig);
fig=fig+1; plot_multicase_blmax(dbh,blmaxi,pftcon,fig);
%fig=fig+1; plot_multicase_bfrmax(dbh,bfrmaxi,pftcon,fig);
%fig=fig+1; plot_multicase_bcr(state,pftcon,fig);
fig=fig+1; plot_multicase_bsap(dbh,bsapd,pftcon,fig);
fig=fig+1; plot_multicase_bsapid(bsapd,bsapi,pftcon,fig);
fig=fig+1; plot_multicase_dbhe(dbh,dbhe,pftcon,fig);
fig=fig+1; plot_multicase_blmaxdi(blmaxi,blmaxd,pftcon,fig);
fig=fig+1; plot_multicase_bsapbag(dbh,bsapi./bagi,pftcon,fig);
fig=fig+1; plot_multicase_blmax_o_dbagdh(dbh,blmax_o_dbagdh,pftcon,fig);

for ip=1:n_pfts
    fig=fig+1; plot_singlecase_cfractions(dbh(ip,:),blmaxi(ip,:), ...
                    bfrmax(ip,:),bcr(ip,:),bsapi(ip,:),bdead(ip,:), ...
                    pftcon,ip,fig)
    pause;
end
