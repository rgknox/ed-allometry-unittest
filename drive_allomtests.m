%==========================================================================
% drive_allomtests.m
%==========================================================================

clear all;
close all;

addpath('functions/');
addpath('tools/cbrewer/');

%addpath('~/local/Matlab/cbrewer/');

% Initialize the allometry library and the pointers
[f_h,f_bag,f_blmax,f_h2d, ...
    f_bsap,f_bcr,f_bfrmax,f_bdead] = allom_lib_v3;

% Some testing constants
ndbh = 2000;
maxdbh = 150.0;

% Define the parameter spaces
% =========================================================================

[traitp] = gen_param_instance;

nc = numel(traitp.wood_density);

% =========================================================================
% Generate a vector of diameters that starts at the smallest known diameter
% and extends to 150cm

% =========================================================================
% Initialize Output Arrays

blmaxi  = zeros(nc,ndbh);
blmaxd  = zeros(nc,ndbh);

bfrmax = zeros(nc,ndbh);
hi     = zeros(nc,ndbh);
bagi   = zeros(nc,ndbh);
bagd   = zeros(nc,ndbh);
dbh    = zeros(nc,ndbh);
bcr    = zeros(nc,ndbh);
bsapi   = zeros(nc,ndbh);
bsapd  = zeros(nc,ndbh);
bdead  = zeros(nc,ndbh);
dbhe   = zeros(nc,ndbh);

blmax_o_dbagdh = zeros(nc,ndbh);


for ic=1:nc
    [traitp.dbh_min(ic),~]  = f_h2d(traitp.h_min(ic),traitp,ic);
    [traitp.dbh_maxh(ic),~] = f_h2d(traitp.h_max(ic),traitp,ic);
    dbh(ic,:)               = linspace(traitp.dbh_min(ic),maxdbh,ndbh);
end

% =========================================================================
% Define Tests
%
% 1) all relationships are monotonic increasing or flat
% 2) real numbers
% 3) no zeros in any non-derivatives
% 4) minimal error between integrated and explicit functions
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
for ic=1:nc
    
    % Initialize Both Integrated and Absolute Quantities
    % ---------------------------------------------------------------------
    d0 = dbh(ic,1);
    [h0,dhdd0]           = f_h(d0,traitp,ic);
    % Height Integrated
    hi(ic,1)          = h0;
    % AGB
    [bagi(ic,1),~]   = f_bag(d0,h0,traitp,ic);
    % AGB
    [bagd(ic,1),dbagd0]    = f_bag(d0,h0,traitp,ic);
    % blmax (integrated)
    [blmaxi(ic,1),~]  = f_blmax(d0,h0,traitp,ic);
    % blmax (integrated)
    [blmaxd(ic,1),~]  = f_blmax(d0,0,traitp,ic);
    % bfrmax
    [bfrmax(ic,1),~] = f_bfrmax(d0,blmaxi(ic,1),0,traitp,ic);
    % bcr
    [bcr(ic,1),~]    = f_bcr(d0,bagi(ic,1),0,traitp,ic);
    % bsap (integrated)
    [bsapi(ic,1),~]   = f_bsap(d0,hi(ic,1),blmaxi(ic,1),0,0,traitp,ic);
    % bsap (direct)
    [bsapd(ic,1),~]   = f_bsap(d0,hi(ic,1),blmaxi(ic,1),0,0,traitp,ic);
    % bdead
    [bdead(ic,1),~]  = f_bdead(bagi(ic,1),bcr(ic,1),blmaxi(ic,1),bsapi(ic,1),0,0,0,0);
    % Reverse (effective diameter)
    dbhe(ic,1)       = dbh(ic,1);
    
    
    blmax_o_dbagdh(ic,1)  = blmaxi(ic,1)./(dbagd0/dhdd0);
    
    % Walk through a range of diameters
    for id=2:ndbh
        dp  = dbh(ic,id-1);
        dc  = dbh(ic,id);
        dd = dc-dp;
        
        % Height
        [~,dhdd] = f_h(dp,traitp,ic);
        hi(ic,id) = hi(ic,id-1) + dhdd*dd;
        
        % Effective diameter( hd2(h) )
        %if(ic==3)
        [dbhe(ic,id),~] = f_h2d(hi(ic,id),traitp,ic);
        %end
        
        % AGB (integrated)
        [~,dbagdd]  = f_bag(dp,hi(ic,id-1),traitp,ic);
        bagi(ic,id) = bagi(ic,id-1) + dbagdd*dd;
        
        % AGB (direct)
        [bagd(ic,id),dbagdd] = f_bag(dc,hi(ic,id),traitp,ic);
        
        % blmax (integrated)
        [~,dblmaxdd]  = f_blmax(dp,hi(ic,id-1),traitp,ic);
        blmaxi(ic,id) = blmaxi(ic,id-1) + dblmaxdd*dd;
        
        % blmax (direct)
        [blmaxd(ic,id),~] = f_blmax(dc,hi(ic,id),traitp,ic);
        
        % bfrmax
        [~,dbfrmaxdd] = f_bfrmax(dp,blmaxi(ic,id-1),dblmaxdd,traitp,ic);
        bfrmax(ic,id) = bfrmax(ic,id-1) + dbfrmaxdd*dd;
        
        % bcr
        [~,dbcrdd]            = f_bcr(dp,bagi(ic,id-1),dbagdd,traitp,ic);
        bcr(ic,id) = bcr(ic,id-1) + dbcrdd*dd;
        
        % bsap (integrated)
        [~,dbsapdd]  = f_bsap(dp,hi(ic,id-1),blmaxi(ic,id-1),dblmaxdd,dhdd,traitp,ic);
        bsapi(ic,id) = bsapi(ic,id-1) + dbsapdd*dd;
        
        % bsap (direct)
        [bsapd(ic,id),~] = f_bsap(dp,hi(ic,id),blmaxi(ic,id),0,0,traitp,ic);
        
        if(hi(ic,id)>=traitp.h_max(ic))
            blmax_o_dbagdh(ic,id) = NaN;
        else
            blmax_o_dbagdh(ic,id)  = blmaxi(ic,id)./(dbagdd./dhdd);
        end
        
        %display([bagi(ic,id)+bcr(ic,id)-blmaxi(ic,id)-bsapi(ic,id),bagi(ic,id),bcr(ic,id),blmaxi(ic,id),bsapi(ic,id)])
        %pause;
        
        % bdead
        [~,dbdeaddd]            = f_bdead(bagi(ic,id-1),bcr(ic,id-1), ...
                                blmaxi(ic,id-1),bsapi(ic,id-1), ...
                                dbagdd,dbcrdd,dblmaxdd,dbsapdd);
        bdead(ic,id) = bdead(ic,id-1) + dbdeaddd*dd;
        
    end
    
    
end



% Test 1: Make direct plots of each
fig=1; plot_multicase_h(dbh,hi,traitp,fig);
fig=fig+1; plot_multicase_bag(dbh,bagi,traitp,fig);
fig=fig+1; plot_multicase_blmax(dbh,blmaxi,traitp,fig);
%fig=fig+1; plot_multicase_bfrmax(dbh,bfrmaxi,traitp,fig);
%fig=fig+1; plot_multicase_bcr(state,traitp,fig);
fig=fig+1; plot_multicase_bsap(dbh,bsapd,traitp,fig);
fig=fig+1; plot_multicase_bsapid(bsapd,bsapi,traitp,fig);
fig=fig+1; plot_multicase_dbhe(dbh,dbhe,traitp,fig);
fig=fig+1; plot_multicase_blmaxdi(blmaxi,blmaxd,traitp,fig);
fig=fig+1; plot_multicase_bsapbag(dbh,bsapi./bagi,traitp,fig);
fig=fig+1; plot_multicase_blmax_o_dbagdh(dbh,blmax_o_dbagdh,traitp,fig);

for ic=1:nc
    fig=fig+1; plot_singlecase_cfractions(dbh(ic,:),blmaxi(ic,:), ...
                    bfrmax(ic,:),bcr(ic,:),bsapi(ic,:),bdead(ic,:), ...
                    traitp,ic,fig)
    pause;
end
