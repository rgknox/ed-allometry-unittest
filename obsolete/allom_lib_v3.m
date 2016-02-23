
%==========================================================================
%
% allom_lib_v3.m
% A library of functions that calculate plant allometry and their
% derivatives.  All relationships are related to diameter, all
% derivatives are assumed to be with respect to change in diameter with
% units of [cm]. In some cases multiple areguments will be provided, yet
% those arguments in those cases are trivially determined from diameter.
%
% Each function presented in this library is written to also return the
% derivative with respect to diameter using a logical switch "dswitch"
% (derivative-switch). With one exception, the h2d function returns the
% change in diameter with respect to height.
%
% The name convention of the functions follows the form d2...  Which
% indicates "diameter to ...".  Allometries for the following variables are
% calculated:
% h:  height [m]
% bag:  biomass above ground [kgC]  (aka AGB)
% blmax:  biomass in leaves when leaves are "on allometry"
%         this also is the "pre-trimmed" value, which is the maximum
%         or potential leaf mass a tree may have [kgC]
% bcr: biomass in coarse roots [kgC] (belowground sap+dead wood, no fines)
% bfrmax: biomass in fine roots when "on allometry" [kgC]
% bsap: biomass in sapwood (above and below) [kgC]
% bdead: biomass (above and below) in the structural pool [kgC]
%
% "on allometry" assumes the plant has leaves flushed, and has had
% sufficient carbon to meet maintenance turnover.
%
% The following model traits are used:
% h_mode, string, height allometry function type
% blmax_mode, string, maximum leaf allometry function type
% bfrmax_mode, string, maximum root allometry function type
% bag_mode, string, AGB allometry function type
% bcr_mode, string, coarse root allometry function type
% bsap_mode, string, sapwood allometry function type
% wood_density, real, mean stem wood specific gravity (heart,sap,bark)
% latosa_int, real, leaf area to sap area ratio, intercept [m2/cm2]
% latosa_slp, real, leaf area to sap area ratio, slope on diameter
% [m2/cm2/cm]
% c2b, real, carbon to biomass ratio (~2.0)
% d2h1, real, parameter 1 for d2h allometry (intercept)
% d2h2, real, parameter 2 for d2h allometry (slope)
% d2h3, real, parameter 3 for d2h allometry (optional)
% eclim, real, climatological influence parameter for d2h allometry
% bl_min, real, leaf biomass of smallest tracked plant [kgC]
% froot_leaf, real, fine root biomass per leaf biomass ratio [kgC/kgC]
% ag_biomass, real, the fraction of stem above ground [-]
% dbh_adult, real, the lower dbh bound where adult allometry starts (10m)
% h_max, real, maximum height of a functional type/group
% d2bl1, real, parameter 1 for d2bl allometry (intercept)
% d2bl2, real, parameter 2 for d2bl allometry (slope)
% d2bl3, real, parameter 3 for d2bl allometry (optional)
% h_min, real, the height associated with newly recruited plant [m]
% dbh_min, real, the dbh associated with a newly recruited plant [cm]
% dbh_max, real, the diameter associated with maximum height [cm]
%               diagnosed from hmax using non-asymptotic functions
%
%==========================================================================

function [p_h,p_bag,p_blmax,p_h2d, ...
    p_bsap,p_bcr,p_bfrmax,p_bdead] = allom_lib_v3

p_h                   = @f_h;
p_bag                 = @f_bag;
p_blmax               = @f_blmax;
p_h2d                 = @f_h2d;
p_bcr                 = @f_bcr;
p_bsap                = @f_bsap;
p_bfrmax              = @f_bfrmax;
p_bdead               = @f_bdead;
end

% =========================================================================
% Parameter Checks and Defaults (subroutine)
% =========================================================================





% =========================================================================
% Generic height to diameter interface
% =========================================================================

function [d,dddh] = f_h2d(h,traitp,ipf)
switch traitp.h_mode{ipf}
    case {'chave14'}
        [d,dddh] = h2d_chave2014(h,traitp,ipf);
    case {'poorter06'}
        [d,dddh] = h2d_poorter2006(h,traitp,ipf);
    case {'2par_pwr'}
        [d,dddh] = h2d_2pwr(h,traitp,ipf);
    case {'obrien'}
        [d,dddh] = h2d_obrien(h,traitp,ipf);
    otherwise
        display('Unknown D-2-H Allometry');
        pause;
        return;
end
end

% =========================================================================
% Generic height interface
% =========================================================================

function [h,dhdd] = f_h(d,traitp,ipf)
switch traitp.h_mode{ipf}
    case {'chave14'}
        [h,dhdd] = d2h_chave2014(d,traitp,ipf);
    case {'poorter06'}
        [h,dhdd] = d2h_poorter2006(d,traitp,ipf);
    case {'2par_pwr'}
        [h,dhdd] = d2h_2pwr(d,traitp,ipf);
    case {'obrien'}
        [h,dhdd] = d2h_obrien(d,traitp,ipf);
    otherwise
        display('Unknown H Allometry');
        pause;
        return;
end
end

% =========================================================================
% Generic AGB interface
% =========================================================================

function [bag,dbagdd] = f_bag(d,h,traitp,ipf)
switch traitp.bag_mode{ipf}
    case {'chave14'}
        [bag,dbagdd] = dh2bag_chave2014(d,h,traitp,ipf);
    case {'2par_pwr'}
        % Switch for woodland dbh->drc
        [bag,dbagdd] = d2bag_2pwr(d,traitp,ipf);
    case {'salda'}
        [bag,dbagdd] = dh2bag_salda(d,h,traitp,ipf);
    otherwise
        display('Unknown BAG Allometry');
        pause;
        return;
end
end

% =========================================================================
% Generic diameter to maximum leaf biomass interface
% =========================================================================

function [blmax,dblmaxdd] = f_blmax(d,h,traitp,ipf)

switch traitp.blmax_mode{ipf}
    case{'salda'}
        [blmax,dblmaxdd] = d2blmax_salda(d,traitp,ipf);
    case{'2par_pwr'}
        [blmax,dblmaxdd] = d2blmax_2pwr(d,traitp,ipf);
    case{'2par_pwr_asym'}
        [blmax,dblmaxdd] = dh2blmax_2pwr_spline_asym(d,h,traitp,ipf);
    case{'2par_pwr_hcap'}
        [blmax,dblmaxdd] = d2blmax_2pwr_hcap(d,traitp,ipf);
    otherwise
        display(sprintf('Uknown blmax Allometry: %s',traitp.blmax_mod{ipf}));
        pause;
        return;
end
end

% =========================================================================
% Generic sapwood biomass interface
% =========================================================================

function [bsap,dbsapdd] = f_bsap(d,h,blmax,dblmaxdd,dhdd,traitp,ipf)

switch traitp.bsap_mode{ipf}
    % ---------------------------------------------------------------------
    % Currently both sapwood area proportionality methods use the same
    % machinery.  The only differences are related to the parameter
    % checking at the beginning.  For constant proportionality, the slope
    % of the la:sa to diameter line is zero.
    % ---------------------------------------------------------------------
    case{'constant','dlinear'}
        [bsap,dbsapdd] = bsap_dlinear(d,h,blmax,dblmaxdd,dhdd,traitp,ipf);
    otherwise
        display('Uknown BSAP Allometry');
        pause;
        return;
end
end

% =========================================================================
% Generic coarse root biomass interface
% =========================================================================

function [bcr,dbcrdd] = f_bcr(d,bag,dbagdd,traitp,ipf)

switch traitp.bcr_mode{ipf}
    case{'constant'}
        [bcr,dbcrdd] = bcr_const(d,bag,dbagdd,traitp,ipf);
    otherwise
        display('Uknown D-2-BCR Allometry');
        pause;
        return;
end
end

% =========================================================================
% Generic maximum fine root biomass interface
% =========================================================================

function [bfrmax,dbfrmaxdd] = f_bfrmax(d,blmax,dblmaxdd,traitp,ipf)

switch traitp.bfrmax_mode{ipf}
    case{'constant'}
        [bfrmax,dbfrmaxdd] = bfrmax_const(d,blmax,dblmaxdd,traitp,ipf);
    otherwise
        display('Uknown D-2-BFR Allometry');
        pause;
        return;
end
end

% =========================================================================
% Dead biomass interface
% =========================================================================

function [bdead,dbdeaddd] = f_bdead(bag,bcr,blmax,bsap, ...
    dbagdd,dbcrdd,dblmaxdd,dbsapdd) %#ok<INUSL>

% bdead is diagnosed as the mass balance from all other pools
% and therefore, no options are necessary
% We are ignoring blmax right now, because it is insignificant in large
% trees and may cause negatives in treelets and saplings

%bdead = bag+bcr-blmax-bsap;
bdead = bag+bcr-bsap;
dbdeaddd = dbagdd+dbcrdd-dbsapdd;

end

% =========================================================================
% Specific bfrmax relationships
% =========================================================================

function [bfrmax,dbfrmaxdd] = bfrmax_const(d,blmax,dblmaxdd,traitp,ipf) %#ok<INUSL>

froot_leaf = traitp.froot_leaf(ipf); % fine root biomass [kgC] per leaf biomass [kgC]

bfrmax = blmax.*froot_leaf;
% dbfr/dd = dbfrmax/dblmax * dblmax/dd
dbfrmaxdd = froot_leaf.*dblmaxdd;

end


% =========================================================================
% Specific bcr relationships
% =========================================================================

function [bcr,dbcrdd] = bcr_const(d,bag,dbagdd,traitp,ipf) %#ok<INUSL>

ag_biomass = traitp.ag_biomass(ipf); % Above-ground biomass fraction

% btot = bag + bcr;
% bag = btot*ag_biomass;
% bag/ag_biomass = bag + bcr;
% bcr = bag*(1/ag_biomass-1)

bcr = bag*(1/ag_biomass-1);
% Derivative
% dbcr/dd = dbcr/dbag * dbag/dd
dbcrdd = (1/ag_biomass-1).*dbagdd;

end


% =========================================================================
% Specific d2bsap relationships
% =========================================================================


function [bsap,dbsapdd] = bsap_dlinear(dbh,h,blmax,dblmaxdd,dhdd,traitp,ipf)

% -------------------------------------------------------------------------
% Calculate sapwood biomass based on leaf area to sapwood area
% proportionality.  In this function, the leaftosapwood area is a function
% of plant size, see Calvo-Alvarado and Bradley Christoferson
% In this case: parameter latosa (from constant proportionality)
%   is the intercept of the diameter function.
%
% For very small plants, the fraction can get very large, so cap the amount
% of sapwood at X% of agb-bleaf;
% -------------------------------------------------------------------------

max_agbfrac = 0.75;

gtokg  = 1000.0; % gram to kg conversion
cm2tom2 = 10000.0; % cm^2 to m^2 conversion
mg2kg   = 1000.0;  % megagrams to kilograms conversion [kg/Mg]
latosa_int = traitp.latosa_int(ipf); % leaf area to sap area ratio [m2/cm2]
latosa_slp = traitp.latosa_slp(ipf); % slope of diameter to latosa line [m2/cm2/cm]
sla          = traitp.slatop(ipf); % Specific leaf area [m2/gC]
wood_density = traitp.wood_density(ipf); % Wood density (specific grav) [Mg/m3]
c2b          = traitp.c2b(ipf);  % carbon to biomass ratio

% -------------------------------------------------------------------------
% Calculate sapwood biomass per linear height and kgC of leaf [m-1]
% Units: (1/latosa)* slatop*    gtokg    *   cm2tom2     / c2b   * mg2kg  * dens
%        [cm2/m2]  *[m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3]
%                   ->[cm2/gC]
%                          ->[cm2/kgC]
%                                         ->[m2/kgC]
%                                                        ->[m2/kg]
%                                                                 ->[m2/Mg]
%                                                                         ->[/m]
% -------------------------------------------------------------------------

latosa = latosa_int + dbh*latosa_slp;

hbl2bsap = sla*gtokg*wood_density*mg2kg/(latosa*c2b*cm2tom2);
[bag,dbagdd] = f_bag(dbh,h,traitp,ipf);

bsap = min(max_agbfrac*bag,hbl2bsap * h * blmax);

% Derivative
% dbldmaxdd is deriv of blmax wrt dbh (use directives to check oop)
% dhdd is deriv of height wrt dbh (use directives to check oop)

if(bsap<max_agbfrac*bag)
%    dbsapdd = min(max_agbfrac*dbagdd,hbl2bsap*(h*dblmaxdd + blmax*dhdd));
    dbsapdd = hbl2bsap*(h*dblmaxdd + blmax*dhdd);
else
    dbsapdd = max_agbfrac*dbagdd;
end

end

% =========================================================================
% Specific d2blmax relationships
% =========================================================================

function [blmax,dblmaxdd] = d2blmax_salda(dbh,traitp,ipf)

d2bl1    = traitp.d2bl1(ipf);   %0.0419;
d2bl2    = traitp.d2bl2(ipf);   %1.56;
d2bl3    = traitp.d2bl3(ipf);   %0.55;
rho      = traitp.wood_density(ipf);
dbh_maxh = traitp.dbh_maxh(ipf);
c2b       = traitp.c2b(ipf);
dbh_adult = traitp.dbh_adult(ipf);
dbh_sap   = traitp.dbh_sap(ipf);
d2bl1_sap = traitp.d2bl1_sap(ipf); %0.0201;
d2bl2_sap = traitp.d2bl2_sap(ipf); %3.1791;

% =========================================================================
% From King et al. 1990 at BCI for saplings
%
% log(bl) = a2 + b2*log(h)
% bl = exp(a2) * h^b2
%
% and:
%
% log(d) = a1 + b1*log(h)
% d = exp(a1) * h^b1
% h = (1/exp(a1)) * d^(1/b1)
%
% bl = exp(a2) * ((1/exp(a1)) * d^(1/b1))^b2
% bl = exp(a2) * (1/exp(a1))^b2 * d^(b2/b1)
% bl = d2bl1 * d ^ d2bl2
% where: d2bl1 = exp(a2) * (1/exp(a1))^b2
%        d2bl2 = (b2/b1)
% For T. tuberculata (canopy tree):
% a1 = -0.0704, b1 = 0.67
% a2 = -4.056,  b2 = 2.13
% ========================================================================

% Saldarriaga has a height cap on leaf biomass
if(dbh<dbh_sap)
    
    % Follow King's log-log linear regression
    blmax    = d2bl1_sap * dbh.^d2bl2_sap ./c2b;
    dblmaxdd = d2bl1_sap * d2bl2_sap * dbh.^(d2bl2_sap-1.0)./c2b;
    
elseif(dbh>=dbh_sap && dbh<dbh_adult)
    
    % Use cubic spline interolation from dbh_sap to dbh_adult
    
    blsap = d2bl1_sap * dbh_sap.^d2bl2_sap./c2b;
    blad  = d2bl1 * dbh_adult.^d2bl2 * rho.^d2bl3;
    dblsapdd = d2bl1_sap * d2bl2_sap * dbh_sap.^(d2bl2_sap-1.0)./c2b;
    dbladdd  = d2bl1 * d2bl2 * dbh_adult.^(d2bl2-1.0) * rho.^d2bl3;
    
    [blmax,dblmaxdd] = cspline(dbh_sap,dbh_adult,blsap,blad,dblsapdd,dbladdd,dbh);
    
elseif(dbh>=dbh_adult && dbh<dbh_maxh)
    blmax = d2bl1 * dbh.^d2bl2 * rho.^d2bl3;
    dblmaxdd = d2bl1*d2bl2 * dbh.^(d2bl2-1.0) * rho.^d2bl3;
else
    blmax    = d2bl1 * dbh_maxh.^d2bl2 * rho.^d2bl3;
    dblmaxdd = 0.0;
end

end


function [blmax,dblmaxdd] = d2blmax_2pwr(dbh,traitp,ipf)

d2bl1     = traitp.d2bl1(ipf);
d2bl2     = traitp.d2bl2(ipf);
c2b       = traitp.c2b(ipf);
bl_min    = traitp.bl_min(ipf);
dbh_min   = traitp.dbh_min(ipf);
dbh_adult = traitp.dbh_adult(ipf);

% There are two portions of the curve. Trees that are adult stature and
% larger and thus applicable to the published regressions, and those
% smaller than that are thus applicable to an extrapolation to known leaf
% biomass found in saplings

if(dbh>=dbh_adult)
    blmax    = d2bl1.*dbh.^d2bl2 ./ c2b;
    dblmaxdd = d2bl1.*d2bl2.*dbh.^(d2bl2-1.0) ./ c2b;
else
    bleaf_ad = d2bl1.*dbh_adult.^d2bl2 ./c2b;
    a2_small = log(bleaf_ad./bl_min)./log(dbh_adult./dbh_min);
    a1_small = bleaf_ad.*c2b/(dbh_adult.^a2_small);
    blmax    = a1_small.*dbh.^a2_small ./ c2b;
    dblmaxdd = a1_small.*a2_small.*dbh.^(a2_small-1.0) ./ c2b;
end
end

function [blmax,dblmaxdd] = dh2blmax_2pwr_spline_asym(dbh,h,traitp,ipf)

d2bl1     = traitp.d2bl1(ipf);
d2bl2     = traitp.d2bl2(ipf);
c2b       = traitp.c2b(ipf);
dbh_adult = traitp.dbh_adult(ipf);

dbh_sap   = traitp.dbh_sap(ipf);
d2bl1_sap = traitp.d2bl1_sap(ipf); %0.0201;
d2bl2_sap = traitp.d2bl2_sap(ipf); %3.1791;

[dbh_eff,ddbhedh]  = f_h2d(h,traitp,ipf);

% In this version of the 2 parameter power fit of diameter to maximum leaf
% biomass; leaf biomass is capped based on maximum height.  However, leaf
% allometry is based on diameter. So the actual asymptoted height
% calculates the diameter that would had generated that height in a
% condition where no asymptote had occured.  This is called dbh_eff.

% There are two portions of the curve. Trees that are adult stature and
% larger and thus applicable to the published regressions, and those
% smaller than that are thus applicable to an extrapolation to known leaf
% biomass found in saplings.


if(dbh_eff<dbh_sap)
    
    % Follow King's log-log linear regression
    blmax = d2bl1_sap * dbh_eff.^d2bl2_sap ./c2b;
    % Follow King's log-log linear regression
    dblmaxdd = d2bl1_sap * d2bl2_sap * dbh_eff.^(d2bl2_sap-1.0)./c2b;
    
elseif(dbh_eff>=dbh_sap && dbh_eff<dbh_adult)
    
    % Use cubic spline interolation from dbh_sap to dbh_adult
    
    blsap = d2bl1_sap * dbh_sap.^d2bl2_sap./c2b;
    blad  = d2bl1.*dbh_adult.^d2bl2./c2b;
    dblsapdd = d2bl1_sap * d2bl2_sap * dbh_sap.^(d2bl2_sap-1.0)./c2b;
    [h_adult,~] = f_h(dbh_adult,traitp,ipf);
    [~,dddh]    = f_h2d(h_adult,traitp,ipf);
    [~,dhdd]    = f_h(dbh_adult,traitp,ipf);
    ddeffdd     = dddh * dhdd;
    dbladdd = d2bl1.*d2bl2.*dbh_adult.^(d2bl2-1.0).*ddeffdd./c2b;
    [blmax,dblmaxdd] = cspline(dbh_sap,dbh_adult,blsap,blad,dblsapdd,dbladdd,dbh_eff);
    
else
    blmax       = d2bl1.*dbh_eff.^d2bl2./c2b;
    [~,dhdd]    = f_h(dbh,traitp,ipf);
    ddeffdd     = ddbhedh * dhdd;
    dblmaxdd    = d2bl1.*d2bl2.*dbh_eff.^(d2bl2-1.0).*ddeffdd./c2b;
    
end
end


function [blmax,dblmaxdd] = d2blmax_2pwr_hcap(dbh,traitp,ipf)

d2bl1     = traitp.d2bl1(ipf);
d2bl2     = traitp.d2bl2(ipf);
c2b       = traitp.c2b(ipf);
bl_min    = traitp.bl_min(ipf);
dbh_min   = traitp.dbh_min(ipf);
dbh_adult = traitp.dbh_adult(ipf);
dbh_maxh  = traitp.dbh_maxh(ipf);

if(dbh>=dbh_adult && dbh<=dbh_maxh)
    blmax    = d2bl1.*dbh.^d2bl2./c2b;
    dblmaxdd = d2bl2.*d2bl1.*dbh.^(d2bl2-1)./c2b;
elseif(dbh>dbh_maxh)
    blmax    = d2bl1.*dbh_maxh.^d2bl2./c2b;
    dblmaxdd = 0.0;
else
    bleaf_ad = d2bl1.*dbh_adult.^d2bl2 ./c2b;
    d2bl2_small = log(bleaf_ad./bl_min)./log(dbh_adult./dbh_min);
    d2bl1_small = bleaf_ad.*c2b/(dbh_adult.^d2bl2_small);
    blmax    = d2bl1_small.*dbh.^d2bl2_small ./ c2b;
    dblmaxdd = d2bl1_small.*d2bl2_small.*dbh.^(d2bl2_small-1)./ c2b;
end
end

% =========================================================================
% Diameter to height (D2H) functions
% =========================================================================

function [h,dhdd] = d2h_chave2014(dbh,traitp,ipf)

% "d2h_chave2014"
% "dbh to height via Chave et al. 2014"

% This function calculates tree height based on tree diameter and the
% environmental stress factor "E", as per Chave et al. 2015 GCB
% As opposed to previous allometric models in ED, in this formulation
% we do not impose a hard cap on tree height.  But, maximum_height
% is an important parameter, but instead of imposing a hard limit, in
% the new methodology, it will be used to trigger a change in carbon
% balance accounting.  Such that a tree that hits its maximum height will
% begin to route available NPP into seed and defense respiration.
%
% The stress function is based on the geographic location of the site.  If
% a user decides to use Chave2015 allometry, the E factor will be read in
% from a global gridded dataset and assigned for each ED patch (note it
% will be the same for each ED patch, but this distinction will help in
% porting ED into different models (patches are pure ED).  It
% assumes that the site is within the pan-tropics, and is a linear function
% of climatic water deficit, temperature seasonality and precipitation
% seasonality.  See equation 6b of Chave et al.
% The relevant equation for height in this function is 6a of the same
% manuscript, and is intended to pair with diameter to relate with
% structural biomass as per equation 7 (in which H is implicit).
%
% Chave et al. Improved allometric models to estimate the abovegroud
% biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
%
% =========================================================================

%eclim: Chave's climatological influence parameter "E"

d2h1     = traitp.d2h1(ipf);      % (alias)
d2h2     = traitp.d2h2(ipf);      % (alias)
d2h3     = traitp.d2h3(ipf);      % (alias)
eclim    = traitp.eclim(ipf);  % (alias)
dbh_maxh = traitp.dbh_maxh(ipf);
ddbh     = 0.1; % 1-mm

% For the non derivative solution, if the tree is large and
% close to any cap that is imposed, then we need to perform a
% step-integration because the asymptotic function makes the indefinite
% integral incredibly messy. Thus we use an Euler step, yes ugly,
% but it is a simple function so don't over-think it (RGK).

k      = 0.25;

if(dbh>0.5*dbh_maxh)
    
    dbh0=0.5*dbh_maxh;
    h  = exp( d2h1 - eclim + d2h2*log(dbh0) + d2h3*log(dbh0).^2.0 );
    while(dbh0<dbh)

        fl = 1./(1+exp(-k*(dbh0-dbh_maxh)));
        ae = d2h1-eclim;
        dhpdd = exp(ae).*( d2h3.*2.0.*dbh0.^(d2h2-1.0).*log(dbh0).* ...
            exp(d2h3.*log(dbh0).^2) + d2h2.*dbh0.^(d2h2-1.0).* ...
            exp(d2h3.*log(dbh0).^2.0) );
        dhdd = dhpdd*(1-fl);
        dbh0 = dbh0+ddbh;
        h    = h+ddbh*dhdd;
    end
%    display('A request for a height calculation near hmax with chave');
%    display('allometry required an explicit euler integration');
%    display('this is innefficient, and was not thought to had been');
%    display('necessary for production runs');
else
    h  = exp( d2h1 - eclim + d2h2*log(dbh) + d2h3*log(dbh).^2.0 );
end

% Deriviative

k      = 0.25;  % Stiffness of the logistic cap
% Find dbh_max where the non asymoted h crosses h_max
% Find the root for 0 = a + bx + cx^2
% where:  a = a1-E_chave-log(h_max)
%         b = a2
%         c = a3
%         dbh_max = exp(x)
% solution: x = (-(br^2 - 4*ar*cr)^(1/2)-br)/(2*cr)
%           x = (+(br^2 - 4*ar*cr)^(1/2)-br)/(2*cr)
%           x1 = exp( (-(b^2 - 4*a*c)^(1/2)-b)/(2*c));
%    dbh_maxh = exp(((d2h2^2 - ...
%        4*(-log(h_max)+d2h1-eclim)*d2h3)^(1/2)-d2h2)/(2*d2h3));
% Logistic function
fl = 1./(1+exp(-k*(dbh-dbh_maxh)));
% Derivative of logistic function wrt dbh
%dfldd = (k.*exp(-k.*(dbh+offset-dbh_max))) ...
%          ./(1+exp(-k*(dbh+offset-dbh_max))).^2;
ae = d2h1-eclim;
dhpdd = exp(ae).*( d2h3.*2.0.*dbh.^(d2h2-1.0).*log(dbh).* ...
    exp(d2h3.*log(dbh).^2) + d2h2.*dbh.^(d2h2-1.0).* ...
    exp(d2h3.*log(dbh).^2.0) );
dhdd = dhpdd*(1-fl);

end


function [h,dhdd] = d2h_poorter2006(dbh,traitp,ipf)

% "d2h_poorter2006"
% "dbh to height via Poorter et al. 2006, these routines use natively
%  asymtotic functions"
%
% Poorter et al calculated height diameter allometries over a variety of
% species in Bolivia, including those that could be classified in guilds
% that were Partial-shade-tolerants, long-lived pioneers, shade-tolerants
% and short-lived pioneers.  There results between height and diameter
% found that small stature trees had less of a tendency to asymotote in
% height and showed more linear relationships, and the largest stature
% trees tended to show non-linear relationships that asymtote.
%
% h = h_max*(1-exp(-a*dbh.^b))
%
% Poorter L, L Bongers and F Bongers.  Architecture of 54 moist-forest tree
% species: traits, trade-offs, and functional groups.  Ecology 87(5), 2006.
%
% =========================================================================

d2h1    = traitp.d2h1(ipf);  %(alias)
d2h2    = traitp.d2h2(ipf);  %(alias)
d2h3    = traitp.d2h3(ipf);  %(alias)

h = d2h1.*(1 - exp(d2h2.*dbh.^d2h3));
%h = h_max - h_max (exp(a.*dbh.^b));
%f(x) = -h_max*exp(g(x))
%g(x) = a*dbh^b
%d/dx f(g(x) = f'(g(x))*g'(x) = -a1*exp(a2*dbh^a3) * a3*a2*dbh^(a3-1)
%
%    h = -d2h3.*exp(d2h1.*dbh.^d2h2).*d2h2*d2h1*dbh.^(d2h2-1);
dhdd = -d2h1*exp(d2h2*dbh.^d2h3) .* d2h3.*d2h2.*dbh.^(d2h3-1);
end



function [h,dhdd] = d2h_2pwr(dbh,traitp,ipf)

% =========================================================================
% "d2h_2pwr"
% "dbh to height via 2 parameter power function"
% where height h is related to diameter by a linear relationship on the log
% transform where log(a) is the intercept and b is the slope parameter.
%
% log(h) = log(a) + b*log(d)
% h      = exp(log(a)) * exp(log(d))^b
% h      = a*d^b
%
% This functional form is used often in temperate allometries
% Therefore, no base reference is cited.  Although, the reader is pointed
% towards Dietze et al. 2008, King 1991, Ducey 2012 and many others for
% reasonable parameters.  Note that this subroutine is intended only for
% trees well below their maximum height, ie initialization.
%
% args
% =========================================================================
% dbh: diameter at breast height
% traitp.d2h1: the intercept parameter (however exponential of the fitted log trans)
% traitp.d2h2: the slope parameter
% return:
% h: total tree height [m]
% =========================================================================

d2h1     = traitp.d2h1(ipf);
d2h2     = traitp.d2h2(ipf);
dbh_maxh = traitp.dbh_maxh(ipf);
k      = 0.25;  % Stiffness of the logistic cap
ddbh     = 0.1; % 1-mm

% For the non derivative solution, if the tree is large and
% close to any cap that is imposed, then we need to perform a
% step-integration because the asymptotic function makes the indefinite
% integral incredibly messy. Thus we use an Euler step, yes ugly,
% but it is a simple function so don't over-think it (RGK).
if(dbh>0.5*dbh_maxh)
    dbh0=0.5*dbh_maxh;
    h = d2h1.*dbh0.^d2h2;
    while(dbh0<dbh)
        fl = 1./(1+exp(-k*(dbh0-dbh_maxh)));
        dhdd = (d2h2*d2h1).*dbh0.^(d2h2-1).*(1-fl);
        h    = h+ddbh*dhdd;
        dbh0 = dbh0+ddbh;
    end
%    display('A request for a height calculation near hmax with chave');
%    display('allometry required an explicit euler integration');
%    display('this is innefficient, and was not thought to had been');
%    display('necessary');
else
    h = d2h1.*dbh.^d2h2;
end

% The diameter at maximum height
%    dbh_maxh = (hmax./d2h1).^(1./d2h2);
% Logistic function
fl = 1./(1+exp(-k*(dbh-dbh_maxh)));
% Derivative of logistic function wrt dbh
dhdd = (d2h2*d2h1).*dbh.^(d2h2-1).*(1-fl);
end

function [h,dhdd] = d2h_obrien(dbh,traitp,ipf)

% =========================================================================
% From King et al. 1990 at BCI for saplings
% log(d) = a + b*log(h)
% d = exp(a) * h^b
% h = (1/exp(a)) * d^(1/b)
% h = d2h1*d^d2h2  where d2h1 = 1/exp(a)  d2h2 = 1/b
% d = (h/d2h1)^(1/d2h2);
% For T. tuberculata (canopy tree) a1 = -0.0704, b1 = 0.67
% =========================================================================


%a = 0.64;
%b = 0.37;
d2h1     = traitp.d2h1(ipf);
d2h2     = traitp.d2h2(ipf);
%h_max    = traitp.h_max(ipf);
dbh_hmax = traitp.dbh_maxh(ipf);

dbh_sap   = traitp.dbh_sap(ipf);
dbh_adult = traitp.dbh_adult(ipf);
d2h1_sap  = traitp.d2h1_sap(ipf);
d2h2_sap  = traitp.d2h2_sap(ipf);

% -------------------------------------------------------------------------
% Users: note that if sapling/treelet allometries are not available for
% your pft of interest, you set dbh_adult = 0 and dbh_sap = 0 in your
% parameter file.
% -------------------------------------------------------------------------

if(dbh<dbh_sap)
    
    h    = d2h1_sap * dbh^d2h2_sap;
    dhdd = d2h1_sap * d2h2_sap * dbh^(d2h2_sap-1.0);
    
elseif(dbh>=dbh_sap && dbh<dbh_adult)
    
    dhdd_sap = d2h1_sap * d2h2_sap * dbh_sap^(d2h2_sap-1.0);
    dhdd_ad  = d2h1 * 10.0^d2h2 * dbh_adult^(d2h1-1.0);
    h_sap    = d2h1_sap * dbh_sap^d2h2_sap;
    h_adult  = 10.0^(log10(dbh_adult) * d2h1 + d2h2);
    [h,dhdd] = cspline(dbh_sap,dbh_adult,h_sap,h_adult,dhdd_sap,dhdd_ad,dbh);
    
elseif(dbh>=dbh_hmax)
    h    = 10.^(log10(dbh_hmax)*d2h1+d2h2);
    dhdd = 0;
else
    h    = 10.^(log10(dbh)*d2h1+d2h2);
    dhdd = d2h1.*10.^d2h2.*dbh.^(d2h1-1.0);
end

end



% =========================================================================
% Diameter 2 above-ground biomass
% =========================================================================

function [bag,dbagdd] = dh2bag_chave2014(dbh,h,traitp,ipf)

% =========================================================================
% This function calculates tree structural biomass from tree diameter,
% height and wood density.
%
% Chave et al. Improved allometric models to estimate the abovegroud
% biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
%
% Input arguments:
% dbh: Diameter at breast height [cm]
% rho:  wood specific gravity (dry wood mass per green volume)
% height: total tree height [m]
% a1: structural biomass allometry parameter 1 (intercept)
% a2: structural biomass allometry parameter 2 (slope)
% Output:
% bag:   Total above ground biomass [kgC]
%
% =========================================================================

d2bag1 = traitp.d2bag1(ipf);
d2bag2 = traitp.d2bag2(ipf);
wood_density = traitp.wood_density(ipf);
c2b    = traitp.c2b(ipf);

bag   = (d2bag1 * (wood_density*dbh.^2*h).^d2bag2)/c2b;

[~,dhdd] = f_h(dbh,traitp,ipf);
dbagdd1  = (d2bag1.*wood_density.^d2bag2)./c2b;
dbagdd2  = d2bag2.*dbh.^(2.*d2bag2).*h.^(d2bag2-1.0).*dhdd;
dbagdd3  = h.^d2bag2.*2.*d2bag2.*dbh.^(2.*d2bag2-1);
dbagdd   = dbagdd1.*(dbagdd2 + dbagdd3);

end


function [bag,dbagdd] = d2bag_2pwr(diam,traitp,ipf)

% =========================================================================
% This function calculates tree above ground biomass according to 2
% parameter power functions. (slope and intercepts of a log-log
% diameter-agb fit:
%
% These relationships are typical for temperate/boreal plants in North
% America.  Parameters are available from Chojnacky 2014 and Jenkins 2003
%
% Note that we are using an effective diameter here, as Chojnacky 2014
% and Jenkins use diameter at the base of the plant for "woodland" species
% The diameters should be converted prior to this routine if drc.
%
% Input arguments:
% diam: effective diameter (dbh or drc) in cm
% FOR SPECIES THAT EXPECT DCM, THIS NEEDS TO BE PRE-CALCULATED!!!!
% Output:
% agb:   Total above ground biomass [kgC]
%
% =========================================================================
% Aceraceae, Betulaceae, Fagaceae and Salicaceae comprised nearly
% three-fourths of the hardwood species (Table 3)
%
% Fabaceae and Juglandaceae had specific gravities .0.60 and were
% combined, as were Hippocastanaceae and Tilaceae with specific gravities
% near 0.30. The remaining 9 families, which included mostly species with
% specific gravity 0.45–0.55, were initially grouped to construct a general
% hardwood taxon for those families having few published biomass equa-
% tions; however, 3 warranted separation, leaving 6 families for the general
% taxon.
% bag = exp(b0 + b1*ln(diameter))/c2b;
% =========================================================================

d2bag1 = traitp.d2bag1(ipf);
d2bag2 = traitp.d2bag2(ipf);
c2b    = traitp.c2b(ipf);

%max_dbh = traitp.maxdbh(ipf);
%if(diam>1.10*max_dbh)
%    display('-----------------------------------------------------');
%    display('Tree diameter is 10% larger than diameter where height');
%    display('hits maximum.  However, you specified an AGB allometry');
%    display('that does not assume capping. Please consider ');
%    display('re-evaluating your allometric assumptions, growth');
%    display('formulations or maximum height');
%    display('------------------------------------------------------');
%end

bag    = (d2bag1 * diam.^d2bag2)./c2b;
dbagdd = (d2bag2.*d2bag1.*diam.^(d2bag2-1.0))./c2b;

end

function [bag,dbagdd] = dh2bag_salda(dbh,h,traitp,ipf)

% In the interest of reducing the number of model parameters, and since
% Saldarriaga 1988 seems as though it is being deprecated, we will use
% hard-wired parameter for dh2bag_salda, which would had required 4
% variable parameters.

d2bag1       = 0.06896;
d2bag2       = 0.572;
d2bag3       = 1.94;
d2bag4       = 0.931;
c2b          = traitp.c2b(ipf);
wood_density = traitp.wood_density(ipf);

bag = d2bag1*(h^d2bag2)*(dbh^d2bag3)*(wood_density.^d2bag4)./c2b;

% bag     = a1 * h^a2 * d^a3 * r^a4
% dbag/dd = a1*r^a4 * d/dd (h^a2*d^a3)
% dbag/dd = a1*r^a4 * [ d^a3 *d/dd(h^a2) + h^a2*d/dd(d^a3) ]
% dbag/dd = a1*r^a4 * [ d^a3 * a2*h^(a2-1)dh/dd + h^a2*a3*d^(a3-1)]

term1 = d2bag1*(wood_density.^d2bag4)./c2b;
term2 = (h^d2bag2)*d2bag3*dbh^(d2bag3-1.0);
[~,dhdd] = f_h(dbh,traitp,ipf);
term3 = d2bag2*h^(d2bag2-1)*(dbh^d2bag3)*dhdd;
dbagdd   = term1*(term2+term3);

end

function [dbhe,ddbhedh] = h2d_chave2014(h,traitp,ipf)

d2h1  = traitp.d2h1(ipf);      % (alias)
d2h2  = traitp.d2h2(ipf);      % (alias)
d2h3  = traitp.d2h3(ipf);      % (alias)
eclim = traitp.eclim(ipf);     % (alias)
ar    = d2h1-eclim;

eroot = (-d2h2 + sqrt(d2h2^2 + 4*log(h.*exp(-ar))*d2h3))./(2*d2h3);
dbhe = exp(eroot);

% How about just inverting the derivative at phi? without asyms
eroot = (-d2h2 + sqrt(d2h2^2 + 4*log(h.*exp(-ar))*d2h3))./(2*d2h3);
dbh1 = exp(eroot);
dhpdd = exp(ar).*( d2h3.*2.0.*dbh1.^(d2h2-1.0).*log(dbh1).* ...
    exp(d2h3.*log(dbh1).^2) + d2h2.*dbh1.^(d2h2-1.0).* ...
    exp(d2h3.*log(dbh1).^2.0) );

ddbhedh = 1/dhpdd;

%    term1 = exp(-d2h2./(2*d2h3));
%    term2 = exp(d2h2.^2./(4.*d2h3.^2));
%    term3 = exp(-ar/d2h3);
%    term4 = h.^(1/d2h3-1.0)./(d2h3);
%    dbh   = term1*term2*term3*term4;
end


function [dbh,ddbhdh] = h2d_poorter2006(h,traitp,ipf)

% -------------------------------------------------------------------------
% Note that the height to diameter conversion in poorter is only necessary
% when initializing saplings.  In other methods, height to diameter is
% useful when defining the dbh at which point to asymptote, as maximum
% height is the user set parameter.  This function should not need to set a
% dbh_max parameter for instance, but it may end up doing so anyway, even
% if it is not used, do to poor filtering.  The poorter et al. d2h and h2d
% functions are already asymptotic, and the parameter governing maximum
% height is the d2h1 parameter.  Note as dbh gets very large, the
% exponential goes to zero, and the maximum height approaches d2h1.
% However, the asymptote has a much different shape than the logistic, so
% therefore in the Poorter et al functions, we do not set d2h1 == h_max.
% That being said, if an h_max that is greater than d2h1 is passed to this
% function, it will return a complex number. During parameter
% initialization, a check will be placed that forces:
% h_max = d2h1*0.98
% -------------------------------------------------------------------------


d2h1    = traitp.d2h1(ipf);  %(alias)
d2h2    = traitp.d2h2(ipf);  %(alias)
d2h3    = traitp.d2h3(ipf);  %(alias)

% -------------------------------------------------------------------------
% h = a1.*(1 - exp(a2.*dbh.^a3))
% h = a1 - a1*exp(a2*dbh^a3)
% a1-h = a1*exp(a2*dbh^a3)
% (a1-h)/a1 = exp(a2*dbh^a3)
% log(1-h/a1) = a2*dbh^a3
% [log(1-h/a1)/a2]^(1/a3) = dbh
%
% derivative dd/dh
% dd/dh = [log((a1-h)/a1)/a2]^(1/a3)'
%       = (1/a3)*[log((a1-h)/a1)/a2]^(1/a3-1)* [(log(a1-h)-log(a1))/a2]'
%       = (1/a3)*[log((a1-h)/a1)/a2]^(1/a3-1) * (1/(a2*(h-a1))
% dd/dh = -((log(1-h/a1)/a2).^(1/a3-1))/(a2*a3*(a1-h))
% -------------------------------------------------------------------------

dbh = (log(1.0-h/d2h1)./d2h2).^(1.0/d2h3);
ddbhdh = -((log(1-h/d2h1)/d2h2).^(1/d2h3-1))./(d2h2*d2h3*(d2h1-h));
end


function [dbh,ddbhdh] = h2d_2pwr(h,traitp,ipf)

d2h1 = traitp.d2h1(ipf);
d2h2 = traitp.d2h2(ipf);

%h = a1.*dbh.^a2;
dbh = (h/d2h1).^(1/d2h2);
%    dbh = (1/a1).^(1/a2).*h^(1/a2);
ddbhdh = (1/d2h2).*(1/d2h1).^(1/d2h2).*h^(1/d2h2-1.0);
end


function [dbh,ddbhdh] = h2d_obrien(h,traitp,ipf)

d2h1     = traitp.d2h1(ipf);
d2h2     = traitp.d2h2(ipf);
h_max    = traitp.h_max(ipf);
dbh_sap   = traitp.dbh_sap(ipf);
dbh_adult = traitp.dbh_adult(ipf);
d2h1_sap = traitp.d2h1_sap(ipf);  %1.0729;
d2h2_sap = traitp.d2h2_sap(ipf);  %1.4925;

% =========================================================================
% From King et al. 1990 at BCI for saplings
% log(d) = a + b*log(h)
% d = exp(a) * h^b
% h = (1/exp(a)) * d^(1/b)
% h = d2h1*d^d2h2  where d2h1 = 1/exp(a)  d2h2 = 1/b
% d = (h/d2h1)^(1/d2h2);
% For T. tuberculata (canopy tree) a1 = -0.0704, b1 = 0.67
% =========================================================================

h_sap   = d2h1_sap * dbh_sap^d2h2_sap;
h_adult = 10.0^(log10(dbh_adult) * d2h1 + d2h2);

if(h<h_sap)
    
    dbh    = (h/d2h1_sap)^(1.0/d2h2_sap);
    ddbhdh = (1/d2h2_sap)*(h/d2h1_sap)^(1.0/d2h2_sap-1.0);
    
elseif(h>=h_sap && h<h_adult)
    
    dddh_sap = (1.0/d2h2_sap)*(h_sap/d2h1_sap)^(1.0/d2h2_sap-1.0);
    dddh_ad  = 1.0/(d2h1*10^d2h2*dbh_adult^(d2h1-1.0));
    [dbh,ddbhdh] = cspline(h_sap,h_adult,dbh_sap,dbh_adult,dddh_sap,dddh_ad,h);
    
elseif(h<h_max)

    dbh    = 10.^((log10(h)-d2h2)./d2h1);
    ddbhdh = 1.0/(d2h1.*10.^d2h2.*dbh.^(d2h1-1.0));
else
    dbh    = 10.^((log10(h_max)-d2h2)./d2h1);
    ddbhdh = 1e20;  % Something super high, because it should be infinite
end

end


function [y,dydx] = cspline(x1,x2,y1,y2,dydx1,dydx2,x)

% ==================================================================================
% This subroutine performs a cubic spline interpolation between known
% endpoints.  The endpoints have known coordinats and slopes
% ==================================================================================

% Arguments
%real(r8),intent(in) :: x1     ! Lower endpoint independent
%real(r8),intent(in) :: x2     ! Upper endpoint independent
%real(r8),intent(in) :: y1     ! Lower endpoint dependent
%real(r8),intent(in) :: y2     ! Upper endpoint dependent
%real(r8),intent(in) :: dydx1  ! Lower endpoint slope
%real(r8),intent(in) :: dydx2  ! Upper endpoint slope
%real(r8),intent(in) :: x      ! Independent
%real(r8),intent(out) :: y     ! Dependent
%real(r8),intent(out) :: dydx  ! Slope

% Temps
%real(r8) :: t
%real(r8) :: a
%real(r8) :: b

t = (x-x1)/(x2-x1);
a = dydx1*(x2-x1) - (y2-y1);
b = -dydx2*(x2-x1) + (y2-y1);

y    = (1.0-t)*y1 + t*y2 + t*(1.0-t)*(a*(1.0-t) + b*t);
dydx = (y2-y1)/(x2-x1) + (1.0-2.0*t)*(a*(1.0-t)+b*t)/(x2-x1) + t*(1.0-t)*(b-a)/(x2-x1);

end

