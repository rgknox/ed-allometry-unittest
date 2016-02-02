function [pdat] = gen_param_instance

% Manually Create a Group of Parameters

pdat.tag          = {'def-chave14-lecure', 'def-chave14-calvocarapa','def-poorter06-lescure','2par-pine',     '2par-maple',    'def-SaldaObrien' , 'King-SaldaObrien'};
pdat.h_mode       = {'chave14',            'chave14',                'poorter06',            '2par_pwr',      '2par_pwr',      'obrien'          , 'obrien'};
pdat.blmax_mode   = {'2par_pwr_asym',      '2par_pwr_asym',          '2par_pwr_asym',        'salda',         'salda',         'salda'           , 'salda'};
pdat.bfrmax_mode  = {'constant',           'constant',               'constant',             'constant',      'constant',      'constant'        , 'constant'};
pdat.bag_mode     = {'chave14',            'chave14',                'chave14',              '2par_pwr',      '2par_pwr',      'salda'           , 'salda'};
pdat.bcr_mode     = {'constant',           'constant',               'constant',             'constant',      'constant',      'constant'        , 'constant'};
pdat.bsap_mode    = {'constant',           'constant',               'constant',             'constant',      'constant',      'constant'        , 'constant'};
pdat.wood_density = [0.7,                  0.7,                      0.7,                    0.39,            0.56,            0.7               , 0.7];
pdat.latosa_int   = [0.75,                 0.75,                     0.75,                   0.75,            0.75,            0.75              , 0.75];
pdat.latosa_slp   = [0.0,                  0.0,                      0.0,                    0.0,             0.0,             0.0               , 0.0];
pdat.c2b          = [2.0,                  2.0,                      2.0,                    2.0,             2.0,             2.0               , 2.0];
pdat.eclim        = [-0.15,                0.0,                      -999,                   -999,            -999,            -999              , -999];
pdat.bl_min       = [0.4,                  0.4,                      0.4,                    0.4,             0.4,             -999               , -999];

pdat.froot_leaf   = [1.0,                  1.0,                      1.0,                    1.0,             1.0,             1.0               , 1.0];
pdat.ag_biomass   = [0.7,                  0.7,                      0.7,                    0.7,             0.7,             0.7               , 0.7];

pdat.h_max        = [35.0,                 35.0,                     61.0,                    35.0,           35.0,            35.0              , 35.0];
pdat.h_min        = [0.5,                  2.5,                      0.5,                    0.5,             0.5,             2.5               , 2.5];

pdat.slatop       = [0.024,                0.024,                    0.024,                  0.008,           0.03,            0.024             , 0.024];

pdat.dbh_adult    = [0.01,                 0.01,                     0.01,                   0.01,            0.01,            0.01,             10.0];
pdat.dbh_sap      = [0.0,                  0.0,                      0.0,                    0.0,             0.0,             0.0,              2.0];

pdat.d2h1         = [0.893,                0.893,                    61.7,                   exp(0.58119),    exp(0.31414),    0.64,             0.64];
pdat.d2h2         = [0.760,                0.760,                    -0.0352,                0.56063,         0.71813,         0.37,             0.37];
pdat.d2h3         = [-0.0340,              -0.0340,                  0.694,                  -999,            -999,            -999,             -999];

pdat.d2h1_sap     = [-999,                 -999,                     -999,                   -999,            -999,            -999,             1.0729];
pdat.d2h2_sap     = [-999,                 -999,                     -999,                   -999,            -999,            -999,             1.4925];

pdat.d2bl1        = [0.00873,              0.012,                    0.00873,                0.0419,          0.0419,          0.0419            , 0.0419];              
pdat.d2bl2        = [2.1360,               2.089,                    2.1360,                 1.56,            1.56,            1.56              , 1.56];
pdat.d2bl3        = [-999,                 -999,                     -999,                   0.55,            0.55,            0.55              , 0.55];         

pdat.d2bl1_sap    = [-999,                 -999,                     -999,                   -999,            -999,            -999,             0.0201];
pdat.d2bl2_sap    = [-999,                 -999,                     -999,                   -999,            -999,            -999,             3.1791];


pdat.d2bag1       = [0.0673,               0.0673,                  0.0673,                  exp(-2.6177),    exp(-1.801),     -999              , -999];
pdat.d2bag2       = [0.976,                0.976,                   0.976,                   2.4638,          2.3852,          -999              , -999];


end