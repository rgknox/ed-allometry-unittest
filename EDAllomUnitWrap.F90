module EDAllomUnitWrap



  type, public :: EDecophyscon_type
     real(r8), pointer :: max_dbh            (:) ! maximum dbh at which height growth ceases... 
     real(r8), pointer :: freezetol          (:) ! minimum temperature tolerance... 
     real(r8), pointer :: wood_density       (:) ! wood density  g cm^-3  ...  
     real(r8), pointer :: alpha_stem         (:) ! live stem turnover rate. y-1 
     real(r8), pointer :: hgt_min            (:) ! sapling height m 
     real(r8), pointer :: cushion            (:) ! labile carbon storage target as multiple of leaf pool. 
     real(r8), pointer :: leaf_stor_priority (:) ! leaf turnover vs labile carbon use prioritisation. ! (1=lose leaves, 0=use store). 
     real(r8), pointer :: leafwatermax       (:) ! amount of water allowed on leaf   surfaces
     real(r8), pointer :: rootresist         (:)
     real(r8), pointer :: soilbeta           (:)
     real(r8), pointer :: crown              (:) ! fraction of the height of the plant that is occupied by crown. For fire model. 
     real(r8), pointer :: bark_scaler        (:) ! scaler from dbh to bark thickness. For fire model. 
     real(r8), pointer :: crown_kill         (:) ! scaler on fire death. For fire model. 
     real(r8), pointer :: initd              (:) ! initial seedling density 
     real(r8), pointer :: sd_mort            (:) ! rate of death of seeds produced from reproduction. 
     real(r8), pointer :: seed_rain          (:) ! seeds that come from outside the gridbox.  
     real(r8), pointer :: BB_slope           (:) ! ball berry slope parameter
     real(r8), pointer :: root_long          (:) ! root longevity (yrs)
     real(r8), pointer :: clone_alloc        (:) ! fraction of carbon balance allocated to clonal reproduction. 
     real(r8), pointer :: seed_alloc         (:) ! fraction of carbon balance allocated to seeds. 
     
     real(r8), pointer ::  c2b               (:) ! carbon to biomass ratio
     real(r8), pointer ::  eclim             (:) ! climatological taper influence parameter from Chave
     real(r8), pointer ::  bl_min            (:) ! leaf biomass at d_min [kgC]
     real(r8), pointer ::  h_max             (:) ! maximum possible height [m]
     real(r8), pointer ::  h_min             (:) ! height of plant at d_min [m]
     real(r8), pointer ::  slatop            (:) ! specific leaf area (SLA) at crown top [m2/gC]
     real(r8), pointer ::  d_adult           (:) ! diameter of an adult tree, roughly the lower diameter
                                                 ! bound for the census where the allometric relations hold
                                                 ! typically 10cm [cm]
     real(r8), pointer ::  d_sap             (:) ! the maximum diameter of a sapling, below which sapling
                                                 ! allometric proportions are used [cm]
     real(r8), pointer ::  f2l_ratio         (:) ! fine-root to leaf biomass ratio when "on-allometry" [gC/gC]
     real(r8), pointer ::  agb_fraction      (:) ! fraction of agb compared to total biomass [-]
     real(r8), pointer ::  latosa_int        (:) ! leaf area to sapwood area intercept (previously sapwood_ratio
                                                 ! [m2/cm2]
     real(r8), pointer ::  latosa_slp        (:) ! slope of the latosa relationship [m2/cm3]
     real(r8), pointer ::  d2h1_ad           (:) ! parameter 1 for adult diam to height allom
     real(r8), pointer ::  d2h2_ad           (:) ! parameter 2 for adult diam to height allom
     real(r8), pointer ::  d2h3_ad           (:) ! parameter 3 for adult diam to height allom
     real(r8), pointer ::  d2h1_sap          (:) ! parameter 1 for sapling diam to height allom
     real(r8), pointer ::  d2h2_sap          (:) ! parameter 2 for sapling diam to height allom
     real(r8), pointer ::  d2bl1_ad          (:) ! parameter 1 for adult diam to leaf biom allom
     real(r8), pointer ::  d2bl2_ad          (:) ! parameter 2 for adult diam to leaf biom allom
     real(r8), pointer ::  d2bl3_ad          (:) ! parameter 3 for adult diam to leaf biom allom
     real(r8), pointer ::  d2bl1_sap         (:) ! parameter 1 for sapling diam to leaf biom allom
     real(r8), pointer ::  d2bl2_sap         (:) ! parameter 2 for sapling diam to leaf boim allom
     real(r8), pointer ::  d2bag1            (:) ! parameter 1 for all diam to AGB allom
     real(r8), pointer ::  d2bag2            (:) ! parameter 2 for all diam to AGB allom
     

  end type EDecophyscon_type
  
  type(EDecophyscon_type), public :: EDecophyscon ! ED ecophysiological constants structure

contains
  
  subroutine EDecophysconInit(EDpftvarcon_inst, numpft)
    !
    ! !USES:
    use EDPftvarcon, only : EDPftvarcon_type
    !
    ! !ARGUMENTS:
    type(EDpftVarCon_type) , intent(in) :: EDpftvarcon_inst
    integer                , intent(in) :: numpft
    !
    ! !LOCAL VARIABLES:
    integer :: m, ib
    !------------------------------------------------------------------------

    allocate( EDecophyscon%max_dbh            (0:numpft)); EDecophyscon%max_dbh            (:) = nan
    allocate( EDecophyscon%freezetol          (0:numpft)); EDecophyscon%freezetol          (:) = nan
    allocate( EDecophyscon%wood_density       (0:numpft)); EDecophyscon%wood_density       (:) = nan           
    allocate( EDecophyscon%alpha_stem         (0:numpft)); EDecophyscon%alpha_stem         (:) = nan             
    allocate( EDecophyscon%hgt_min            (0:numpft)); EDecophyscon%hgt_min            (:) = nan                
    allocate( EDecophyscon%cushion            (0:numpft)); EDecophyscon%cushion            (:) = nan                
    allocate( EDecophyscon%leaf_stor_priority (0:numpft)); EDecophyscon%leaf_stor_priority (:) = nan     
    allocate( EDecophyscon%leafwatermax       (0:numpft)); EDecophyscon%leafwatermax       (:) = nan           
    allocate( EDecophyscon%rootresist         (0:numpft)); EDecophyscon%rootresist         (:) = nan             
    allocate( EDecophyscon%soilbeta           (0:numpft)); EDecophyscon%soilbeta           (:) = nan               
    allocate( EDecophyscon%crown              (0:numpft)); EDecophyscon%crown              (:) = nan                  
    allocate( EDecophyscon%bark_scaler        (0:numpft)); EDecophyscon%bark_scaler        (:) = nan            
    allocate( EDecophyscon%crown_kill         (0:numpft)); EDecophyscon%crown_kill         (:) = nan             
    allocate( EDecophyscon%initd              (0:numpft)); EDecophyscon%initd              (:) = nan                  
    allocate( EDecophyscon%sd_mort            (0:numpft)); EDecophyscon%sd_mort            (:) = nan                
    allocate( EDecophyscon%seed_rain          (0:numpft)); EDecophyscon%seed_rain          (:) = nan              
    allocate( EDecophyscon%BB_slope           (0:numpft)); EDecophyscon%BB_slope           (:) = nan               
    allocate( EDecophyscon%root_long          (0:numpft)); EDecophyscon%root_long          (:) = nan                
    allocate( EDecophyscon%seed_alloc         (0:numpft)); EDecophyscon%seed_alloc         (:) = nan               
    allocate( EDecophyscon%clone_alloc        (0:numpft)); EDecophyscon%clone_alloc        (:) = nan              
    allocate( EDecophyscon%sapwood_ratio      (0:numpft)); EDecophyscon%sapwood_ratio      (:) = nan            

    allocate( EDecophyscon%c2b                (0:numpft)); EDecophyscon%c2b                (:) = nan
    allocate( EDecophyscon%eclim              (0:numpft)); EDecophyscon%eclim              (:) = nan
    allocate( EDecophyscon%bl_min             (0:numpft)); EDecophyscon%bl_min             (:) = nan
    allocate( EDecophyscon%h_max              (0:numpft)); EDecophyscon%h_max              (:) = nan
    allocate( EDecophyscon%h_min              (0:numpft)); EDecophyscon%h_min              (:) = nan
    allocate( EDecophyscon%slatop             (0:numpft)); EDecophyscon%slatop             (:) = nan
    allocate( EDecophyscon%d_adult            (0:numpft)); EDecophyscon%d_adult            (:) = nan
    allocate( EDecophyscon%d_sap              (0:numpft)); EDecophyscon%d_sap              (:) = nan
    allocate( EDecophyscon%f2l_ratio          (0:numpft)); EDecophyscon%f2l_ratio          (:) = nan
    allocate( EDecophyscon%agb_fraction       (0:numpft)); EDecophyscon%agb_fraction       (:) = nan
    allocate( EDecophyscon%latosa_int         (0:numpft)); EDecophyscon%latosa_int         (:) = nan
    allocate( EDecophyscon%latosa_slp         (0:numpft)); EDecophyscon%latosa_slp         (:) = nan
    allocate( EDecophyscon%d2h1_ad            (0:numpft)); EDecophyscon%d2h1_ad            (:) = nan
    allocate( EDecophyscon%d2h2_ad            (0:numpft)); EDecophyscon%d2h2_ad            (:) = nan
    allocate( EDecophyscon%d2h3_ad            (0:numpft)); EDecophyscon%d2h3_ad            (:) = nan
    allocate( EDecophyscon%d2h1_sap           (0:numpft)); EDecophyscon%d2h1_sap           (:) = nan
    allocate( EDecophyscon%d2h2_sap           (0:numpft)); EDecophyscon%d2h2_sap           (:) = nan
    allocate( EDecophyscon%d2bl1_ad           (0:numpft)); EDecophyscon%d2bl1_ad           (:) = nan
    allocate( EDecophyscon%d2bl2_ad           (0:numpft)); EDecophyscon%d2bl2_ad           (:) = nan
    allocate( EDecophyscon%d2bl3_ad           (0:numpft)); EDecophyscon%d2bl3_ad           (:) = nan
    allocate( EDecophyscon%d2bl1_sap          (0:numpft)); EDecophyscon%d2bl1_sap          (:) = nan
    allocate( EDecophyscon%d2bl2_sap          (0:numpft)); EDecophyscon%d2bl2_sap          (:) = nan
    allocate( EDecophyscon%d2bag1             (0:numpft)); EDecophyscon%d2bag1             (:) = nan
    allocate( EDecophyscon%d2bag2             (0:numpft)); EDecophyscon%d2bag2             (:) = nan

    do m = 0,numpft
       EDecophyscon%max_dbh(m)               = EDPftvarcon_inst%max_dbh(m)
       EDecophyscon%freezetol(m)             = EDPftvarcon_inst%freezetol(m)
       EDecophyscon%wood_density(m)          = EDPftvarcon_inst%wood_density(m)
       EDecophyscon%alpha_stem(m)            = EDPftvarcon_inst%alpha_stem(m)
       EDecophyscon%hgt_min(m)               = EDPftvarcon_inst%hgt_min(m)
       EDecophyscon%cushion(m)               = EDPftvarcon_inst%cushion(m)
       EDecophyscon%leaf_stor_priority(m)    = EDPftvarcon_inst%leaf_stor_priority(m)
       EDecophyscon%leafwatermax(m)          = EDPftvarcon_inst%leafwatermax(m)
       EDecophyscon%rootresist(m)            = EDPftvarcon_inst%rootresist(m)
       EDecophyscon%soilbeta(m)              = EDPftvarcon_inst%soilbeta(m)
       EDecophyscon%crown(m)                 = EDPftvarcon_inst%crown(m)
       EDecophyscon%bark_scaler(m)           = EDPftvarcon_inst%bark_scaler(m)
       EDecophyscon%crown_kill(m)            = EDPftvarcon_inst%crown_kill(m)
       EDecophyscon%initd(m)                 = EDPftvarcon_inst%initd(m)
       EDecophyscon%sd_mort(m)               = EDPftvarcon_inst%sd_mort(m)
       EDecophyscon%seed_rain(m)             = EDPftvarcon_inst%seed_rain(m)
       EDecophyscon%bb_slope(m)              = EDPftvarcon_inst%bb_slope(m)
       EDecophyscon%root_long(m)             = EDPftvarcon_inst%root_long(m)
       EDecophyscon%seed_alloc(m)            = EDPftvarcon_inst%seed_alloc(m)
       EDecophyscon%clone_alloc(m)           = EDPftvarcon_inst%clone_alloc(m)
       EDecophyscon%sapwood_ratio(m)         = EDPftvarcon_inst%sapwood_ratio(m)
    end do

  end subroutine EDecophysconInit
    




end module EDAllomUnitWrap
