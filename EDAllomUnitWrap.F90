!
! gfortran -shared -fPIC -g -o EDAllomUnitWrap.o EDAllomUnitWrap.F90
!
!
module EDAllomUnitWrap

   use iso_c_binding, only: c_double, c_int, c_char

   integer,parameter :: strlen_a = 20 ! token stinglength for allometry parameter names
   
  type, public :: EDecophyscon_type
     real(kind=8), pointer :: max_dbh            (:) ! maximum dbh at which height growth ceases... 
     real(kind=8), pointer :: wood_density       (:) ! wood density  g cm^-3  ...  
     real(kind=8), pointer :: hgt_min            (:) ! sapling height m 
     
     real(kind=8), pointer ::  c2b               (:) ! carbon to biomass ratio
     real(kind=8), pointer ::  eclim             (:) ! climatological taper influence parameter from Chave
     real(kind=8), pointer ::  bl_min            (:) ! leaf biomass at d_min [kgC]
     real(kind=8), pointer ::  h_max             (:) ! maximum possible height [m]
     real(kind=8), pointer ::  h_min             (:) ! height of plant at d_min [m]
     real(kind=8), pointer ::  slatop            (:) ! specific leaf area (SLA) at crown top [m2/gC]
     real(kind=8), pointer ::  d_adult           (:) ! diameter of an adult tree, roughly the lower diameter
                                                 ! bound for the census where the allometric relations hold
                                                 ! typically 10cm [cm]
     real(kind=8), pointer ::  d_sap             (:) ! the maximum diameter of a sapling, below which sapling
                                                 ! allometric proportions are used [cm]
     real(kind=8), pointer ::  f2l_ratio         (:) ! fine-root to leaf biomass ratio when "on-allometry" [gC/gC]
     real(kind=8), pointer ::  agb_fraction      (:) ! fraction of agb compared to total biomass [-]
     real(kind=8), pointer ::  latosa_int        (:) ! leaf area to sapwood area intercept (previously sapwood_ratio
                                                 ! [m2/cm2]
     real(kind=8), pointer ::  latosa_slp        (:) ! slope of the latosa relationship [m2/cm3]
     real(kind=8), pointer ::  d2h1_ad           (:) ! parameter 1 for adult diam to height allom
     real(kind=8), pointer ::  d2h2_ad           (:) ! parameter 2 for adult diam to height allom
     real(kind=8), pointer ::  d2h3_ad           (:) ! parameter 3 for adult diam to height allom
     real(kind=8), pointer ::  d2h1_sap          (:) ! parameter 1 for sapling diam to height allom
     real(kind=8), pointer ::  d2h2_sap          (:) ! parameter 2 for sapling diam to height allom
     real(kind=8), pointer ::  d2bl1_ad          (:) ! parameter 1 for adult diam to leaf biom allom
     real(kind=8), pointer ::  d2bl2_ad          (:) ! parameter 2 for adult diam to leaf biom allom
     real(kind=8), pointer ::  d2bl3_ad          (:) ! parameter 3 for adult diam to leaf biom allom
     real(kind=8), pointer ::  d2bl1_sap         (:) ! parameter 1 for sapling diam to leaf biom allom
     real(kind=8), pointer ::  d2bl2_sap         (:) ! parameter 2 for sapling diam to leaf boim allom
     real(kind=8), pointer ::  d2bag1            (:) ! parameter 1 for all diam to AGB allom
     real(kind=8), pointer ::  d2bag2            (:) ! parameter 2 for all diam to AGB allom
     
     character(len=strlen_a), pointer :: hallom_mode (:)  ! height allometry function type
     character(len=strlen_a), pointer :: lallom_mode (:)  ! maximum leaf allometry function type
     character(len=strlen_a), pointer :: fallom_mode (:)  ! maximum root allometry function type
     character(len=strlen_a), pointer :: aallom_mode (:)  ! AGB allometry function type
     character(len=strlen_a), pointer :: callom_mode (:)  ! coarse root allometry function type
     character(len=strlen_a), pointer :: sallom_mode (:)  ! sapwood allometry function type

  end type EDecophyscon_type
 
  type pftptr_var
     real(kind=8), dimension(:), pointer :: var_rp
     integer, dimension(:), pointer :: var_ip
     character(len = strlen_a), dimension(:), pointer :: var_cp
     character(len = strlen_a) :: var_name
     integer :: vtype
  end type pftptr_var

  type EDecophyscon_ptr
     type(pftptr_var), allocatable :: var(:)
  end type EDecophyscon_ptr
  

  type(EDecophyscon_type), public :: EDecophyscon ! ED ecophysiological constants structure
  type(EDecophyscon_ptr), public :: EDecophysptr  ! Pointer structure for obj-oriented id
  
  integer :: numparm ! Number of different PFT parameters
  integer :: numpft

contains
  

   subroutine EDEcophysconPySet(name,ipft,rval,ival,cval) 

      implicit none
      ! Arguments
      integer(c_int),intent(in) :: ipft
      character(kind=c_char), intent(in) :: name(*)
      real(c_double),intent(in) :: rval
      integer(c_int),intent(in) :: ival
      character(kind=c_char), intent(in) :: cval(*)
      ! Locals
      logical :: npfound
      integer :: ip
      
      print*,ichar(name)

!      print*,"NAME: ",name(1:2)," IPFT: ",ipft," RVAL: ",rval," IVAL: ",ival," CVAL: ",cval

      ip=0
      npfound = .true.
      do ip=1,numparm

         if (trim(name(:)) == trim(EDecophysptr%var(ip)%var_name ) ) then
            npfound = .false.
            if(EDecophysptr%var(ip)%vtype == 1) then ! real
               EDecophysptr%var(ip)%var_rp(ipft) = rval
            elseif(EDecophysptr%var(ip)%vtype == 2) then ! real
               EDecophysptr%var(ip)%var_ip(ipft) = ival
            else
               EDecophysptr%var(ip)%var_cp(ipft) = cval
            end if
         end if

         if(npfound)then
            print*,"The parameter you loaded DNE: ",name(:)
            stop
         end if

      end do
      return
   end subroutine EDEcophysconPySet


  subroutine EDEcophysconAlloc(numpft_in)
    !

    ! !ARGUMENTS:
    integer(c_int), intent(in) :: numpft_in
    ! LOCALS:
    integer                    :: iv   ! The parameter incrementer
    !------------------------------------------------------------------------

    numpft = numpft_in

    allocate( EDecophysptr%var (100) )

    iv=0

    allocate( EDecophyscon%max_dbh            (1:numpft)); EDecophyscon%max_dbh            (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "max_dbh"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%max_dbh
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%wood_density       (1:numpft)); EDecophyscon%wood_density       (:) = nan           
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "wood_density"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%wood_density
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%hgt_min            (1:numpft)); EDecophyscon%hgt_min            (:) = nan                
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "hgt_min"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%hgt_min
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%c2b                (1:numpft)); EDecophyscon%c2b                (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "c2b"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%c2b
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%eclim              (1:numpft)); EDecophyscon%eclim              (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "eclim"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%eclim
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%bl_min             (1:numpft)); EDecophyscon%bl_min             (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "bl_min"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%bl_min
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%h_max              (1:numpft)); EDecophyscon%h_max              (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "h_max"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%h_max
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%h_min              (1:numpft)); EDecophyscon%h_min              (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "h_min"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%h_min
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%slatop             (1:numpft)); EDecophyscon%slatop             (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "slatop"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%slatop
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d_adult            (1:numpft)); EDecophyscon%d_adult            (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d_adult"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d_adult
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d_sap              (1:numpft)); EDecophyscon%d_sap              (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d_sap"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d_sap
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%f2l_ratio          (1:numpft)); EDecophyscon%f2l_ratio          (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "f2l_ratio"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%f2l_ratio
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%agb_fraction       (1:numpft)); EDecophyscon%agb_fraction       (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "agb_fraction"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%agb_fraction
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%latosa_int         (1:numpft)); EDecophyscon%latosa_int         (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "eclim"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%eclim
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%latosa_slp         (1:numpft)); EDecophyscon%latosa_slp         (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "latosa_slp"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%latosa_slp
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2h1_ad            (1:numpft)); EDecophyscon%d2h1_ad            (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2h1_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2h1_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2h2_ad            (1:numpft)); EDecophyscon%d2h2_ad            (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2h2_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2h2_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2h3_ad            (1:numpft)); EDecophyscon%d2h3_ad            (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2h3_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2h3_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2h1_sap           (1:numpft)); EDecophyscon%d2h1_sap           (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2h1_sap"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2h1_sap
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2h2_sap           (1:numpft)); EDecophyscon%d2h2_sap           (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2h2_sap"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2h2_sap
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bl1_ad           (1:numpft)); EDecophyscon%d2bl1_ad           (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bl1_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bl1_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bl2_ad           (1:numpft)); EDecophyscon%d2bl2_ad           (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bl2_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bl2_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bl3_ad           (1:numpft)); EDecophyscon%d2bl3_ad           (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bl3_ad"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bl3_ad
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bl1_sap          (1:numpft)); EDecophyscon%d2bl1_sap          (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bl1_sap"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bl1_sap
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bl2_sap          (1:numpft)); EDecophyscon%d2bl2_sap          (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bl2_sap"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bl2_sap
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bag1             (1:numpft)); EDecophyscon%d2bag1             (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bag1"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bag1
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%d2bag2             (1:numpft)); EDecophyscon%d2bag2             (:) = nan
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "d2bag2"
    EDecophysptr%var(iv)%var_rp  => EDecophyscon%d2bag2
    EDecophysptr%var(iv)%vtype = 1

    allocate( EDecophyscon%hallom_mode(1:numpft)); EDecophyscon%hallom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "hallom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%hallom_mode
    EDecophysptr%var(iv)%vtype = 3

    allocate( EDecophyscon%lallom_mode(1:numpft)); EDecophyscon%lallom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "lallom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%lallom_mode
    EDecophysptr%var(iv)%vtype = 3

    allocate( EDecophyscon%fallom_mode(1:numpft)); EDecophyscon%fallom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "fallom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%fallom_mode
    EDecophysptr%var(iv)%vtype = 3

    allocate( EDecophyscon%aallom_mode(1:numpft)); EDecophyscon%aallom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "aallom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%aallom_mode
    EDecophysptr%var(iv)%vtype = 3

    allocate( EDecophyscon%callom_mode(1:numpft)); EDecophyscon%callom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "callom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%callom_mode
    EDecophysptr%var(iv)%vtype = 3

    allocate( EDecophyscon%sallom_mode(1:numpft)); EDecophyscon%sallom_mode             (:) = ""
    iv = iv + 1
    EDecophysptr%var(iv)%var_name = "sallom_mode"
    EDecophysptr%var(iv)%var_cp  => EDecophyscon%sallom_mode
    EDecophysptr%var(iv)%vtype = 3

    
    numparm = iv


    print*,"ALLOCATED ",numparm," PARAMETERS, FOR ",numpft," PFTs"


    return
 end subroutine EDEcophysconAlloc

end module EDAllomUnitWrap
