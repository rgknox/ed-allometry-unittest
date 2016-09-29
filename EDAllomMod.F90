!===============================================================================
!
! EDAllomMod.F90
!
! A library of functions that calculate plant allometry and their
! derivatives.  Most relationships are related to diameter [cm].  All
! derivatives with respect to change in diameter have same units.
! In some cases multiple areguments will be provided, yet
! those arguments in those cases are trivially determined from diameter.
!
! Each function presented in this library is written to also return the
! derivative with respect to diameter using a logical switch "dswitch"
! (derivative-switch). With one exception, the h2d function returns the
! change in diameter with respect to height.
!
! The name convention of the functions follows the form d2...  Which
! indicates "diameter to ...".  Allometries for the following variables are
! calculated:
! h:  height [m]
! bag:  biomass above ground [kgC]  (aka AGB)
! blmax:  biomass in leaves when leaves are "on allometry"
!         this also is the "pre-trimmed" value, which is the maximum
!         or potential leaf mass a tree may have [kgC]
! bcr: biomass in coarse roots [kgC] (belowground sap+dead wood, no fines)
! bfrmax: biomass in fine roots when "on allometry" [kgC]
! bsap: biomass in sapwood (above and below) [kgC]
! bdead: biomass (above and below) in the structural pool [kgC]
!
! "on allometry" assumes the plant has leaves flushed, and has had
! sufficient carbon to meet maintenance turnover.
!
! The following model traits are used:
! hallom_mode, integer, height allometry function type
! lallom_mode, integer, maximum leaf allometry function type
! fallom_mode, integer, maximum root allometry function type
! aallom_mode, integer, AGB allometry function type
! callom_mode, integer, coarse root allometry function type
! sallom_mode, integer, sapwood allometry function type
! wood_density, real, mean stem wood specific gravity (heart,sap,bark)
! latosa_int, real, leaf area to sap area ratio, intercept [m2/cm2]
! latosa_slp, real, leaf area to sap area ratio, slope on diameter
! [m2/cm2/cm]
! c2b, real, carbon to biomass ratio (~2.0)
! d2h1_ad, real, parameter 1 for d2h allometry (intercept)
! d2h2_ad, real, parameter 2 for d2h allometry (slope)
! d2h3_ad, real, parameter 3 for d2h allometry (optional)
! eclim, real, climatological influence parameter for d2h allometry
! bl_min, real, leaf biomass of smallest tracked plant [kgC]
! l2f_ratio, real, fine root biomass per leaf biomass ratio [kgC/kgC]
! agb_fraction, real, the fraction of stem above ground [-]
! d_adult, real, the lower dbh bound where adult allometry starts (10m)
! h_max, real, maximum height of a functional type/group
! d2bl1_ad, real, parameter 1 for d2bl allometry (intercept)
! d2bl2_ad, real, parameter 2 for d2bl allometry (slope)
! d2bl3_ad, real, parameter 3 for d2bl allometry (optional)
! h_min, real, the height associated with newly recruited plant [m]
! dbh_min, real, the dbh associated with a newly recruited plant [cm]
! dbh_max, real, the diameter associated with maximum height [cm]
!               diagnosed from hmax using non-asymptotic functions
!
!===============================================================================

module EDAllomMod

! If this is a unit-test, these globals will be provided by a wrapper
#ifdef ALLOMUNITTEST
  use EDAllomUnitWrap, only: EDecophyscon !,iulog
  use iso_c_binding, only: r8 => c_double
  use iso_c_binding, only: li => c_int    ! li is "local int"
#endif
!#ifndef
!  use EDEcophysConType , only : EDecophyscon
!  use clm_varctl       , only : iulog
!  use shr_kind_mod     , only : r8 => shr_kind_r8;
!  use shr_kind_mod     , only : li => shr_kind_in;  ! li is "local int"
!#endif
  
  implicit none

  private
  public :: h2d_allom
  public :: h_allom
  public :: bag_allom
  public :: blmax_allom
  public :: bsap_allom
  public :: bcr_allom
  public :: bfrmax_allom
  public :: bdead_allom

contains


  ! ============================================================================
  ! Parameter Checks and Defaults (subroutine)
  ! ============================================================================
  
  
  


  ! ============================================================================
  ! Generic height to diameter interface
  ! ============================================================================

  subroutine h2d_allom(h,ipft,d,dddh)

    
    
    real(r8),intent(in)    :: h     ! height of plant [m]
    integer(li),intent(in) :: ipft  ! PFT index
    real(r8),intent(out)   :: d     ! plant diameter [cm]
    real(r8),intent(out)   :: dddh  ! change in diameter per height [cm/m]

    select case(EDecophyscon%hallom_mode(ipft))
    case (1) ! chave 2014
       call h2d_chave2014(h,ipft,d,dddh)
    case (2)  ! poorter 2006
       call h2d_poorter2006(h,ipft,d,dddh)
    case (3) ! 2 parameter power function
       call h2d_2pwr(h,ipft,d,dddh)
    case (4) ! Obrien et al. 199X BCI
       call h2d_obrien(h,ipft,d,dddh)
    case (5) ! Martinez-Cano
       call h2d_martcano(h,ipft,d,dddh)
    case DEFAULT
!       write(iulog,*) 'Unknown H2D Allometry: ',EDecophyscon%hallom_mode(ipft)
       stop
    end select
    return
  end subroutine h2d_allom

  ! ============================================================================
  ! Generic height interface
  ! ============================================================================

  subroutine h_allom(d,ipft,h,dhdd)

     ! Arguments
     real(r8),intent(in)     :: d     ! plant diameter [cm]
     integer(li),intent(in)  :: ipft  ! PFT index
     real(r8),intent(out)    :: h     ! plant height [m]
     real(r8),intent(out)    :: dhdd  ! change in height per diameter [m/cm]

     ! Locals
     integer                 :: hallom_mode
     real(r8)                :: h_sap
     real(r8)                :: h_ad
     real(r8)                :: dhdd_sap
     real(r8)                :: dhdd_ad
     real(r8)                :: p1
     real(r8)                :: p2
     real(r8)                :: p3
                  
     associate( &
                d_adult  => EDecophyscon%d_adult(ipft),  &
                d_sap    => EDecophyscon%d_sap(ipft),    &
                dbh_hmax => EDecophyscon%max_dbh(ipft), &
                eclim    => EDecophyscon%eclim(ipft) )
      
       ! --------------------------------------------------------------------------
       ! We may split up allometric equations by height
       ! using one method for saplings, and another for adults
       ! If the two methods are different, use a spline to interpolate between
       ! the two.  It is assumed that the ability to form a reasonable spline between
       ! the sapling diameter and the adult diameter has been checked already.
       ! --------------------------------------------------------------------------
       
       if ( d<d_sap .or. d>=d_adult) then
          
          if(d<d_sap) then
             p1 = EDecophyscon%d2h1_sap(ipft)
             p2 = EDecophyscon%d2h2_sap(ipft)
             p3 = EDecophyscon%d2h3_sap(ipft)
             hallom_mode = trim(EDecophyscon%hallom_sap_mode(ipft))
          else
             p1 = EDecophyscon%d2h1_ad(ipft)
             p2 = EDecophyscon%d2h2_ad(ipft)
             p3 = EDecophyscon%d2h3_ad(ipft)
             hallom_mode = trim(EDecophyscon%hallom_ad_mode(ipft))
          end if
          
          select case(hallom_mode)
          case (1)   ! "chave14")
             call d2h_chave2014(d,p1,p2,p3,eclim,dbh_hmax,h,dhdd)
          case (2)   ! "poorter06"
             call d2h_poorter2006(d,p1,p2,p3,h,dhdd)
          case (3)   ! "2parameter power function h=a*d^b "
             call d2h_2pwr(d,p1,p2,dbh_hmax,h,dhdd)
          case (4)   ! "obrien"
             call d2h_obrien(d,p1,p2,p3,dbh_hmax,h,dhdd)
          case (5)   ! Martinez-Cano
             call d2h_martcano(d,p1,p2,p3,h,dhdd)
          case DEFAULT
             stop
          end select
          
       else
          
          ! ---------------------------------------------------------------------
          ! Interpolate between the two methods using a spline
          ! ---------------------------------------------------------------------
          
          ! Calculate the boundaries for the sapling
          p1 = EDecophyscon%d2h1_sap(ipft)
          p2 = EDecophyscon%d2h2_sap(ipft)
          p3 = EDecophyscon%d2h3_sap(ipft)
          hallom_mode = trim(EDecophyscon%hallom_sap_mode(ipft))
          
          select case(hallom_mode)
          case (1)   !"chave14") 
             call d2h_chave2014(d_sap,p1,p2,p3,dbh_hmax,h_sap,dhdd_sap)
          case (2)   ! "poorter06")
             call d2h_poorter2006(d_sap,p1,p2,p3,h_sap,dhdd_sap)
          case (3) ! "2par_pwr")
             call d2h_2pwr(d_sap,p1,p2,dbh_hmax,h_sap,dhdd_sap)
          case (4) ! "obrien")
             call d2h_obrien(d_sap,p1,p2,p3,dbh_hmax,h_sap,dhdd_sap)
          case (5) ! Martinez-Cano
             call d2h_martcano(d_sap,p1,p2,p3,h_sap,dhdd_sap)
          case DEFAULT
             stop
          end select
          
          ! Calculate the boundaries for the adult
          p1 = EDecophyscon%d2h1_ad(ipft)
          p2 = EDecophyscon%d2h2_ad(ipft)
          p3 = EDecophyscon%d2h3_ad(ipft)
          hallom_mode = trim(EDecophyscon%hallom_ad_mode(ipft))
          
          select case(hallom_mode)
          case (1)   !"chave14") 
             call d2h_chave2014(d_adult,p1,p2,p3,dbh_hmax,h_ad,dhdd_ad)
          case (2)   ! "poorter06")
             call d2h_poorter2006(d_adult,p1,p2,p3,h_ad,dhdd_ad)
          case (3) ! "2par_pwr")
             call d2h_2pwr(d_adult,p1,p2,dbh_hmax,h_ad,dhdd_ad)
          case (4) ! "obrien")
             call d2h_obrien(d_adult,p1,p2,p3,dbh_hmax,h_ad,dhdd_ad)
          case (5) ! Martinez-Cano
             call d2h_martcano(d_adult,p1,p2,p3,h_ad,dhdd_ad)
          case DEFAULT
             stop
          end select
          
          ! Use cubic spline interolation from d_sap to d_adult
          call cspline(d_sap,d_adult,h_sap,h_ad,dhdd_sap,dhdd_ad,d,h,dhdd)
       end if
       
       
     end associate
     return
  end subroutine h_allom
  
  ! ============================================================================
  ! Generic AGB interface
  ! ============================================================================
  
  subroutine bag_allom(d,h,ipft,bag,dbagdd)
    
    
    real(r8),intent(in)    :: d       ! plant diameter [cm]
    real(r8),intent(in)    :: h       ! plant height [m]
    integer(li),intent(in) :: ipft    ! PFT index
    real(r8),intent(out)   :: bag     ! plant height [m]
    real(r8),intent(out)   :: dbagdd  ! change in agb per diameter [kgC/cm]

    select case(EDecophyscon%aallom_mode(ipft))
    case (1) !"chave14") 
       call dh2bag_chave2014(d,h,ipft,bag,dbagdd)
    case (2) !"2par_pwr")
        ! Switch for woodland dbh->drc
       call d2bag_2pwr(d,ipft,bag,dbagdd)
    case (3) !"salda")
       call dh2bag_salda(d,h,ipft,bag,dbagdd)
    case DEFAULT
!       write(iulog,*) 'Unknown D-2-BAG Allometry: ', &
!            EDecophyscon%aallom_mode(ipft)
       stop
    end select
    return
  end subroutine bag_allom

  ! ============================================================================
  ! Generic diameter to maximum leaf biomass interface
  ! ============================================================================

  subroutine blmax_allom(d,h,ipft,blmax,dblmaxdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out)   :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    
    select case(EDecophyscon%lallom_mode(ipft))
    case(1) !"salda")
       call d2blmax_salda(d,ipft,blmax,dblmaxdd)
    case(2) !"2par_pwr")
       call d2blmax_2pwr(d,ipft,blmax,dblmaxdd)
    case(3) !"2par_pwr_asym")
       call dh2blmax_2pwr_spline_asym(d,h,ipft,blmax,dblmaxdd)
    case(4) !"2par_pwr_hcap")
       call d2blmax_2pwr_hcap(d,ipft,blmax,dblmaxdd)
    case DEFAULT
!       write(iulog,*) 'Unknown D-2-BLMAX Allometry: ', &
!            EDecophyscon%lallom_mode(ipft)
       stop
    end select
    return
  end subroutine blmax_allom

  ! ============================================================================
  ! Generic sapwood biomass interface
  ! ============================================================================
  
  subroutine bsap_allom(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! chage in blmax per diam [kgC/cm]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out)   :: dbsapdd   ! change leaf bio per d [kgC/cm]

    select case(EDecophyscon%sallom_mode(ipft))
       ! ---------------------------------------------------------------------
       ! Currently both sapwood area proportionality methods use the same
       ! machinery.  The only differences are related to the parameter
       ! checking at the beginning.  For constant proportionality, the slope
       ! of the la:sa to diameter line is zero.
       ! ---------------------------------------------------------------------
    case(1,2) !"constant","dlinear")
       call bsap_dlinear(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)
    case DEFAULT
!       write(iulog,*) 'Unknown D-2-BLMAX Allometry: ', &
!            EDecophyscon%lallom_mode(ipft)
       stop
    end select
    return
  end subroutine bsap_allom

  ! ============================================================================
  ! Generic coarse root biomass interface
  ! ============================================================================
  
  subroutine bcr_allom(d,bag,dbagdd,ipft,bcr,dbcrdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: bag       ! above ground biomass [kgC]
    real(r8),intent(in)    :: dbagdd    ! change in agb per diameter [kgC/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bcr       ! coarse root biomass [kgC]
    real(r8),intent(out)   :: dbcrdd    ! change croot bio per diam [kgC/cm]

    select case(EDecophyscon%callom_mode(ipft))
    case(1) !"constant")
       call bcr_const(d,bag,dbagdd,ipft,bcr,dbcrdd)
    case DEFAULT
!       write(iulog,*) 'Unknown D-2-BCR Allometry: ', &
!            EDecophyscon%callom_mode(ipft)
       stop
    end select
    return
  end subroutine bcr_allom

  ! ============================================================================
  ! Generic maximum fine root biomass interface
  ! ============================================================================

  subroutine bfrmax_allom(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: blmax     ! max leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bfrmax    ! max fine-root root biomass [kgC]
    real(r8),intent(out)   :: dbfrmaxdd ! change frmax bio per diam [kgC/cm]
    
    select case(EDecophyscon%fallom_mode(ipft))
    case(1) ! "constant")
       call bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
    case DEFAULT
!       write(iulog,*) 'Unknown D-2-BFRMAX Allometry: ', &
!            EDecophyscon%fallom_mode(ipft)
       stop
    end select
    return
  end subroutine bfrmax_allom

  ! ============================================================================
  ! Dead biomass interface
  ! ============================================================================

  subroutine bdead_allom(bag,bcr,blmax,bsap,dbagdd,dbcrdd,dblmaxdd,dbsapdd, &
       bdead,dbdeaddd)

    
    real(r8),intent(in)  :: bag       ! agb [kgC]
    real(r8),intent(in)  :: bcr       ! coarse root biomass [kgC]
    real(r8),intent(in)  :: blmax     ! max leaf biomass [kgC]
    real(r8),intent(in)  :: bsap      ! sapwood biomass [kgC]

    real(r8),intent(in)  :: dbagdd    ! change in agb per d [kgC/cm]
    real(r8),intent(in)  :: dbcrdd    ! change in croot per d [kgC/cm]
    real(r8),intent(in)  :: dblmaxdd  ! change in blmax per d [kgC/cm]
    real(r8),intent(in)  :: dbsapdd   ! change in bsap per d [kgC/cm]
    
    real(r8),intent(out) :: bdead     ! dead biomass (heartw/struct) [kgC]
    real(r8),intent(out) :: dbdeaddd  ! change in bdead per d [kgC/cm]

    ! bdead is diagnosed as the mass balance from all other pools
    ! and therefore, no options are necessary
    ! We are ignoring blmax right now, because it is insignificant in large
    ! trees and may cause negatives in treelets and saplings
    
    !bdead = bag+bcr-blmax-bsap
    bdead = bag+bcr-bsap
    dbdeaddd = dbagdd+dbcrdd-dbsapdd

    return
  end subroutine bdead_allom
  
  ! ============================================================================
  ! Specific bfrmax relationships
  ! ============================================================================
  
  subroutine bfrmax_const(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: blmax     ! max leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bfrmax    ! max fine-root root biomass [kgC]
    real(r8),intent(out)   :: dbfrmaxdd ! change frmax bio per diam [kgC/cm]
    
    associate( f2l_ratio => EDecophyscon%f2l_ratio(ipft) ) 
      
      bfrmax = blmax*f2l_ratio
      
      ! dbfr/dd = dbfrmax/dblmax * dblmax/dd
      dbfrmaxdd = f2l_ratio*dblmaxdd
    end associate
    return
  end subroutine bfrmax_const


  ! ============================================================================
  ! Specific bcr relationships
  ! ============================================================================
  
  subroutine bcr_const(d,bag,dbagdd,ipft,bcr,dbcrdd)
    
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: bag       ! above ground biomass [kg]
    real(r8),intent(in)    :: dbagdd    ! change in agb per diameter [kg/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bcr       ! coarse root biomass [kg]
    real(r8),intent(out)   :: dbcrdd    ! change croot bio per diam [kg/cm]

    associate( agb_fraction => EDecophyscon%agb_fraction(ipft) )

      ! btot = bag + bcr
      ! bag = btot*agb_fraction
      ! bag/agb_fraction = bag + bcr
      ! bcr = bag*(1/agb_fraction-1)
      bcr = bag*(1.0_r8/agb_fraction-1.0_r8)

      ! Derivative
      ! dbcr/dd = dbcr/dbag * dbag/dd
      dbcrdd = (1.0_r8/agb_fraction-1.0_r8)*dbagdd
    end associate
    return
  end subroutine bcr_const

  
  ! ============================================================================
  ! Specific d2bsap relationships
  ! ============================================================================

  subroutine bsap_dlinear(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)

    ! -------------------------------------------------------------------------
    ! Calculate sapwood biomass based on leaf area to sapwood area
    ! proportionality.  In this function, the leaftosapwood area is a function
    ! of plant size, see Calvo-Alvarado and Bradley Christoferson
    ! In this case: parameter latosa (from constant proportionality)
    !   is the intercept of the diameter function.
    !
    ! For very small plants, the fraction can get very large, so cap the amount
    ! of sapwood at X! of agb-bleaf
    ! -------------------------------------------------------------------------

    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    real(r8),intent(in)    :: h         ! plant height [m]
    real(r8),intent(in)    :: blmax     ! plant leaf biomass [kgC]
    real(r8),intent(in)    :: dblmaxdd  ! change in blmax per diam [kgC/cm]
    real(r8),intent(in)    :: dhdd      ! change in height per diameter [m/cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: bsap      ! plant leaf biomass [kgC]
    real(r8),intent(out)   :: dbsapdd   ! change leaf bio per diameter [kgC/cm]

    real(r8)               :: latosa    ! m2 leaf area per cm2 sap area
    real(r8)               :: hbl2bsap  ! sapwood biomass per lineal height
                                          ! and kg of leaf
    real(r8)               :: bag       ! aboveground biomass [kgC]
    real(r8)               :: dbagdd    ! change in agb per diam [kgC/cm]

    ! Constrain sapwood to be no larger than 75% of total agb
    real(r8),parameter :: max_agbfrac = 0.75_r8 
    real(r8),parameter :: gtokg  = 1000.0_r8   ! gram to kg conversion
    real(r8),parameter :: cm2tom2 = 10000.0_r8 ! cm**2 to m**2 conversion
    real(r8),parameter :: mg2kg   = 1000.0_r8  ! Mg to kg conversion [kg/Mg]

    associate ( latosa_int => EDecophyscon%latosa_int(ipft), &
         latosa_slp => EDecophyscon%latosa_slp(ipft), &
         sla        => EDecophyscon%slatop(ipft), &
         wood_density => EDecophyscon%wood_density(ipft), &
         c2b          => EDecophyscon%c2b(ipft))

      ! ------------------------------------------------------------------------
      ! Calculate sapwood biomass per linear height and kgC of leaf [m-1]
      ! Units: 
      ! (1/latosa)* slatop*    gtokg    *   cm2tom2     / c2b   * mg2kg  * dens
      ! [cm2/m2]*[m2/gC]*[1000gC/1kgC]*[1m2/10000cm2] /[kg/kgC]*[kg/Mg]*[Mg/m3]
      !            ->[cm2/gC]
      !                      ->[cm2/kgC]
      !                                   ->[m2/kgC]
      !                                                 ->[m2/kg]
      !                                                          ->[m2/Mg]
      !                                                                  ->[/m]
      ! ------------------------------------------------------------------------
      
      latosa = latosa_int + d*latosa_slp

      hbl2bsap = sla*gtokg*wood_density*mg2kg/(latosa*c2b*cm2tom2)

      call bag_allom(d,h,ipft,bag,dbagdd)

      bsap = min(max_agbfrac*bag,hbl2bsap * h * blmax)

      ! Derivative
      ! dbldmaxdd is deriv of blmax wrt dbh (use directives to check oop)
      ! dhdd is deriv of height wrt dbh (use directives to check oop)
      
      if (bsap<max_agbfrac*bag) then
         dbsapdd = hbl2bsap*(h*dblmaxdd + blmax*dhdd)
      else
         dbsapdd = max_agbfrac*dbagdd
      end if
    end associate
    return
  end subroutine bsap_dlinear

  ! ============================================================================
  ! Specific d2blmax relationships
  ! ============================================================================
  
  subroutine d2blmax_salda(d,ipft,blmax,dblmaxdd)
    
    
    real(r8),intent(in)    :: d         ! plant diameter [cm]
    integer(li),intent(in) :: ipft      ! PFT index
    real(r8),intent(out)   :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out)   :: dblmaxdd  ! change leaf bio per diam [kgC/cm]
    real(r8) :: blad,blsap,dbladdd,dblsapdd

    associate( &
         d2bl1_ad    => EDecophyscon%d2bl1_ad(ipft),  &   !0.0419
         d2bl2_ad    => EDecophyscon%d2bl2_ad(ipft),  &    !1.56
         d2bl3_ad    => EDecophyscon%d2bl3_ad(ipft),  &   !0.55
         rho         => EDecophyscon%wood_density(ipft), &
         dbh_maxh => EDecophyscon%max_dbh(ipft),     &
         c2b       => EDecophyscon%c2b(ipft),         &
         d_adult => EDecophyscon%d_adult(ipft),       &
         d_sap   => EDecophyscon%d_sap(ipft),         &
         d2bl1_sap => EDecophyscon%d2bl1_sap(ipft),   &   !0.0201
         d2bl2_sap => EDecophyscon%d2bl2_sap(ipft)) !3.1791

      ! ======================================================================
      ! From King et al. 1990 at BCI for saplings
      !
      ! log(bl) = a2 + b2*log(h)
      ! bl = exp(a2) * h**b2
      !
      ! and:
      !
      ! log(d) = a1 + b1*log(h)
      ! d = exp(a1) * h**b1
      ! h = (1/exp(a1)) * d**(1/b1)
      !
      ! bl = exp(a2) * ((1/exp(a1)) * d**(1/b1))**b2
      ! bl = exp(a2) * (1/exp(a1))**b2 * d**(b2/b1)
      ! bl = d2bl1_ad * d ** d2bl2_ad
      ! where: d2bl1_ad = exp(a2) * (1/exp(a1))**b2
      !        d2bl2_ad = (b2/b1)
      ! For T. tuberculata (canopy tree):
      ! a1 = -0.0704, b1 = 0.67
      ! a2 = -4.056,  b2 = 2.13
      ! ======================================================================

      ! Saldarriaga has a height cap on leaf biomass
      if (d<d_sap) then
         
         ! Follow King's log-log linear regression
         blmax    = d2bl1_sap * d**d2bl2_sap /c2b
         dblmaxdd = d2bl1_sap * d2bl2_sap * d**(d2bl2_sap-1.0_r8)/c2b
      
      elseif (d>=d_sap .and. d<d_adult) then
    
         ! Use cubic spline interolation from d_sap to d_adult
         blsap = d2bl1_sap * d_sap**d2bl2_sap/c2b
         blad  = d2bl1_ad * d_adult**d2bl2_ad * rho**d2bl3_ad
         dblsapdd = d2bl1_sap * d2bl2_sap * d_sap**(d2bl2_sap-1.0_r8)/c2b
         dbladdd  = d2bl1_ad * d2bl2_ad * &
              d_adult**(d2bl2_ad-1.0_r8) * rho**d2bl3_ad
    
         call cspline(d_sap,d_adult,blsap,blad,dblsapdd,dbladdd,d,blmax,dblmaxdd)
    
      elseif(d>=d_adult .and. d<dbh_maxh) then
         blmax = d2bl1_ad * d**d2bl2_ad * rho**d2bl3_ad
         dblmaxdd = d2bl1_ad*d2bl2_ad * d**(d2bl2_ad-1.0_r8) * rho**d2bl3_ad
      else
         blmax    = d2bl1_ad * dbh_maxh**d2bl2_ad * rho**d2bl3_ad
         dblmaxdd = 0.0
      end if
      return
    end associate
  end subroutine d2blmax_salda
  
  subroutine d2blmax_2pwr(d,ipft,blmax,dblmaxdd)
    
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    integer(li),intent(in)       :: ipft      ! PFT index
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out) :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    real(r8) :: a1_small,a2_small,bleaf_ad

    associate( &
         d2bl1_ad  => EDecophyscon%d2bl1_ad(ipft), &
         d2bl2_ad  => EDecophyscon%d2bl2_ad(ipft), &
         c2b       => EDecophyscon%c2b(ipft), &
         bl_min    => EDecophyscon%bl_min(ipft), &
         dbh_min   => EDecophyscon%dbh_min(ipft), &
         d_adult   => EDecophyscon%d_adult(ipft))
      
      ! There are two portions of the curve. Trees that are adult stature and
      ! larger and thus applicable to the published regressions, and those
      ! smaller than that are thus applicable to an extrapolation to known leaf
      ! biomass found in saplings
      
      if (d>=d_adult) then
         blmax    = d2bl1_ad*d**d2bl2_ad / c2b
         dblmaxdd = d2bl1_ad*d2bl2_ad *d**(d2bl2_ad-1.0_r8) / c2b
      else
         bleaf_ad = d2bl1_ad*d_adult**d2bl2_ad / c2b
         a2_small = log(bleaf_ad/bl_min)/log(d_adult/dbh_min)
         a1_small = bleaf_ad*c2b/(d_adult**a2_small)
         blmax    = a1_small*d**a2_small / c2b
         dblmaxdd = a1_small*a2_small*d**(a2_small-1.0_r8) / c2b
      end if
      return
    end associate
  end subroutine d2blmax_2pwr

  subroutine dh2blmax_2pwr_spline_asym(d,h,ipft,blmax,dblmaxdd)

    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    real(r8),intent(in)  :: h         ! plant height [m]
    integer(li),intent(in)       :: ipft      ! PFT index
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out) :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]
    real(r8) :: dbh_eff
    real(r8) :: ddbhedh,ddedh,ddeffdd 
    real(r8) :: blsap,blad,dblsapdd,dbladdd
    real(r8) :: h_adult,dh_adultdd
    real(r8) :: dd_adultdh,dhdd
    real(r8) :: dj,hj ! junk variables

    associate( &
         d2bl1_ad  => EDecophyscon%d2bl1_ad(ipft), &
         d2bl2_ad  => EDecophyscon%d2bl2_ad(ipft), &
         c2b       => EDecophyscon%c2b(ipft), &
         d_adult   => EDecophyscon%d_adult(ipft), &
         d_sap     => EDecophyscon%d_sap(ipft), &
         d2bl1_sap => EDecophyscon%d2bl1_sap(ipft), & !0.0201;
         d2bl2_sap => EDecophyscon%d2bl2_sap(ipft)) !3.1791;)

      call h2d_allom(h,ipft,dbh_eff,ddbhedh)

      ! In this version of the 2 parameter power fit of diameter to maximum leaf
      ! biomass; leaf biomass is capped based on maximum height.  However, leaf
      ! allometry is based on diameter. So the actual asymptoted height
      ! calculates the diameter that would had generated that height in a
      ! condition where no asymptote had occured.  This is called dbh_eff.
      
      ! There are two portions of the curve. Trees that are adult stature and
      ! larger and thus applicable to the published regressions, and those
      ! smaller than that are thus applicable to an extrapolation to known leaf
      ! biomass found in saplings.
      

      if (dbh_eff<d_sap) then
    
         ! Follow King's log-log linear regression
         blmax = d2bl1_sap * dbh_eff**d2bl2_sap /c2b
         ! Follow King's log-log linear regression
         dblmaxdd = d2bl1_sap * d2bl2_sap * dbh_eff**(d2bl2_sap-1.0_r8)/c2b
         
      elseif (dbh_eff>=d_sap .and. dbh_eff<d_adult) then
    
         ! Use cubic spline interolation from d_sap to d_adult
         blsap = d2bl1_sap * d_sap**d2bl2_sap/c2b
         blad  = d2bl1_ad * d_adult**d2bl2_ad/c2b
         dblsapdd = d2bl1_sap * d2bl2_sap * d_sap**(d2bl2_sap-1.0_r8)/c2b

         call h_allom(d_adult,ipft,h_adult,dh_adultdd)
         call h2d_allom(h_adult,ipft,dj,dd_adultdh)

         ddeffdd = dd_adultdh * dh_adultdd
         dbladdd = d2bl1_ad*d2bl2_ad*d_adult**(d2bl2_ad-1.0_r8)*ddeffdd/c2b
         
         call cspline(d_sap,d_adult,blsap,blad,dblsapdd, &
              dbladdd,dbh_eff,blmax,dblmaxdd)

      else

         blmax       = d2bl1_ad*dbh_eff**d2bl2_ad/c2b
         call h_allom(d,ipft,hj,dhdd)
         call h2d_allom(hj,ipft,dj,ddedh)
         ddeffdd     = ddedh * dhdd
         dblmaxdd    = d2bl1_ad*d2bl2_ad*dbh_eff**(d2bl2_ad-1.0_r8)*ddeffdd/c2b
    
      end if
      return
    end associate
  end subroutine dh2blmax_2pwr_spline_asym

  subroutine d2blmax_2pwr_hcap(d,ipft,blmax,dblmaxdd)
    
    
    real(r8),intent(in)  :: d         ! plant diameter [cm]
    integer(li),intent(in)       :: ipft      ! PFT index
    real(r8),intent(out) :: blmax     ! plant leaf biomass [kg]
    real(r8),intent(out) :: dblmaxdd  ! change leaf bio per diameter [kgC/cm]

    real(r8) :: bleaf_ad,d2bl2_ad_small,d2bl1_ad_small

    associate( &
         d2bl1_ad => EDecophyscon%d2bl1_ad(ipft), &
         d2bl2_ad => EDecophyscon%d2bl2_ad(ipft), &
         c2b      => EDecophyscon%c2b(ipft), &
         bl_min   => EDecophyscon%bl_min(ipft), &
         dbh_min  => EDecophyscon%dbh_min(ipft), &
         d_adult  => EDecophyscon%d_adult(ipft), &
         dbh_maxh => EDecophyscon%max_dbh(ipft) )

      if ( d>=d_adult .and. d<=dbh_maxh) then
         blmax    = d2bl1_ad*d**d2bl2_ad/c2b
         dblmaxdd = d2bl2_ad*d2bl1_ad*d**(d2bl2_ad-1)/c2b
      elseif (d>dbh_maxh) then
         blmax    = d2bl1_ad*dbh_maxh**d2bl2_ad/c2b
         dblmaxdd = 0.0
      else
         bleaf_ad = d2bl1_ad*d_adult**d2bl2_ad /c2b
         d2bl2_ad_small = log(bleaf_ad/bl_min)/log(d_adult/dbh_min)
         d2bl1_ad_small = bleaf_ad*c2b/(d_adult**d2bl2_ad_small)
         blmax    = d2bl1_ad_small*d**d2bl2_ad_small/c2b
         dblmaxdd = d2bl1_ad_small*d2bl2_ad_small*d**(d2bl2_ad_small-1)/c2b
      end if
      return
    end associate
  end subroutine d2blmax_2pwr_hcap

  ! =========================================================================
  ! Diameter to height (D2H) functions
  ! =========================================================================

  subroutine d2h_chave2014(d,p1,p2,p3,eclim,dbh_hmax,h,dhdd)

    ! "d2h_chave2014"
    ! "d to height via Chave et al. 2014"
    
    ! This function calculates tree height based on tree diameter and the
    ! environmental stress factor "E", as per Chave et al. 2015 GCB
    ! As opposed to previous allometric models in ED, in this formulation
    ! we do not impose a hard cap on tree height.  But, maximum_height
    ! is an important parameter, but instead of imposing a hard limit, in
    ! the new methodology, it will be used to trigger a change in carbon
    ! balance accounting.  Such that a tree that hits its maximum height will
    ! begin to route available NPP into seed and defense respiration.
    !
    ! The stress function is based on the geographic location of the site.  If
    ! a user decides to use Chave2015 allometry, the E factor will be read in
    ! from a global gridded dataset and assigned for each ED patch (note it
    ! will be the same for each ED patch, but this distinction will help in
    ! porting ED into different models (patches are pure ED).  It
    ! assumes that the site is within the pan-tropics, and is a linear function
    ! of climatic water deficit, temperature seasonality and precipitation
    ! seasonality.  See equation 6b of Chave et al.
    ! The relevant equation for height in this function is 6a of the same
    ! manuscript, and is intended to pair with diameter to relate with
    ! structural biomass as per equation 7 (in which H is implicit).
    !
    ! Chave et al. Improved allometric models to estimate the abovegroud
    ! biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
    !
    ! =========================================================================
    
    !eclim: Chave's climatological influence parameter "E"

    
    real(r8),intent(in)  :: d     ! plant diameter [cm]
    real(r8),intent(in)  :: p1       ! parameter a 
    real(r8),intent(in)  :: p2       ! parameter b
    real(r8),intent(in)  :: p3       ! parameter c
    real(r8),intent(in)  :: eclim    ! climatological parameter "E"
    real(r8),intent(in)  :: dbh_hmax ! dbh at maximum height [cm]

    real(r8),intent(out) :: h     ! plant height [m]
    real(r8),intent(out) :: dhdd  ! change in height per diameter [m/cm]

    real(r8) :: dbh0,fl,ae
    real(r8) :: dhpdd

    real(r8),parameter :: ddbh = 0.1_r8 ! 1-mm 
    real(r8),parameter :: k    = 0.25_r8

    ! For the non derivative solution, if the tree is large and
    ! close to any cap that is imposed, then we need to perform a
    ! step-integration because the asymptotic function makes the indefinite
    ! integral incredibly messy. Thus we use an Euler step, yes ugly,
    ! but it is a simple function so don't over-think it
    
    if (d>0.5_r8*dbh_maxh) then
       dbh0=0.5_r8*dbh_maxh
       h  = exp( p1 - eclim + p2*log(dbh0) + p3*log(dbh0)**2.0_r8 )
       do while (dbh0<d) 
          fl = 1.0_r8/(1.0_r8+exp(-k*(dbh0-dbh_maxh)))
          ae = p1-eclim
          dhpdd = exp(ae)* & 
                ( p3*2.0_r8*dbh0**(p2-1.0_r8)*log(dbh0)* &
                exp(p3*log(dbh0)**2.0_r8) + p2*dbh0**(p2-1.0_r8)* &
                exp(p3*log(dbh0)**2.0_r8) )
          dhdd = dhpdd*(1.0_r8-fl)
          dbh0 = dbh0+ddbh
          h    = h+ddbh*dhdd
       end do
       
         !    display("A request for a height calculation near hmax with chave")
         !    display("allometry required an explicit euler integration")
         !    display("this is innefficient, and was not thought to had been")
         !    display("necessary for production runs")
      else
         h  = exp( p1 - eclim + p2*log(d) + p3*log(d)**2.0 )
      end if

      ! Deriviative

      ! Find dbh_max where the non asymoted h crosses h_max
      ! Find the root for 0 = a + bx + cx**2
      ! where:  a = a1-E_chave-log(h_max)
      !         b = a2
      !         c = a3
      !         dbh_max = exp(x)
      ! solution: x = (-(br**2 - 4*ar*cr)**(1/2)-br)/(2*cr)
      !           x = (+(br**2 - 4*ar*cr)**(1/2)-br)/(2*cr)
      !           x1 = exp( (-(b**2 - 4*a*c)**(1/2)-b)/(2*c))
      !    dbh_maxh = exp(((p2**2 - ...
      !        4*(-log(h_max)+p1-eclim)*p3)**(1/2)-p2)/(2*p3))
      ! Logistic function
      fl = 1.0_r8/(1.0_r8+exp(-k*(d-dbh_maxh)))
      ! Derivative of logistic function wrt d
      !dfldd = (k.*exp(-k.*(dbh+offset-dbh_max))) ...
      !          /(1+exp(-k*(dbh+offset-dbh_max)))**2
      ae = p1-eclim
      dhpdd = exp(ae)*( p3*2.0_r8*d**(p2-1.0_r8)*log(d)* &
            exp(p3*log(d)**2.0_r8) + p2*d**(p2-1.0_r8)* &
            exp(p3*log(d)**2.0_r8) )
      dhdd = dhpdd*(1.0_r8-fl)
      
      return

  end subroutine d2h_chave2014

  subroutine d2h_poorter2006(d,p1,p2,p3,h,dhdd)
    
    ! "d2h_poorter2006"
    ! "d to height via Poorter et al. 2006, these routines use natively
    !  asymtotic functions"
    !
    ! Poorter et al calculated height diameter allometries over a variety of
    ! species in Bolivia, including those that could be classified in guilds
    ! that were Partial-shade-tolerants, long-lived pioneers, shade-tolerants
    ! and short-lived pioneers.  There results between height and diameter
    ! found that small stature trees had less of a tendency to asymotote in
    ! height and showed more linear relationships, and the largest stature
    ! trees tended to show non-linear relationships that asymtote.
    !
    ! h = h_max*(1-exp(-a*d**b))
    !
    ! Poorter L, L Bongers and F Bongers.  Architecture of 54 moist-forest tree
    ! species: traits, trade-offs, and functional groups.  Ecology 87(5), 2006.
    !
    ! =========================================================================

    
    real(r8),intent(in)  :: d     ! plant diameter [cm] 
    real(r8),intent(in)     :: p1       ! parameter a = h_max
    real(r8),intent(in)     :: p2       ! parameter b
    real(r8),intent(in)     :: p3       ! parameter c
    real(r8),intent(out) :: h     ! plant height [m]
    real(r8),intent(out) :: dhdd  ! change in height per diameter [m/cm]
    
    h = p1*(1.0_r8 - exp(p2*d**p3))
    !h = h_max - h_max (exp(a*d**b))
    !f(x) = -h_max*exp(g(x))
    !g(x) = a*d**b
    !d/dx f(g(x) = f'(g(x))*g'(x) = -a1*exp(a2*d**a3) * a3*a2*d**(a3-1)
    !
    dhdd = -p1*exp(p2*d**p3) * p3*p2*d**(p3-1.0_r8)
    return
    
  end subroutine d2h_poorter2006

  ! ===========================================================================

  subroutine d2h_2pwr(d,p1,p2,dbh_hmax,h,dhdd)

    ! =========================================================================
    ! "d2h_2pwr"
    ! "d to height via 2 parameter power function"
    ! where height h is related to diameter by a linear relationship on the log
    ! transform where log(a) is the intercept and b is the slope parameter.
    !
    ! log(h) = log(a) + b*log(d)
    ! h      = exp(log(a)) * exp(log(d))**b
    ! h      = a*d**b
    !
    ! This functional form is used often in temperate allometries
    ! Therefore, no base reference is cited.  Although, the reader is pointed
    ! towards Dietze et al. 2008, King 1991, Ducey 2012 and many others for
    ! reasonable parameters.  Note that this subroutine is intended only for
    ! trees well below their maximum height, ie initialization.
    !
    ! =========================================================================
    ! From King et al. 1990 at BCI for saplings
    ! log(d) = a + b*log(h)
    ! d = exp(a) * h**b
    ! h = (1/exp(a)) * d**(1/b)
    ! h = p1*d**p2  where p1 = 1/exp(a)  p2 = 1/b
    ! d = (h/p1)**(1/p2)
    ! For T. tuberculata (canopy tree) a = -0.0704, b = 0.67
    ! =========================================================================

    ! args
    ! =========================================================================
    ! d: diameter at breast height
    ! p1: the intercept parameter 
    !                       (however exponential of the fitted log trans)
    ! p2: the slope parameter
    ! return:
    ! h: total tree height [m]
    ! =========================================================================

    
    real(r8),intent(in)     :: d        ! plant diameter [cm]
    real(r8),intent(in)     :: p1       ! parameter a
    real(r8),intent(in)     :: p2       ! parameter b
    real(r8),intent(in)     :: dbh_maxh ! dbh where max height occurs [cm]
    real(r8),intent(out)    :: h        ! plant height [m]
    real(r8),intent(out)    :: dhdd     ! change in height per diameter [m/cm]

    real(r8) :: dbh0,fl

    real(r8),parameter :: ddbh = 0.1_r8 ! 1-mm 
    real(r8),parameter :: k    = 0.25_r8  ! sharpness coef for logistic

    ! For the non derivative solution, if the tree is large and
    ! close to any cap that is imposed, then we need to perform a
    ! step-integration because the asymptotic function makes the indefinite
    ! integral incredibly messy. Thus we use an Euler step, yes ugly.
    
    if (d>0.5_r8*dbh_maxh) then
       dbh0=0.5_r8*dbh_maxh
       h = p1*dbh0**p2
       do while(dbh0<d)
          fl = 1.0_r8/(1.0_r8+exp(-k*(dbh0-dbh_maxh)))
          dhdd = (p2*p1)*dbh0**(p2-1.0_r8)*(1.0_r8-fl)
          h    = h+ddbh*dhdd
          dbh0 = dbh0+ddbh
       end do
       !    display("A request for a height calculation near hmax with chave")
       !    display("allometry required an explicit euler integration")
       !    display("this is innefficient, and was not thought to had been")
       !    display("necessary")
    else
       h = p1*d**p2
    end if
    
    ! The diameter at maximum height
    !    dbh_maxh = (hmax/p1)**(1/p2)
    ! Logistic function
    fl = 1.0_r8/(1.0_r8+exp(-k*(d-dbh_maxh)))
    ! Derivative of logistic function wrt d
    dhdd = (p2*p1)*d**(p2-1.0_r8)*(1.0_r8-fl)

    return
 end subroutine d2h_2pwr

  ! ============================================================================

  subroutine d2h_obrien(d,p1,p2,dbh_hmax,h,dhdd)
    
    real(r8),intent(in)    :: d        ! plant diameter [cm]
    real(r8),intent(in)    :: p1       ! parameter a
    real(r8),intent(in)    :: p2       ! parameter b
    real(r8),intent(in)    :: dbh_hmax ! dbh where max height occurs [cm]
    real(r8),intent(out)   :: h        ! plant height [m]
    real(r8),intent(out)   :: dhdd     ! change in height per diameter [m/cm]
    
    !p1 = 0.64
    !p2 = 0.37
    if(d>=dbh_hmax) then
       h    = 10.0_r8**(log10(dbh_hmax)*p1+p2)
       dhdd = 0.0_r8
    else
       h    = 10.0_r8**(log10(d)*p1+p2)
       dhdd = p1*10.0_r8**p2*d**(p1-1.0_r8)
    end if
    
    return
    
  end subroutine d2h_obrien

  ! ===========================================================================
  
  subroutine d2h_martcano(d,p1,p2,p3,h,dhdd)
     
     ! =========================================================================
     ! "d2h_martcano"
     ! "d to height via 3 parameter Michaelis-Menten following work at BCI
     ! by Martinez-Cano et al. 2016
     ! 
     ! h = (a*d**b)/(c+d**b)
     !
     ! h' = [(a*d**b)'(c+d**b) - (c+d**b)'(a*d**b)]/(c+d**b)**2
     ! dhdd = h' = [(ba*d**(b-1))(c+d**b) - (b*d**(b-1))(a*d**b)]/(c+d**b)**2
     !
     ! args
     ! =========================================================================
     ! d: diameter at breast height
     ! h: total tree height [m]
     ! =========================================================================
     
     real(r8),intent(in)  :: d     ! plant diameter [cm]   
     real(r8),intent(in)  :: p1       ! parameter a
     real(r8),intent(in)  :: p2       ! parameter b 
     real(r8),intent(in)  :: p3       ! parameter c
     real(r8),intent(out) :: h     ! plant height [m]
     real(r8),intent(out) :: dhdd  ! change in height per diameter [m/cm]
     
     
     h = (p1*d**p2)/(p3+d**p2)
     
     dhdd = ((p2*p1*d**(p2-1._r8))*(p3+d**p2) - &
             (p2*d**(p2-1._r8))*(p1*d**p2)        )/ &
             (p3+d**p2)**2._r8
       
     return
  end subroutine d2h_martcano
  

  ! =========================================================================
  ! Diameter 2 above-ground biomass
  ! =========================================================================
  
  subroutine dh2bag_chave2014(d,h,ipft,bag,dbagdd)

    ! =========================================================================
    ! This function calculates tree structural biomass from tree diameter,
    ! height and wood density.
    !
    ! Chave et al. Improved allometric models to estimate the abovegroud
    ! biomass of tropical trees.  Global Change Biology. V20, p3177-3190. 2015.
    !
    ! Input arguments:
    ! d: Diameter at breast height [cm]
    ! rho:  wood specific gravity (dry wood mass per green volume)
    ! height: total tree height [m]
    ! a1: structural biomass allometry parameter 1 (intercept)
    ! a2: structural biomass allometry parameter 2 (slope)
    ! Output:
    ! bag:   Total above ground biomass [kgC]
    !
    ! =========================================================================

    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: h       ! plant height [m]
    integer(li),intent(in)       :: ipft    ! PFT index
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out) :: dbagdd  ! change in agb per diameter [kgC/cm]

    real(r8) :: hj,dhdd
    real(r8) :: dbagdd1,dbagdd2,dbagdd3

    associate( &
         d2bag1 => EDecophyscon%d2bag1(ipft), &
         d2bag2 => EDecophyscon%d2bag2(ipft), &
         wood_density => EDecophyscon%wood_density(ipft), &
         c2b => EDecophyscon%c2b(ipft) )

      bag   = (d2bag1 * (wood_density*d**2.0_r8*h)**d2bag2)/c2b
      call h_allom(d,ipft,hj,dhdd)

!      call d2h_chave2014(d,ipft,hj,dhdd) ! This hj should ~= h
!                                         ! but perhaps not as precise

!      fl = 1.0_r8/(1.0_r8+exp(-k*(d-dbh_maxh)))
!      ae = p1-eclim
!      dhpdd = exp(ae)*( p3*2.0_r8*d**(p2-1.0_r8)*log(d)* &
!           exp(p3*log(d)**2.0_r8) + p2*d**(p2-1.0_r8)* &
!           exp(p3*log(d)**2.0_r8) )
!      dhdd = dhpdd*(1.0_r8-fl)

      dbagdd1  = (d2bag1*wood_density**d2bag2)/c2b
      dbagdd2  = d2bag2*d**(2.0_r8*d2bag2)*h**(d2bag2-1.0_r8)*dhdd
      dbagdd3  = h**d2bag2*2.0_r8*d2bag2*d**(2.0_r8*d2bag2-1.0_r8)
      dbagdd   = dbagdd1*(dbagdd2 + dbagdd3)

      return
    end associate
  end subroutine dh2bag_chave2014

  subroutine d2bag_2pwr(d,ipft,bag,dbagdd)

    ! =========================================================================
    ! This function calculates tree above ground biomass according to 2
    ! parameter power functions. (slope and intercepts of a log-log
    ! diameter-agb fit:
    !
    ! These relationships are typical for temperate/boreal plants in North
    ! America.  Parameters are available from Chojnacky 2014 and Jenkins 2003
    !
    ! Note that we are using an effective diameter here, as Chojnacky 2014
    ! and Jenkins use diameter at the base of the plant for "woodland" species
    ! The diameters should be converted prior to this routine if drc.
    !
    ! Input arguments:
    ! diam: effective diameter (d or drc) in cm
    ! FOR SPECIES THAT EXPECT DCM, THIS NEEDS TO BE PRE-CALCULATED!!!!
    ! Output:
    ! agb:   Total above ground biomass [kgC]
    !
    ! =========================================================================
    ! Aceraceae, Betulaceae, Fagaceae and Salicaceae comprised nearly
    ! three-fourths of the hardwood species (Table 3)
    !
    ! Fabaceae and Juglandaceae had specific gravities .0.60 and were
    ! combined, as were Hippocastanaceae and Tilaceae with specific gravities
    ! near 0.30. The remaining 9 families, which included mostly species with
    ! specific gravity 0.45–0.55, were initially grouped to construct a general
    ! hardwood taxon for those families having few published biomass equa-
    ! tions however, 3 warranted separation, leaving 6 families for the general
    ! taxon.
    ! bag = exp(b0 + b1*ln(diameter))/c2b
    ! =========================================================================

    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    integer(li),intent(in)       :: ipft    ! PFT index
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out) :: dbagdd  ! change in agb per diameter [kgC/cm]
    
    associate( &
         d2bag1 => EDecophyscon%d2bag1(ipft), &
         d2bag2 => EDecophyscon%d2bag2(ipft), &
         c2b    => EDecophyscon%c2b(ipft))
    
      !max_dbh = EDecophyscon%maxdbh(ipft)
      !if(diam>1.10*max_dbh)
      !    display("-----------------------------------------------------")
      !    display("Tree diameter is 10! larger than diameter where height")
      !    display("hits maximum.  However, you specified an AGB allometry")
      !    display("that does not assume capping. Please consider ")
      !    display("re-evaluating your allometric assumptions, growth")
      !    display("formulations or maximum height")
      !    display("------------------------------------------------------")
      !end
    
      bag    = (d2bag1 * d**d2bag2)/c2b
      dbagdd = (d2bag2*d2bag1*d**(d2bag2-1.0_r8))/c2b

      return
    end associate
    
  end subroutine d2bag_2pwr
  
  
  subroutine dh2bag_salda(d,h,ipft,bag,dbagdd)

    ! In the interest of reducing the number of model parameters, and since
    ! Saldarriaga 1988 seems as though it is being deprecated, we will use
    ! hard-wired parameter for dh2bag_salda, which would had required 4
    ! variable parameters.

    
    real(r8),intent(in)  :: d       ! plant diameter [cm]
    real(r8),intent(in)  :: h       ! plant height [m]
    integer(li),intent(in)       :: ipft    ! PFT index
    real(r8),intent(out) :: bag     ! plant height [m]
    real(r8),intent(out) :: dbagdd  ! change in agb per diameter [kgC/cm]

    real(r8) :: term1,term2,term3,hj,dhdd

    real(r8),parameter :: d2bag1       = 0.06896_r8
    real(r8),parameter :: d2bag2       = 0.572_r8
    real(r8),parameter :: d2bag3       = 1.94_r8
    real(r8),parameter :: d2bag4       = 0.931_r8

    associate( &
         c2b          => EDecophyscon%c2b(ipft), &
         wood_density => EDecophyscon%wood_density(ipft))
      
      bag = d2bag1*(h**d2bag2)*(d**d2bag3)*(wood_density**d2bag4)/c2b

      ! bag     = a1 * h**a2 * d**a3 * r**a4
      ! dbag/dd = a1*r**a4 * d/dd (h**a2*d**a3)
      ! dbag/dd = a1*r**a4 * [ d**a3 *d/dd(h**a2) + h**a2*d/dd(d**a3) ]
      ! dbag/dd = a1*r**a4 * [ d**a3 * a2*h**(a2-1)dh/dd + h**a2*a3*d**(a3-1)]
      
      term1 = d2bag1*(wood_density**d2bag4)/c2b
      term2 = (h**d2bag2)*d2bag3*d**(d2bag3-1.0_r8)
      
      call h_allom(d,ipft,hj,dhdd)
      term3 = d2bag2*h**(d2bag2-1)*(d**d2bag3)*dhdd
      dbagdd   = term1*(term2+term3)
      return
    end associate

  end subroutine dh2bag_salda

  ! ============================================================================
  ! height to diameter conversions
  ! Note that these conversions will only back-calculate the actual diameter
  ! for plants that have not started to experience height capping or an
  ! asymptote.  In these cases they may be called effective diameter.
  ! ============================================================================
  
  subroutine h2d_chave2014(h,ipft,de,ddedh)

    
    real(r8),intent(in)  :: h       ! plant height [m]
    integer(li),intent(in)       :: ipft    ! PFT index
    real(r8),intent(out) :: de      ! effective plant diameter [cm]
    real(r8),intent(out) :: ddedh   ! effective change in d per height [cm/m]

    real(r8) :: ar, eroot, dbh1,dhpdd

    associate( &
         d2h1_ad => EDecophyscon%d2h1_ad(ipft), &      ! (alias)
         d2h2_ad => EDecophyscon%d2h2_ad(ipft), &
         d2h3_ad => EDecophyscon%d2h3_ad(ipft), &
         eclim => EDecophyscon%eclim(ipft))

      ar    = d2h1_ad-eclim
      eroot = (-d2h2_ad + sqrt(d2h2_ad**2.0_r8 + 4.0_r8*log(h*exp(-ar))*d2h3_ad)) & 
           /(2.0_r8*d2h3_ad)

      de = exp(eroot)

      ! Invert the derivative at d? without asymtote
      dhpdd = exp(ar)*( d2h3_ad*2.0_r8*de**(d2h2_ad-1.0_r8)*log(de)* &
           exp(d2h3_ad*log(de)**2) + d2h2_ad*de**(d2h2_ad-1.0_r8)* &
           exp(d2h3_ad*log(de)**2.0_r8) )

      ddedh = 1.0_r8/dhpdd

      !    term1 = exp(-d2h2_ad/(2*d2h3_ad))
      !    term2 = exp(d2h2_ad**2/(4*d2h3_ad**2))
      !    term3 = exp(-ar/d2h3_ad)
      !    term4 = h**(1/d2h3_ad-1.0_r8)/(d2h3_ad)
      !    d   = term1*term2*term3*term4
      return
    end associate
  end subroutine h2d_chave2014



subroutine h2d_poorter2006(h,ipft,d,dddh)
  
  ! -------------------------------------------------------------------------
  ! Note that the height to diameter conversion in poorter is only necessary
  ! when initializing saplings.  In other methods, height to diameter is
  ! useful when defining the dbh at which point to asymptote, as maximum
  ! height is the user set parameter.  This function should not need to set a
  ! dbh_max parameter for instance, but it may end up doing so anyway, even
  ! if it is not used, do to poor filtering.  The poorter et al. d2h and h2d
  ! functions are already asymptotic, and the parameter governing maximum
  ! height is the d2h1_ad parameter.  Note as dbh gets very large, the
  ! exponential goes to zero, and the maximum height approaches d2h1_ad.
  ! However, the asymptote has a much different shape than the logistic, so
  ! therefore in the Poorter et al functions, we do not set d2h1_ad == h_max.
  ! That being said, if an h_max that is greater than d2h1_ad is passed to this
  ! function, it will return a complex number. During parameter
  ! initialization, a check will be placed that forces:
  ! h_max = d2h1_ad*0.98
  ! -------------------------------------------------------------------------

  
  real(r8),intent(in)  :: h      ! plant height [m]
  integer(li),intent(in)       :: ipft   ! PFT index
  real(r8),intent(out) :: d      ! plant diameter [cm]
  real(r8),intent(out) :: dddh   ! change in d per height [cm/m]
  
  associate( &
       d2h1_ad => EDecophyscon%d2h1_ad(ipft), &  !(alias)
       d2h2_ad => EDecophyscon%d2h2_ad(ipft), &  !(alias)
       d2h3_ad => EDecophyscon%d2h3_ad(ipft))  !(alias)
  
    ! -------------------------------------------------------------------------
    ! h = a1*(1 - exp(a2*d**a3))
    ! h = a1 - a1*exp(a2*d**a3)
    ! a1-h = a1*exp(a2*d**a3)
    ! (a1-h)/a1 = exp(a2*d**a3)
    ! log(1-h/a1) = a2*d**a3
    ! [log(1-h/a1)/a2]**(1/a3) = d
    !
    ! derivative dd/dh
    ! dd/dh = [log((a1-h)/a1)/a2]**(1/a3)'
    !       = (1/a3)*[log((a1-h)/a1)/a2]**(1/a3-1)* [(log(a1-h)-log(a1))/a2]'
    !       = (1/a3)*[log((a1-h)/a1)/a2]**(1/a3-1) * (1/(a2*(h-a1))
    ! dd/dh = -((log(1-h/a1)/a2)**(1/a3-1))/(a2*a3*(a1-h))
    ! -------------------------------------------------------------------------

    d = (log(1.0_r8-h/d2h1_ad)/d2h2_ad)**(1.0_r8/d2h3_ad)
    dddh = -((log(1-h/d2h1_ad)/d2h2_ad)**(1.0_r8/d2h3_ad-1.0_r8))/ &
         (d2h2_ad*d2h3_ad*(d2h1_ad-h))

    return
  end associate
end subroutine h2d_poorter2006


subroutine h2d_2pwr(h,ipft,d,dddh)

  
  real(r8),intent(in)  :: h      ! plant height [m]
  integer(li),intent(in)       :: ipft   ! PFT index
  real(r8),intent(out) :: d      ! plant diameter [cm]
  real(r8),intent(out) :: dddh   ! change in d per height [cm/m]
  
  associate( &
       d2h1_ad => EDecophyscon%d2h1_ad(ipft), &
       d2h2_ad => EDecophyscon%d2h2_ad(ipft))

    !h = a1*d**a2
    d = (h/d2h1_ad)**(1.0_r8/d2h2_ad)
    !    d = (1/a1)**(1/a2)*h**(1/a2)
    dddh = (1.0_r8/d2h2_ad)*(1.0_r8/d2h1_ad)**(1.0_r8/d2h2_ad) &
         *h**(1.0_r8/d2h2_ad-1.0_r8)

    return
  end associate
end subroutine h2d_2pwr

subroutine h2d_obrien(h,ipft,d,dddh)

  
  real(r8),intent(in)    :: h      ! plant height [m]
  integer(li),intent(in) :: ipft   ! PFT index
  real(r8),intent(out)   :: d      ! plant diameter [cm]
  real(r8),intent(out)   :: dddh   ! change in d per height [cm/m]
  
  real(r8) :: h_sap, h_adult, dddh_ad, dddh_sap

  associate( &
       d2h1_ad  => EDecophyscon%d2h1_ad(ipft), &
       d2h2_ad  => EDecophyscon%d2h2_ad(ipft), &
       h_max    => EDecophyscon%h_max(ipft), &
       d_sap    => EDecophyscon%d_sap(ipft), &
       d_adult  => EDecophyscon%d_adult(ipft), &
       d2h1_sap => EDecophyscon%d2h1_sap(ipft), &  !1.0729
       d2h2_sap => EDecophyscon%d2h2_sap(ipft))  !1.4925

    ! =========================================================================
    ! From King et al. 1990 at BCI for saplings
    ! log(d) = a + b*log(h)
    ! d = exp(a) * h**b
    ! h = (1/exp(a)) * d**(1/b)
    ! h = d2h1_ad*d**d2h2_ad  where d2h1_ad = 1/exp(a)  d2h2_ad = 1/b
    ! d = (h/d2h1_ad)**(1/d2h2_ad)
    ! For T. tuberculata (canopy tree) a1 = -0.0704, b1 = 0.67
    ! =========================================================================

    h_sap   = d2h1_sap * d_sap**d2h2_sap
    h_adult = 10.0_r8**(log10(d_adult) * d2h1_ad + d2h2_ad)

    if (h<h_sap) then
       d    = (h/d2h1_sap)**(1.0_r8/d2h2_sap)
       dddh = (1.0_r8/d2h2_sap)*(h/d2h1_sap)**(1.0_r8/d2h2_sap-1.0_r8)
    elseif (h>=h_sap .and. h<h_adult) then
    
       dddh_sap = (1.0_r8/d2h2_sap)*(h_sap/d2h1_sap)**(1.0_r8/d2h2_sap-1.0_r8)
       dddh_ad  = 1.0_r8/(d2h1_ad*10.0_r8**d2h2_ad*d_adult**(d2h1_ad-1.0_r8))
       call cspline(h_sap,h_adult,d_sap,d_adult,dddh_sap,dddh_ad,h,d,dddh)
    elseif (h<h_max) then
       d    = 10.0_r8**((log10(h)-d2h2_ad)/d2h1_ad)
       dddh = 1.0_r8/(d2h1_ad*10.0_r8**d2h2_ad*d**(d2h1_ad-1.0_r8))
    else
       d    = 10.0_r8**((log10(h_max)-d2h2_ad)/d2h1_ad)
       dddh = 1.0d20  ! Something super high, because it should be infinite
    end if
    return
  end associate
end subroutine h2d_obrien

  ! ============================================================================

  subroutine h2d_martcano(h,ipft,d,dddh)
     
     ! =========================================================================
     ! "d2h_martcano"
     ! "d to height via 3 parameter Michaelis-Menten following work at BCI
     ! by Martinez-Cano et al. 2016
     ! 
     ! h = (a*d**b)/(c+d**b)
     !
     ! d = [(h*c)/(a-h)]**(1/b)
     ! d = [(h*c)**(1/b)] / [(a-h)**(1/b)]
     ! d' = {[(h*c)**(1/b)]' [(a-h)**(1/b)] -  [(a-h)**(1/b)]'[(h*c)**(1/b)]} /
     !       [(a-h)**(1/b)]**2
     ! dddh = d' = {[(1/b)(h*c)**(1/b-1)] [(a-h)**(1/b)] -  
     !              [(1/b)(a-h)**(1/b-1)] [(h*c)**(1/b)]} /
     !             [(a-h)**(1/b)]**2
     ! 
     ! =========================================================================

     real(r8),intent(in)  :: h     ! plant height [m]
     integer(li),intent(in) :: ipft  ! PFT index
     real(r8),intent(out)   :: d     ! plant diameter [cm]
     real(r8),intent(out)   :: dddh  ! change in diameter per height [cm/m]

     associate( &
           a  => EDecophyscon%d2h1_ad(ipft), &
           b  => EDecophyscon%d2h2_ad(ipft), &
           c  => EDecophyscon%d2h3_ad(ipft))

       d = ((h*c)/(a-h))**(1._r8/b)

       dddh =  (((1._r8/b)*(h*c)**(1._r8/b-1._r8))*((a-h)**(1._r8/b)) - & 
                ((1._r8/b)*(a-h)**(1._r8/b-1._r8))* ((h*c)**(1._r8/b)) ) / &
                ((a-h)**(1._r8/b))**2._r8

       
     end associate
  end subroutine h2d_martcano

  ! ===========================================================================

subroutine cspline(x1,x2,y1,y2,dydx1,dydx2,x,y,dydx)

  ! ============================================================================
  ! This subroutine performs a cubic spline interpolation between known
  ! endpoints.  The endpoints have known coordinats and slopes
  ! ============================================================================

  ! Arguments

  real(r8),intent(in) :: x1     ! Lower endpoint independent
  real(r8),intent(in) :: x2     ! Upper endpoint independent
  real(r8),intent(in) :: y1     ! Lower endpoint dependent
  real(r8),intent(in) :: y2     ! Upper endpoint dependent
  real(r8),intent(in) :: dydx1  ! Lower endpoint slope
  real(r8),intent(in) :: dydx2  ! Upper endpoint slope
  real(r8),intent(in) :: x      ! Independent
  real(r8),intent(out) :: y     ! Dependent
  real(r8),intent(out) :: dydx  ! Slope
  
  ! Temps
  real(r8) :: t
  real(r8) :: a
  real(r8) :: b

  t = (x-x1)/(x2-x1)
  a = dydx1*(x2-x1) - (y2-y1)
  b = -dydx2*(x2-x1) + (y2-y1)

  y    = (1.0_r8-t)*y1 + t*y2 + t*(1.0_r8-t)*(a*(1.0_r8-t) + b*t)
  dydx = (y2-y1)/(x2-x1) + (1.0_r8-2.0_r8*t)*(a*(1.0_r8-t)+b*t)/(x2-x1) + t*(1.0_r8-t)*(b-a)/(x2-x1)

end subroutine cspline

end module EDAllomMod
