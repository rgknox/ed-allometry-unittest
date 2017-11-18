!
! gfortran -shared -fPIC -g -o EDAllomUnitWrap.o EDAllomUnitWrap.F90
!
!

module FatesConstantsMod

   use iso_c_binding, only: fates_r8  => c_double
   use iso_c_binding, only: fates_int => c_int

   real(fates_r8), parameter :: kg_per_Megag = 1000.0
   real(fates_r8), parameter :: cm2_per_m2   = 10000.0
   real(fates_r8), parameter :: g_per_kg     = 1000.0

end module FatesConstantsMod


module shr_log_mod

   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int

   contains

   function shr_log_errMsg(source, line) result(ans)
      character(kind=c_char,len=*), intent(in) :: source
      integer(c_int), intent(in) :: line
      character(kind=c_char,len=128) :: ans

      ans = "source: " // trim(source) // " line: "
   end function shr_log_errMsg
   
end module shr_log_mod


module FatesGlobals

   contains

   integer function fates_log()
      fates_log = -1
   end function fates_log
   
   subroutine fates_endrun(msg) 

      implicit none
      character(len=*), intent(in) :: msg    ! string to be printed
      
      stop

   end subroutine fates_endrun

end module FatesGlobals



module EDPftvarcon
   
   use iso_c_binding, only : r8 => c_double
   use iso_c_binding, only : i4 => c_int
   use iso_c_binding, only : c_char

   integer,parameter :: SHR_KIND_CS = 80                     ! short char

   type, public :: EDPftvarcon_inst_type

     real(r8), pointer :: allom_hmode(:)
     real(r8), pointer :: allom_amode(:)
     real(r8), pointer :: allom_lmode(:)
     real(r8), pointer :: allom_smode(:)
     real(r8), pointer :: allom_cmode(:)
     real(r8), pointer :: allom_fmode(:)
     
     real(r8), pointer :: allom_d2h1(:)
     real(r8), pointer :: allom_d2h2(:)
     real(r8), pointer :: allom_d2h3(:)
     
     real(r8), pointer :: allom_dbh_maxheight(:)

     real(r8), pointer :: allom_agb1(:)
     real(r8), pointer :: allom_agb2(:)
     real(r8), pointer :: allom_agb3(:)
     real(r8), pointer :: allom_agb4(:)

     real(r8), pointer :: allom_d2bl1(:)
     real(r8), pointer :: allom_d2bl2(:)
     real(r8), pointer :: allom_d2bl3(:)

     real(r8), pointer :: wood_density(:)

     real(r8), pointer :: c2b(:)
     real(r8), pointer :: allom_latosa_int(:)
     real(r8), pointer :: allom_latosa_slp(:)
     real(r8), pointer :: slatop(:)

     real(r8), pointer :: allom_l2fr(:)
     real(r8), pointer :: allom_agb_frac(:)

     real(r8), pointer :: allom_blca_expnt_diff(:)
     real(r8), pointer :: allom_d2ca_coefficient_min(:)
     real(r8), pointer :: allom_d2ca_coefficient_max(:)


  end type EDPftvarcon_inst_type
 
  type pftptr_var
     real(r8), dimension(:), pointer :: var_rp
     integer(i4), dimension(:), pointer :: var_ip
     character(len=shr_kind_cs) :: var_name
     integer :: vtype
  end type pftptr_var

  type EDPftvarcon_ptr_type
     type(pftptr_var), allocatable :: var(:)
  end type EDPftvarcon_ptr_type
  

  type(EDPftvarcon_inst_type), public :: EDPftvarcon_inst ! ED ecophysiological constants structure
  type(EDPftvarcon_ptr_type),  public :: EDPftvarcon_ptr  ! Pointer structure for obj-oriented id
  
  integer :: numparm ! Number of different PFT parameters
  integer :: numpft

contains
  

   subroutine EDPftvarconPySet(ipft,rval,ival,name)

      implicit none
      ! Arguments
      integer(i4),intent(in) :: ipft
      character(kind=c_char,len=*), intent(in) :: name
      real(r8),intent(in) :: rval
      integer(i4),intent(in) :: ival
      ! Locals
      logical :: npfound
      integer :: ip
      integer :: namelen
      
      namelen = len(trim(name))

      print*,"F90: ARGS: ",trim(name)," IPFT: ",ipft," RVAL: ",rval," IVAL: ",ival

      ip=0
      npfound = .true.
      do ip=1,numparm

         if (trim(name) == trim(EDPftvarcon_ptr%var(ip)%var_name ) ) then
            print*,"F90: Found ",trim(name)," in lookup table"
            npfound = .false.
            if(EDPftvarcon_ptr%var(ip)%vtype == 1) then ! real
               EDPftvarcon_ptr%var(ip)%var_rp(ipft) = rval
            elseif(EDPftvarcon_ptr%var(ip)%vtype == 2) then ! integer
               EDPftvarcon_ptr%var(ip)%var_ip(ipft) = ival
            else
               print*,"F90: STRANGE TYPE"
               stop
            end if
         end if
      end do

      if(npfound)then
         print*,"F90: The parameter you loaded DNE: ",name(:)
         stop
      end if

      ! Performa a check to see if the target array is being filled

      if (trim(name) == 'wood_density' ) then
         if (EDPftvarcon_inst%wood_density(ipft) == rval) then
            print*,"F90: POINTER CHECK PASSES:",rval," = ",EDPftvarcon_inst%wood_density(ipft)
         else
            print*,"F90: POINTER CHECK FAILS:",rval," != ",EDPftvarcon_inst%wood_density(ipft)
            stop
         end if
      end if

      return
   end subroutine EDPftvarconPySet


  subroutine EDPftvarconAlloc(numpft_in)
    !

    ! !ARGUMENTS:
    integer(i4), intent(in) :: numpft_in
    ! LOCALS:
    integer                    :: iv   ! The parameter incrementer
    !------------------------------------------------------------------------

    numpft = numpft_in

    allocate( EDPftvarcon_ptr%var (100) ) ! Make this plenty large

    iv=0

    allocate( EDPftvarcon_inst%allom_dbh_maxheight   (1:numpft)); EDPftvarcon_inst%allom_dbh_maxheight (:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_dbh_maxheight"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_dbh_maxheight
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_hmode(1:numpft)); EDPftvarcon_inst%allom_hmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_hmode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_hmode
    EDPftvarcon_ptr%var(iv)%vtype    = 1
    
    allocate( EDPftvarcon_inst%allom_amode(1:numpft)); EDPftvarcon_inst%allom_amode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_amode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_amode
    EDPftvarcon_ptr%var(iv)%vtype    = 1
    
    allocate( EDPftvarcon_inst%allom_lmode(1:numpft)); EDPftvarcon_inst%allom_lmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_lmode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_lmode
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_smode(1:numpft)); EDPftvarcon_inst%allom_smode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_smode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_smode
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_cmode(1:numpft)); EDPftvarcon_inst%allom_cmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_cmode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_cmode
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_fmode(1:numpft)); EDPftvarcon_inst%allom_fmode(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_fmode"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_fmode
    EDPftvarcon_ptr%var(iv)%vtype    = 1
     
    allocate( EDPftvarcon_inst%allom_d2h1(1:numpft)); EDPftvarcon_inst%allom_d2h1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h1"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2h1
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2h2(1:numpft)); EDPftvarcon_inst%allom_d2h2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h2"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2h2
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2h3(1:numpft)); EDPftvarcon_inst%allom_d2h3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2h3"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2h3
    EDPftvarcon_ptr%var(iv)%vtype    = 1
    
    allocate( EDPftvarcon_inst%allom_agb1(1:numpft)); EDPftvarcon_inst%allom_agb1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb1"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_agb1
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_agb2(1:numpft)); EDPftvarcon_inst%allom_agb2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb2"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_agb2
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_agb3(1:numpft)); EDPftvarcon_inst%allom_agb3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb3"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_agb3
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_agb4(1:numpft)); EDPftvarcon_inst%allom_agb4(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb4"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_agb4
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2bl1(1:numpft)); EDPftvarcon_inst%allom_d2bl1(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl1"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2bl1
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2bl2(1:numpft)); EDPftvarcon_inst%allom_d2bl2(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl2"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2bl2
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2bl3(1:numpft)); EDPftvarcon_inst%allom_d2bl3(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2bl3"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2bl3
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%wood_density(1:numpft)); EDPftvarcon_inst%wood_density(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_wood_density"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%wood_density
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%c2b(1:numpft)); EDPftvarcon_inst%c2b(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_c2b"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%c2b
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_latosa_int(1:numpft)); EDPftvarcon_inst%allom_latosa_int(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_latosa_int"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_latosa_int
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_latosa_slp(1:numpft)); EDPftvarcon_inst%allom_latosa_slp(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_latosa_slp"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_latosa_slp
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%slatop(1:numpft)); EDPftvarcon_inst%slatop(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_slatop"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%slatop
    EDPftvarcon_ptr%var(iv)%vtype    = 1
    
    allocate( EDPftvarcon_inst%allom_l2fr(1:numpft)); EDPftvarcon_inst%allom_l2fr(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_l2fr"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_l2fr
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_agb_frac(1:numpft)); EDPftvarcon_inst%allom_agb_frac(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_agb_frac"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_agb_frac
    EDPftvarcon_ptr%var(iv)%vtype    = 1
    
    allocate( EDPftvarcon_inst%allom_blca_expnt_diff(1:numpft)); EDPftvarcon_inst%allom_blca_expnt_diff(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_blca_expnt_diff"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_blca_expnt_diff
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2ca_coefficient_min(1:numpft)); EDPftvarcon_inst%allom_d2ca_coefficient_min(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2ca_coefficient_min"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2ca_coefficient_min
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    allocate( EDPftvarcon_inst%allom_d2ca_coefficient_max(1:numpft)); EDPftvarcon_inst%allom_d2ca_coefficient_max(:) = nan
    iv = iv + 1
    EDPftvarcon_ptr%var(iv)%var_name = "fates_allom_d2ca_coefficient_max"
    EDPftvarcon_ptr%var(iv)%var_rp   => EDPftvarcon_inst%allom_d2ca_coefficient_max
    EDPftvarcon_ptr%var(iv)%vtype    = 1

    
    numparm = iv


    print*,"F90: ALLOCATED ",numparm," PARAMETERS, FOR ",numpft," PFTs"


    return
 end subroutine EDPftvarconAlloc

end module EDPftvarcon
