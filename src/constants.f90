MODULE constants
USE SECOND_PRECISION,  ONLY : dp, ik
IMPLICIT NONE
PUBLIC

real(dp), parameter   :: Na = 6.022140857d23  ! 1/mol Avogadro constant
real(dp), parameter   :: R  = 8.3144598       ! [J/K/mol] Universal gas constant
real(dp), parameter   :: kb = R/Na            ! [J/K] Boltzmann constant
real(dp), parameter   :: pi = ACOS(-1d0)      ! pi
real(dp), parameter   :: K0 = 273.15          ! [K] Zero degree celcius in K
integer(ik), parameter:: min_s = 60           ! [s] seconds in minute
integer(ik), parameter:: hour_s  = 3600       ! [s] seconds in hour
integer(ik), parameter:: day_s = 24*hour_s

! ----------------------------------------------------------------
! USER-DEFINED TYPES

! Container for input parameters for creating simulated datapoints
! Uses gaussian function as generator
type parametered_input
  real(dp) :: sigma       ! Standard deviation for the Gaussian
  real(dp) :: min_c       ! Minimum value for the parametrized concentration OR constant value if max_c <= min_c
  real(dp) :: max_c       ! Peak value
  real(dp) :: peaktime    ! Time of peak value
  real(dp) :: omega       ! Angular frequency [hours] of modifying sine function
  real(dp) :: amplitude   ! Amplitude of modificaion
  LOGICAL  :: LOGSCALE    ! Is concetration scale logaritmically or linearily.
                          ! True for logaritmically scaled
end type parametered_input

real(dp), parameter :: SIM_TIME = 24.0d0
type timetype
  real(dp):: SIM_TIME_H     = SIM_TIME
  real(dp):: SIM_TIME_S     = SIM_TIME*3600.d0
  real(dp):: sec            = 0
  real(dp):: dt             = 10.0d0
  real(dp):: dt_chem        = 10.0d0
  real(dp):: dt_aero        = 10.0d0
  integer :: ind_netcdf     = 1
  integer :: tm             = 0
  integer :: day            = 0
  integer :: PRINT_INTERVAL = 15
  integer :: FSAVE_INTERVAL = 1
  character(8) :: hms       = "00:00:00"
  logical :: printnow       = .true.
  logical :: savenow        = .true.
  logical :: PRINTACDC      = .false.
end type timetype

! UHMA-wide parameters that need to be constants
integer,parameter ::  uhma_sections= 100                  !size sections!60
integer,parameter ::  uhma_vbs_bins= 699                  !615, & ! 598, & !785, & !802, &!785,& !1291, &    !When input_flag=1: number of organic species (volatility basis set, Donahue et al., 2008) 316 11! luxis file is 189 and pontus is 221 or 220 !
integer,parameter ::  uhma_cond=uhma_vbs_bins             ! number of condensable species (sulphuric acid + organics)
integer,parameter ::  uhma_noncond=0                      ! number of non-condensable species
integer,parameter ::  uhma_compo=uhma_cond+uhma_noncond   ! number of composition classes
integer,parameter ::  uhma_init_modes=2                   ! number of initial particle modes

! holds variables related to particles
type type_particles
  real(dp) :: n_limit,nuc_rate, cluster_vol
  real(dp),dimension(uhma_sections) :: n_conc,radius,rdry,rdry_orig,core,mass,gr, original_radiis
  real(dp),dimension(uhma_sections,uhma_compo) :: vol_conc ! um^3/m^3
  integer, dimension(uhma_sections) :: order
  CHARACTER(*) :: DESCRIPTION = "Type to store particle number concentrations and volumes, masses, etc. and bins"
end type type_particles

! holds variables related to vapors and ambient conditions
type type_ambient
  real(dp)                        :: temp,pres,rh,vap_limit,gf_org,nh3_mix,ELVOC_Nucl,oh_conc,boundaryheight,nuc_coeff, so2_mix
  real(dp),dimension(uhma_compo)  :: density,molarmass
  real(dp),dimension(uhma_cond)   :: c_sat,surf_ten,diff_vol,alpha,molecvol,molecmass,vap_conc,sink,n_crit, vap_conc_min
  real(dp), dimension(uhma_cond)  :: parameter_a, parameter_b
  integer                         :: sulfuric_acid_index
  integer, dimension(uhma_compo)  :: vapor_type, condensing_type, index_in_chemistry
  character(len=60), dimension(uhma_cond) :: vapor_names
end type type_ambient

! holds misc. variables and options
type type_options
  real(dp) :: meas_interval,latitude,longitude,day_of_year,nuc_number
  integer :: year, month, day
  integer, dimension(3) :: nuc_array
  character (len=2) :: dist_approach
  character (len=3) :: nuc_approach,solver
  logical :: nucleation,condensation,coagulation,dry_deposition,snow_scavenge,equilibrate,cluster_divide, vapor_chemistry, vapor_loss,BLH_dilution
  logical :: quasistationary, raoult
end type type_options

! ------------------------------------------------------------
! PROCEDURES
interface operator(+)
  module procedure ADD
end interface operator(+)

type(type_options) :: aerosol_options

CONTAINS

type(timetype) function add(time, sec)
  implicit none
  type(timetype), intent(in) :: time
  real(dp), intent(in)       :: sec
  add     = time
  add%sec  = time%sec + sec
  add%tm  = int(add%sec/60)
  add%day = int(add%sec/3600/24)
  write(add%hms, '(i2.2, ":" i2.2, ":" i2.2)') int(add%sec+0.5d0)/3600, &
        int(MODULO(int(add%sec+0.5d0),3600)/60), MODULO(MODULO(int(add%sec+0.5d0),3600), 60)
  IF (MODULO(int(add%sec+0.5d0), 60*add%PRINT_INTERVAL) == 0) THEN
    add%printnow = .true.
  ELSE
    add%printnow = .false.
  END IF

  IF (MODULO(int(add%sec+0.5d0), 60*add%FSAVE_INTERVAL) == 0) THEN
    add%savenow = .true.
    add%ind_netcdf = add%ind_netcdf + 1
  ELSE
    add%savenow = .false.
  END IF

end function add



end MODULE constants
