MODULE constants
USE SECOND_PRECISION,  ONLY : dp, ik
IMPLICIT NONE
PUBLIC

real(dp), parameter   :: Na = 6.022140857d23  ! 1/mol Avogadro constant
real(dp), parameter   :: R  = 8.3144598       ! [J/K/mol] Universal gas constant
real(dp), parameter   :: kb = R/Na            ! [J/K] Boltzmann constant
real(dp), parameter   :: pi = ACOS(-1d0)      ! pi
real(dp), parameter   :: K0 = 273.15d0        ! [K] Zero degree celcius in K
integer(ik), parameter:: min_s = 60           ! [s] seconds in minute
integer(ik), parameter:: hour_s  = 3600       ! [s] seconds in hour
integer(ik), parameter:: day_s = 24*hour_s

! ----------------------------------------------------------------
! USER-DEFINED TYPES

! Container for input parameters for creating simulated datapoints
! Uses gaussian function as generator. Use ParameterTweaker.py for more complex functions

type input_mod
  ! Mode of operation:
  ! 0 = use values that are read in, poissibly modifying by a factor or a constant
  ! 1 = Use NORMALD to create function in LINEAR mode
  ! 2 = Use NORMALD to create function in LOGARITMIC mode
  INTEGER   :: MODE  = 0
  CHARACTER(16) :: NAME = '-'! Human readable name for modified value
  real(dp)  :: min = 0d0     ! Minimum value for the parametrized concentration OR constant value if max <= min
  real(dp)  :: max = 1d5     ! Peak value
  real(dp)  :: sig = 1d0     ! Standard deviation for the Gaussian=sig of the bell curve
  real(dp)  :: mju = 12d0    ! Time of peak value
  real(dp)  :: fv  = 0d0    ! Angular frequency [hours] of modIFying sine function
  real(dp)  :: ph  = 0d0    ! Angular frequency [hours] of modIFying sine function
  real(dp)  :: am  = 1d0    ! Amplitude of modIFicaion
  real(dp)  :: multi = 1d0   ! Multiplication factor in MODE0
  real(dp)  :: shift = 0d0   ! Constant to be added in MODE0
end type input_mod

real(dp), parameter :: SIM_TIME = 24.0d0

type timetype
  real(dp)      :: SIM_TIME_H     = SIM_TIME
  real(dp)      :: SIM_TIME_S     = SIM_TIME*3600.d0
  real(dp)      :: dt             = 10.0d0
  real(dp)      :: sec            = 0
  real(dp)      :: min            = 0
  real(dp)      :: hrs            = 0
  real(dp)      :: day            = 0
  real(dp)      :: dt_chem        = 10.0d0
  real(dp)      :: dt_aero        = 10.0d0
  integer       :: ind_netcdf     = 1
  integer       :: JD             = 0
  integer       :: PRINT_INTERVAL = 15
  integer       :: FSAVE_INTERVAL = 5
  character(8)  :: hms            = "00:00:00"
  logical       :: printnow       = .true.
  logical       :: savenow        = .true.
  logical       :: PRINTACDC      = .false.
end type timetype

type(timetype)  :: MODELTIME

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
end type type_particles

! holds variables related to vapors and ambient conditions
type type_ambient
  real(dp)                        :: temp,pres,rh,vap_limit,gf_org,nh3_mix,ELVOC_Nucl,oh_conc,boundaryheight,nuc_coeff, so2_mix
  real(dp),dimension(uhma_compo)  :: density,molarmass
  real(dp),dimension(uhma_cond)   :: c_sat,surf_ten,dIFf_vol,alpha,molecvol,molecmass,vap_conc,sink,n_crit, vap_conc_min
  real(dp), dimension(uhma_cond)  :: parameter_a, parameter_b
  integer                         :: sulfuric_acid_index
  integer, dimension(uhma_compo)  :: vapor_type, condensing_type, index_in_chemistry
  character(len=60), dimension(uhma_cond) :: vapor_names
end type type_ambient

! holds misc. variables and options
type type_options
  ! real(dp) :: meas_interval ! Commented out, this is in timetype
  real(dp)              :: latitude,longitude, day_of_year, nuc_number
  real(dp)              :: cons_multipliers(5)
  real(dp)              :: cons_shIFters(5)
  integer               :: year, month, day
  integer, dimension(3) :: nuc_array
  character (len=2)     :: dist_approach
  character (len=3)     :: nuc_approach,solver
  logical               :: nucleation,condensation,coagulation,dry_deposition,snow_scavenge,equilibrate,cluster_divide, vapor_chemistry, vapor_loss,BLH_dilution
  logical               :: quasistationary, raoult
end type type_options

! ------------------------------------------------------------
! PROCEDURES
interface operator(+)
  module procedure ADD
  module procedure PLUS
end interface operator(+)

interface operator(*)
  module procedure MULTIPLICATION
end interface operator(*)

interface operator(.mod.)
  module procedure MOD_CONC
end interface operator(.mod.)

type(type_options) :: aerosol_options

CONTAINS

type(timetype) function ADD(time, sec)
  implicit none
  type(timetype), intent(in)            :: time
  real(dp),       intent(in), optional  :: sec
  ADD = time
  IF (present(sec)) THEN
    ADD%sec = time%sec + sec
  ELSE
    ADD%sec = time%sec + time%dt
  END IF
  ADD%min = ADD%sec/60d0
  ADD%hrs = ADD%sec/3600d0
  ADD%day = ADD%sec/3600d0/24d0
  write(ADD%hms, '(i2.2, ":" i2.2, ":" i2.2)') nint(ADD%sec)/3600, &
    int(MODULO(nint(ADD%sec),3600)/60), MODULO(MODULO(nint(ADD%sec),3600), 60)
  IF (MODULO(nint(ADD%sec), 60*ADD%PRINT_INTERVAL) == 0) THEN
    ADD%printnow = .true.
  ELSE
    ADD%printnow = .false.
  END IF
  IF (MODULO(nint(ADD%sec), 60*ADD%FSAVE_INTERVAL) == 0) THEN
    ADD%savenow = .true.
    ADD%ind_netcdf = ADD%ind_netcdf + 1
  ELSE
    ADD%savenow = .false.
  END IF
end function ADD

REAL(dp) FUNCTION PLUS(c, MODS)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  REAL(dp), INTENT(in)        :: c
  PLUS = c + MODS%min
END FUNCTION PLUS

REAL(dp) FUNCTION MULTIPLICATION(c, MODS)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  REAL(dp), INTENT(in)        :: c
  MULTIPLICATION = c * MODS%multi
END FUNCTION MULTIPLICATION

REAL(dp) FUNCTION MOD_CONC(c, MODS)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  REAL(dp), INTENT(in)                :: c
  if (MODS%MODE == 0) THEN
    MOD_CONC = c * MODS%multi
    MOD_CONC = MOD_CONC + MODS%shift
  ELSE
    MOD_CONC = NORMALD(MODS)
  END IF
END FUNCTION MOD_CONC

!==============================================================================
! Function to return y-value at [time] from a normal distribution function with
! standard deviation [sig], minimum value [min], maximum value [max]
! and time of maximum at [mju]. Relies completely on the correct formulation
! of MODS (input_mod)
!..............................................................................
REAL(dp) FUNCTION NORMALD(MODS, timein)
  IMPLICIT NONE
  type(input_mod)   :: MODS
  type(timetype), OPTIONAL  :: timein
  type(timetype)            :: time
  REAL(dp) :: f, D
  if (PRESENT(timein)) THEN
    time = timein
  ELSE
    time = MODELTIME
  END IF

  D = MODS%mju + sin((MODELTIME%hrs-MODS%mju)*MODS%fv)*MODS%am + MODS%ph
  NORMALD = 1d0/SQRT(2d0*pi*MODS%sig**2) * EXP(- (MODELTIME%hrs-D)**2/(2d0*MODS%sig**2))
  f = 1d0/SQRT(2d0*pi*MODS%sig**2)
  ! Check if minimum is same or more than maximum and if so, use constant concentration (minumum)
  IF (ABS(MODS%max-MODS%min) <1e-9) THEN
    NORMALD = MODS%min
  ELSE
    if (MODS%MODE == 2) THEN
      f = (LOG10(MODS%max-MODS%min+1))/f
      NORMALD = 10**(NORMALD*f)-1 + MODS%min
    ELSEIF  (MODS%MODE == 1) THEN
      f = (MODS%max-MODS%min)/f
      NORMALD = NORMALD*f + MODS%min
    END IF
  END IF
end FUNCTION NORMALD

end MODULE constants
