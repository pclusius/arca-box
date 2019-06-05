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
  INTEGER   :: MODE  = 0
  ! 0 = use values that are read in from column "col", possibly modifying by a factor or a constant
  ! 1 = Use NORMALD to create function in LINEAR mode
  ! 2 = Use NORMALD to create function in LOGARITMIC mode
  integer   :: col = -1     ! Column for input in whatever file the value will be

  real(dp)  :: multi = 1d0  ! Multiplication factor in MODE0
  real(dp)  :: shift = 0d0  ! Constant to be added in MODE0
  real(dp)  :: min = 0d0    ! Minimum value for the parametrized concentration OR constant value if max <= min
  real(dp)  :: max = 1d5    ! Peak value
  real(dp)  :: sig = 1d0    ! Standard deviation for the Gaussian=sig of the bell curve
  real(dp)  :: mju = 12d0   ! Time of peak value
  real(dp)  :: fv  = 0d0    ! Angular frequency [hours] of modifying sine function
  real(dp)  :: ph  = 0d0    ! Angular frequency [hours] of modifying sine function
  real(dp)  :: am  = 1d0    ! Amplitude of modificaion
  CHARACTER(16) :: NAME = 'NONAME'! Human readable name for modified variable
end type input_mod

type timetype
  real(dp)      :: SIM_TIME_H     = 1d0
  real(dp)      :: SIM_TIME_S     = 3600.d0
  real(dp)      :: dt             = 10.0d0
  real(dp)      :: sec            = 0
  real(dp)      :: min            = 0
  real(dp)      :: hrs            = 0
  real(dp)      :: day            = 0
  real(dp)      :: dt_chem        = 10.0d0
  real(dp)      :: dt_aero        = 10.0d0
  real(dp)      :: PRINT_INTERVAL = 15d0  *60d0
  real(dp)      :: FSAVE_INTERVAL = 5d0   *60d0
  integer       :: ind_netcdf     = 1
  integer       :: JD             = 0
  character(8)  :: hms            = "00:00:00"
  logical       :: printnow       = .true.
  logical       :: savenow        = .true.
  logical       :: PRINTACDC      = .false.
end type timetype

type(timetype)  :: MODELTIME
type(input_mod), allocatable :: MODS(:) ! THIS VECTOR HOLDS ALL INPUT AND MODIFICATION PARAMETERS
CHARACTER(len=12), allocatable :: UNITS(:) ! THIS VECTOR HOLDS ALL INPUT UNITS

!type for number of column and rows
type nrowcol
  integer :: rows,cols
end type nrowcol


! ------------------------------------------------------------
! PROCEDURES
interface operator(+)
  module procedure ADD
end interface operator(+)

interface operator(.mod.)
  module procedure MOD_CONC
end interface operator(.mod.)

CONTAINS

! =================================================================================================
! Timetype update function.
! .................................................................................................
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
  IF (MODULO(nint(ADD%sec*100), NINT(ADD%PRINT_INTERVAL*100)) == 0) THEN
    ADD%printnow = .true.
  ELSE
    ADD%printnow = .false.
  END IF
  IF (MODULO(nint(ADD%sec*100), NINT(ADD%FSAVE_INTERVAL*100)) == 0) THEN
    ADD%savenow = .true.
    ADD%ind_netcdf = ADD%ind_netcdf + 1
  ELSE
    ADD%savenow = .false.
  END IF
end function ADD

REAL(dp) FUNCTION MOD_CONC(c, MODS)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  REAL(dp), INTENT(in)        :: c
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
