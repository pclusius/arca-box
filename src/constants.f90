MODULE constants
USE SECOND_PRECISION,  ONLY : dp, sp
IMPLICIT NONE
PUBLIC

real(dp), parameter   :: Na = 6.022140857d23  ! 1/mol Avogadro constant
real(dp), parameter   :: R  = 8.3144598       ! [J/K/mol] Universal gas constant
real(dp), parameter   :: kb = R/Na            ! [J/K] Boltzmann constant
real(dp), parameter   :: pi = ACOS(-1d0)      ! pi
real(dp), parameter   :: K0 = 273.15d0        ! [K] Zero degree celcius in K
integer(sp), parameter:: min_s = 60           ! [s] seconds in minute
integer(sp), parameter:: hour_s  = 3600       ! [s] seconds in hour
integer(sp), parameter:: day_s = 24*hour_s
real(dp), parameter   :: um3_to_m3 = (1D-6)**3 ! used for  vol_concentration
REAL(dp), PARAMETER   :: Mair = 28.96D-3          ! Mean molecular weight of air (kg)

! Saturation vapour pressure of water in Pa
REAL, PARAMETER       :: a0 = 6.107799961,     & ! Parameters to calculate the saturation vapour pressure for water
                         a1 = 4.436518524E-1,  &
                         a2 = 1.428945805E-2,  &
                         a3 = 2.650648471E-4,  &
                         a4 = 3.031240396E-6,  &
                         a5 = 2.034080948E-8,  &
                         a6 = 6.136820929E-11

! ----------------------------------------------------------------
! USER-DEFINED TYPES

! Container for input, name, unit and parameters for creating simulated datapoints

type input_mod
  ! Mode of operation:
  INTEGER   :: MODE  = 0
  ! 0 = use values that are read in from column "col", possibly modifying by a factor or a constant
  ! 1 = Use NORMALD to create function in LINEAR mode
  ! 2 = Use NORMALD to create function in LOGARITMIC mode
  integer   :: col = -1     ! Column for input in whatever file the value will be. if -1, only modifiers are used

  real(dp)  :: multi = 1d0  ! Multiplication factor in MODE0
  real(dp)  :: shift = 0d0  ! Constant to be added in MODE0
  real(dp)  :: min = 0d0    ! Minimum value for the parametrized concentration OR constant value if max <= min
  real(dp)  :: max = 1d5    ! Peak value
  real(dp)  :: sig = 1d0    ! Standard deviation for the Gaussian=sig of the bell curve
  real(dp)  :: mju = 12d0   ! Time of peak value
  real(dp)  :: fv  = 0d0    ! Angular frequency [hours] of modifying sine function
  real(dp)  :: ph  = 0d0    ! Angular frequency [hours] of modifying sine function
  real(dp)  :: am  = 1d0    ! Amplitude of modificaion
  CHARACTER(5)  :: UNIT = '#'      ! Unit for the given number. CASE INSENSITIVE
  CHARACTER(16) :: NAME = 'NONAME' ! Human readable name for modified variable

  ! UNITS FOR INPUT (case insensitive):
  ! #      = number concentration in 1/cm3. DEFAULT ASSUMPTION
  ! 'ppm'  = parts per million (1/1e6)
  ! 'ppb'  = parts per billion (1/1e9)
  ! 'ppt'  = parts per trillion (1/1e12)
  ! 'ppq'  = parts per quadrillion (1/1e15)
  ! 'Pa'   = Pascals
  ! 'hPa'  = Hectopascals
  ! 'kPa'  = Kilopascals
  ! 'mbar' = millibars
  ! 'atm'  = Atmosphere (101325 Pascals)
  ! 'K'    = Kelvin
  ! 'C'    = Celsius

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
  logical       :: printnow       = .false.
  logical       :: savenow        = .true.
  logical       :: PRINTACDC      = .false.
end type timetype

!type for number of column and rows
type nrowcol
  integer :: rows,cols
end type nrowcol

!error type which is used to optimize computation speed at given simulatin precision
TYPE error_type
  LOGICAL :: error_state = .false.  !there is no error at the start
  INTEGER :: error_process  !process where the error occurs
  CHARACTER(150) :: error_specification  !specification on error type (e.g."particle conc" during coagulation)
END TYPE error_type

! Type for storing the (currently only the extra) particles
TYPE inert_particles
  REAL(dp), ALLOCATABLE :: binseries(:,:) ! the time series and bins
  REAL(dp), ALLOCATABLE :: options(:)     ! Additional information, like molar mass, density etc.
  REAL(dp), ALLOCATABLE :: time(:)        ! Time vector of the particles
  REAL(dp), ALLOCATABLE :: sections(:)    ! diameters for the centers of the sections
  CHARACTER(20) :: name  ! Name for the stuff
END TYPE inert_particles


! ------------------------------------------------------------
! PROCEDURES
interface operator(+)
  module procedure ADD
end interface operator(+)

interface operator(.mod.)
  module procedure MOD_CONC
end interface operator(.mod.)
! ------------------------------------------------------------

type(timetype)  :: MODELTIME
type(input_mod), allocatable  :: MODS(:) ! THIS VECTOR HOLDS ALL INPUT AND MODIFICATION PARAMETERS
TYPE(input_mod)               :: ORIGINAL_TEMP(2), ORIGINAL_press(2)
REAL(dp)                      :: J_ACDC_NH3 = 0d0
REAL(dp)                      :: J_ACDC_DMA = 0d0
REAL(dp)                      :: J_NH3_BY_IONS(3) = 0d0
REAL(dp)                      :: acdc_cluster_diam = 2.17d-9

REAL(dp) :: C_AIR_NOW, RESOLVED_BASE, RESOLVED_J
CONTAINS

! =================================================================================================
! Timetype update function. When timestep is added, this will update all other time-related
! variables accordingly
! .................................................................................................
PURE type(timetype) function ADD(time, sec)
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


! =================================================================================================
! Function to calculate concentrations based on the input units and modifyers
! .................................................................................................
PURE REAL(dp) FUNCTION MOD_CONC(c, MODS)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  REAL(dp), INTENT(in)        :: c

  if (MODS%MODE == 0) THEN
    MOD_CONC = c * MODS%multi
    MOD_CONC = MOD_CONC + MODS%shift

  ELSE
    MOD_CONC = NORMALD(MODS)
  END IF

  ! If the concentration is a mixing ratio, defined in MODS%UNIT, it will be converted to number concentration
  if (UCASE(TRIM(MODS%UNIT)) == 'PPM') THEN
    MOD_CONC = MOD_CONC * 1d-6 * C_AIR_NOW
  elseif (UCASE(TRIM(MODS%UNIT)) == 'PPB') THEN
    MOD_CONC = MOD_CONC * 1d-9 * C_AIR_NOW
  elseif (UCASE(TRIM(MODS%UNIT)) == 'PPT') THEN
    MOD_CONC = MOD_CONC * 1d-12 * C_AIR_NOW
  elseif (UCASE(TRIM(MODS%UNIT)) == 'PPQ') THEN
    MOD_CONC = MOD_CONC * 1d-15 * C_AIR_NOW
  elseif (UCASE(TRIM(MODS%UNIT)) == 'HPA') THEN
    MOD_CONC = MOD_CONC * 1d2
  elseif (UCASE(TRIM(MODS%UNIT)) == 'KPA') THEN
    MOD_CONC = MOD_CONC * 1000d0
  elseif (UCASE(TRIM(MODS%UNIT)) == 'ATM') THEN
    MOD_CONC = MOD_CONC * 1.01325d5
  elseif (UCASE(TRIM(MODS%UNIT)) == 'BAR') THEN
    MOD_CONC = MOD_CONC * 1.d5
  elseif (UCASE(TRIM(MODS%UNIT)) == 'MBAR') THEN
    MOD_CONC = MOD_CONC * 1.d2
  END if

END FUNCTION MOD_CONC

!==============================================================================
! Function to return y-value at [time] from a normal distribution function with
! standard deviation [sig], minimum value [min], maximum value [max]
! and time of maximum at [mju]. Relies completely on the correct formulation
! of MODS (input_mod)
!..............................................................................
PURE REAL(dp) FUNCTION NORMALD(MODS, timein)
  IMPLICIT NONE
  type(input_mod), INTENT(in) :: MODS
  type(timetype), OPTIONAL, INTENT(in)  :: timein
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



!==============================================================================
! Change character string to uppercase
!..............................................................................
PURE FUNCTION UCASE(word)
  IMPLICIT NONE
  character(*), INTENT(IN) :: word
  character(len=len(word)) :: UCASE
  integer :: i
  UCASE = word
  FORALL (i=1:len(word), ((ichar(word(i:i))>96) .and. (ichar(word(i:i))<123))) UCASE(i:i) = char(ichar(word(i:i))-32)
END FUNCTION UCASE

end MODULE constants
