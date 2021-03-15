MODULE constants
USE SECOND_PRECISION,  ONLY : dp, sp
IMPLICIT NONE
PUBLIC

real(dp), parameter    :: Na = 6.022140857d23    ! 1/mol Avogadro constant
real(dp), parameter    :: Rg = 8.3144598         ! [J/K/mol] Universal gas constant
real(dp), parameter    :: kb = Rg/Na             ! [J/K] Boltzmann constant
real(dp), parameter    :: pi = ACOS(-1d0)        ! pi
real(dp), parameter    :: K0 = 273.15d0          ! [K] Zero degree celcius in K
integer(sp), parameter :: min_s = 60             ! [s] seconds in minute
integer(sp), parameter :: hour_s  = 3600         ! [s] seconds in hour
integer(sp), parameter :: day_s = 24*hour_s      ! Seconds in a day
real(dp), parameter    :: um3_to_m3 = (1D-6)**3  ! used for  vol_concentration
REAL(dp), PARAMETER    :: Mair = 28.9647D-3      ! Mean molecular weight of air (kg)
REAL(dp), PARAMETER    :: g_0 = 9.80665D0        ! Gravitational acceleration
REAL(dp), PARAMETER    :: Diff_H2SO4_0 = 0.09D-4 ! H2SO4 diffusivity at 0 RH (why no temperature and pressure dependence?)
REAL(dp), PARAMETER    :: Keq1 = 0.13D0          ! H2SO4 diffusivity RH dependence parameters, Hanson & Eisele
REAL(dp), PARAMETER    :: Keq2 = 0.016D0         ! H2SO4 diffusivity RH dependence parameters, Hanson & Eisele
INTEGER, PARAMETER     :: dint = selected_int_kind(16)

! Saturation vapour pressure of water in Pa
REAL, PARAMETER        :: a0 = 6.107799961,     & ! Parameters to calculate the saturation vapour pressure for water
                          a1 = 4.436518524E-1,  &
                          a2 = 1.428945805E-2,  &
                          a3 = 2.650648471E-4,  &
                          a4 = 3.031240396E-6,  &
                          a5 = 2.034080948E-8,  &
                          a6 = 6.136820929E-11


! ======================================================================================================================
! USER-DEFINED TYPES
! ======================================================================================================================

! ===========================================================================
! Type for input, name, unit and parameters for creating simulated datapoints
type input_mod

    INTEGER   :: MODE  = 0 ! Mode of operation:
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
    real(dp)  :: am  = 1d0    ! Amplitude of modification
    CHARACTER(5)  :: UNIT = '-'      ! Unit for the given number. CASE INSENSITIVE
    CHARACTER(16) :: TIED = ''       ! variable is coupled with this variable
    CHARACTER(16) :: NAME = 'NONAME' ! Human readable name for modified variable
    LOGICAL :: ISPROVIDED = .false.  ! Human readable name for modified variable

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
    ! 'bar'  = bars
    ! 'atm'  = Atmosphere (101325 Pascals)
    ! 'K'    = Kelvin
    ! 'C'    = Celsius

end type input_mod

type timetype
    real(dp)      :: SIM_TIME_H     = 1d0       ! [h]
    real(dp)      :: SIM_TIME_S     = 3600.d0   ! [sec]
    real(dp)      :: dt             = 10.0d0    ! [sec]
    real(dp)      :: sec            = 0
    real(dp)      :: min            = 0
    real(dp)      :: hrs            = 0
    real(dp)      :: day            = 0
    real(dp)      :: dt_chm         = 10.0d0    ! [sec]
    real(dp)      :: dt_aer         = 10.0d0    ! [sec]
    real(dp)      :: PRINT_INTERVAL = 15d0*60d0 ! [sec]
    real(dp)      :: FSAVE_INTERVAL = 5d0 *60d0 ! [sec]
    integer       :: ind_netcdf     = 1         ! index for output file
    integer       :: prevPrint_i    = 0         ! index for output file
    integer       :: prevSave_i     = 0         ! index for output file
    integer       :: JD             = 0         ! Julian day of the year, calculated if date is provided
    character(8)  :: hms            = "00:00:00"! Pretty print of time
    logical       :: printnow       = .true.    ! this flag is true when prints are wanted
    logical       :: savenow        = .true.    ! this flag is true when saves are wanted
    logical       :: PRINTACDC      = .false.   ! this flag is true when cluster information is printed
end type timetype


!===============================================================
! error type which is used to optimize computation speed at given simulation precision
!===============================================================
TYPE error_type
    LOGICAL         :: err  = .false.  ! There is no error at the start
    INTEGER         :: proc = 0        ! Process where the error occurs
    INTEGER         :: cch  = 1        ! Process where the error occurs
    INTEGER         :: coa  = 2        ! Process where the error occurs
    INTEGER         :: dep  = 3        ! Process where the error occurs
    CHARACTER(15)   :: pr_name(3) = ['Cond & Chem',&
                                     'Coagulation',&
                                     'Deposition ']
    CHARACTER(150)  :: err_text         ! Specification on error type (e.g."particle conc" during coagulation)
    LOGICAL         :: increase(3) = .false.  ! .true. if speed can be increased for cch,coa or dep
    LOGICAL         :: in_turn(4)  = .true.  ! .true. if process is in turn, 4th is for any process == .true.
END TYPE error_type


!===============================================================
! Type for storing the particles and particle losses
!===============================================================
TYPE particle_grid
    REAL(dp), ALLOCATABLE :: conc_matrix(:,:)   ! the time series of bins
    REAL(dp), ALLOCATABLE :: options(:)         ! Additional information, like molar mass, density etc.
    REAL(dp), ALLOCATABLE :: time(:)            ! Time vector of the particles
    REAL(dp), ALLOCATABLE :: sections(:)        ! diameters for the centers of the sections
    REAL(dp), ALLOCATABLE :: conc_modelbins(:)  ! current concentration fitted to model bins
    CHARACTER(20) :: name                       ! Name for the stuff
END TYPE particle_grid



!===============================================================
! type describing the particle size distribution and composition
!===============================================================
type PSD
    INTEGER  :: PSD_style  !sets type of particle size distribution representation
    INTEGER  :: nr_bins    !initial number of bins, particle sizes for any method
    REAL(dp) :: dp_range(2)    !lower and upper limit of diameter for simulation

    ! FULL STATIONARY METHOD
    REAL(dp), ALLOCATABLE :: diameter_fs(:)          ! particle diameter [m] (nr_bins)
    REAL(dp), ALLOCATABLE :: dp_dry_fs(:)            ! dry particle diameter [m] (nr_bins)
    REAL(dp), ALLOCATABLE :: volume_fs(:)            ! particle volume   [m³ / particle] (nr_bins)
    REAL(dp), ALLOCATABLE :: density_fs(:)           ! particle density  [kg * m⁻³] (n_cond_tot)
    REAL(dp), ALLOCATABLE :: particle_density_fs(:)  ! particle density [kg * m⁻³] (nr_bins)
    REAL(dp), ALLOCATABLE :: particle_mass_fs(:)     ! particle mass [kg] (nr_bins)
    REAL(dp), ALLOCATABLE :: composition_fs(:,:)     ! mass of all species in the particle phase [kg/particle] (nr_bins,n_cond_tot)
    REAL(dp), ALLOCATABLE :: conc_fs(:)              ! particle concentration in each size bin [m⁻³] (nr_bins)

    ! MOVING AVERAGE, FIXED GRID
    REAL(dp), ALLOCATABLE :: diameter_ma(:)          ! particle diameter [m] (nr_bins)
    REAL(dp), ALLOCATABLE :: dp_dry_ma(:)            ! dry particle diameter [m] (nr_bins)
    REAL(dp), ALLOCATABLE :: volume_ma(:)            ! particle volume   [m³ / particle] (nr_bins)
    REAL(dp), ALLOCATABLE :: density_ma(:)           ! particle density  [kg * m⁻³] (n_cond_tot)
    REAL(dp), ALLOCATABLE :: particle_density_ma(:)  ! particle density [kg * m⁻³] (nr_bins)
    REAL(dp), ALLOCATABLE :: particle_mass_ma(:)     ! particle mass [kg] (nr_bins)
    REAL(dp), ALLOCATABLE :: composition_ma(:,:)     ! mass of all species in the particle phase [kg/particle] (nr_bins,n_cond_tot)
    REAL(dp), ALLOCATABLE :: conc_ma(:)              ! particle concentration in each size bin [m⁻³] (nr_bins)
    REAL(dp), ALLOCATABLE :: grid_ma(:)              ! grid that contains the bin borders [m⁻³] (nr_bins+1)

END TYPE PSD



! This datatype contains all parameters for input vapours
type :: vapour_ambient
    integer                         :: n_condorg            ! number of organics, all type 1
    integer                         :: n_condtot            ! number of all vapours, all type 1+2 (including acids, generic etc.)
    integer                         :: ind_H2SO4            ! index for H2SO4
    integer                         :: ind_GENERIC          ! index of GENERIC composition
    integer, allocatable            :: cond_type(:)         ! 1 = Organic vapour 2 = Acid (that's it for the moment)
    real(dp), allocatable           :: molar_mass(:)        ! molar mass [kg]
    real(dp), allocatable           :: psat_a(:), psat_b(:) ! Parameters to calculate sat. vapour press.: log10(Psat) = a-b/T Psat=[atm]
    real(dp), allocatable           :: alpha(:)             ! accomodation coefficient i.e. "sticking coef."= 1.0 []
    real(dp), allocatable           :: density(:)           ! bulk density of liquid [kg/m³]
    real(dp), allocatable           :: surf_tension(:)      ! Surface tension [N/m]
    real(dp), allocatable           :: diff(:)              ! mass diffusivity
    real(dp), allocatable           :: c_speed(:)           ! thermal speed
    real(dp), allocatable           :: molec_dia(:)         ! molecule diameter [m]
    real(dp), allocatable           :: molec_mass(:)        ! molecule mass [kg]
    real(dp), allocatable           :: molec_volume(:)      ! molecule volume [m³]
    real(dp), allocatable           :: diff_vol(:)          ! diffusion volume [m³]
    real(dp), allocatable           :: diff_dia(:)          ! diffusion diameter [m³]
    real(dp), allocatable           :: wet_dia(:)           ! wet molecule diameter [m]
    real(dp), allocatable           :: wet_mass(:)          ! wet molecule mass [kg]
    real(dp), allocatable           :: c_sat(:)             ! saturation concentration 1/m3
    real(dp), allocatable           :: mfractions(:)        ! dimension(tot_spec) mole fractions
    character(len=25), allocatable  :: vapour_names(:)      ! Must match those in chemistry, come from "Vapours file"
end type vapour_ambient


! ------------------------------------------------------------
! PROCEDURES
interface operator(.mod.)
    module procedure MOD_CONC
end interface operator(.mod.)


! ------------------------------------------------------------

type(timetype)                  :: GTIME                      ! Model time, available in all modules and usually not explicitly sent in subroutines etc.
type(input_mod), allocatable    :: MODS(:)                    ! THIS VECTOR HOLDS ALL INPUT AND MODIFICATION PARAMETERS
type(vapour_ambient)            :: VAPOUR_PROP                ! Stores all vapour properties
REAL(dp) :: GC_AIR_NOW, GTEMPK, GPRES, GRH, GCS               ! Global Air concentration, Temperature, Pressure, Relative humidity and H2SO4 condensation sink
REAL(dp)                        :: J_ACDC_NH3_M3 = 0d0        ! Formation rate of NH3 calculated in ACDC [1/s/m³]
REAL(dp)                        :: J_ACDC_DMA_M3 = 0d0        ! Formation rate of DMA calculated in ACDC [1/s/m³]
REAL(dp)                        :: J_TOTAL_M3    = 0d0        ! Formation rates from all methods summed [1/s/m³]
REAL(dp)                        :: clusteracid, clusterbase   ! stores informatin of cluster population
REAL(dp)                        :: dclusteracid, dclusterbase ! stores information of cluster population
REAL(dp)                        :: J_NH3_BY_IONS(3) = 0d0     ! Stores charge dependent formation rates [1/s/m³]
REAL(dp)                        :: acdc_cluster_diam = 1.4d-9 ! diameter of the largest clusters, updated in ACDC
REAL(dp)                        :: RESOLVED_BASE, RESOLVED_J  ! stores output of SOLVE_BASES
REAL(dp)                        :: RESOLVED_J_FACTR           ! stores output of SOLVE_BASES
Logical                         :: NO_NEGATIVE_CONCENTRATIONS = .true.

! speed_up: factor for increasing integration time step for individual prosesses
! (1): chemistry & Condensation & ACDC (2): Coagulation; (3): Deposition
TYPE(error_type)                :: PRC
INTEGER(dint)                   :: speed_up(size(PRC%pr_name,1)) = 1


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
        if (sec < 0) THEN
            IF ((ADD%PRINT_INTERVAL*(ADD%prevPrint_i))-ADD%sec  >= sec   ) THEN
                ADD%printnow = .false.
                ADD%prevPrint_i = ADD%prevPrint_i - 1
            END IF
            if ((ADD%FSAVE_INTERVAL*(ADD%prevSave_i))-ADD%sec  >= sec   ) THEN
                ADD%savenow = .false.
                ADD%prevSave_i = ADD%prevSave_i - 1
            END IF
        END IF

        ADD%sec = time%sec + sec

    ELSE
        ADD%sec = time%sec + time%dt
    END IF
    ADD%min = ADD%sec/60d0
    ADD%hrs = ADD%sec/3600d0
    ADD%day = ADD%sec/3600d0/24d0
    write(ADD%hms, '(i2.2, ":" i2.2, ":" i2.2)') nint(ADD%sec)/3600, &
    int(MODULO(nint(ADD%sec),3600)/60d0), MODULO(MODULO(nint(ADD%sec),3600), 60)
    IF ((ADD%printnow .eqv. .false.) .and. (ADD%sec >= (ADD%PRINT_INTERVAL*(1+ADD%prevPrint_i)))) THEN
        ADD%printnow = .true.
        ADD%prevPrint_i = INT(ADD%sec)/INT(ADD%PRINT_INTERVAL)
    ELSE
        ADD%printnow = .false.
    END IF
    IF ((ADD%savenow .eqv. .false.) .and. (ADD%sec >= (ADD%FSAVE_INTERVAL*(1+ADD%prevSave_i)))) THEN
        ADD%savenow = .true.
        ADD%prevSave_i =  INT(ADD%sec)/INT(ADD%FSAVE_INTERVAL)
    ELSE
        ADD%savenow = .false.
    END IF

end function ADD


! =================================================================================================
! Function to calculate concentrations based on the input units and modifyers
! .................................................................................................
PURE REAL(dp) FUNCTION MOD_CONC(c, MOD)
    IMPLICIT NONE
    type(input_mod), INTENT(in) :: MOD
    REAL(dp), INTENT(in)        :: c

    if (MOD%MODE < 1) THEN
        MOD_CONC = c * MOD%multi
        MOD_CONC = MOD_CONC + MOD%shift
    ELSE
        MOD_CONC = NORMALD(MOD)
    END IF

    ! If the concentration is a mixing ratio, defined in MOD%UNIT, it will be converted to number concentration
    if (UCASE(TRIM(MOD%UNIT)) == 'PPM') THEN
        MOD_CONC = MOD_CONC * 1d-6 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPB') THEN
        MOD_CONC = MOD_CONC * 1d-9 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPT') THEN
        MOD_CONC = MOD_CONC * 1d-12 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPQ') THEN
        MOD_CONC = MOD_CONC * 1d-15 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'HPA') THEN
        MOD_CONC = MOD_CONC * 1d2
    elseif (UCASE(TRIM(MOD%UNIT)) == 'KPA') THEN
        MOD_CONC = MOD_CONC * 1000d0
    elseif (UCASE(TRIM(MOD%UNIT)) == 'ATM') THEN
        MOD_CONC = MOD_CONC * 1.01325d5
    elseif (UCASE(TRIM(MOD%UNIT)) == 'BAR') THEN
        MOD_CONC = MOD_CONC * 1.d5
    elseif (UCASE(TRIM(MOD%UNIT)) == 'MBAR') THEN
        MOD_CONC = MOD_CONC * 1.d2
    END if

END FUNCTION MOD_CONC


PURE REAL(dp) FUNCTION UCONV(c,MOD)
    IMPLICIT NONE
    type(input_mod), INTENT(in) :: MOD
    REAL(dp), INTENT(in)        :: c

    ! If the concentration is a mixing ratio, defined in MODS%UNIT, it will be converted to number concentration
    if (UCASE(TRIM(MOD%UNIT)) == 'PPM') THEN
        UCONV = c * 1d-6 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPB') THEN
        UCONV = c * 1d-9 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPT') THEN
        UCONV = c * 1d-12 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'PPQ') THEN
        UCONV = c * 1d-15 * GC_AIR_NOW
    elseif (UCASE(TRIM(MOD%UNIT)) == 'HPA') THEN
        UCONV = c * 1d2
    elseif (UCASE(TRIM(MOD%UNIT)) == 'KPA') THEN
        UCONV = c * 1000d0
    elseif (UCASE(TRIM(MOD%UNIT)) == 'ATM') THEN
        UCONV = c * 1.01325d5
    elseif (UCASE(TRIM(MOD%UNIT)) == 'BAR') THEN
        UCONV = c * 1.d5
    elseif (UCASE(TRIM(MOD%UNIT)) == 'MBAR') THEN
        UCONV = c * 1.d2
    elseif (UCASE(TRIM(MOD%UNIT)) == '#') THEN
        UCONV = c
    END if

END FUNCTION UCONV

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
        time = GTIME
    END IF

    D       = MODS%mju + sin((GTIME%hrs-MODS%mju)*MODS%fv)*MODS%am + MODS%ph
    f       = 1d0/SQRT(2d0*pi*MODS%sig**2)
    NORMALD = f * EXP(- (GTIME%hrs-D)**2/(2d0*MODS%sig**2))
    ! Check if minimum is same or more than maximum and if so, use constant concentration (minumum)
    IF (ABS(MODS%max-MODS%min) <1e-12) THEN
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


SUBROUTINE SET_ERROR(proc, msg)
    implicit none
    integer, intent(in) :: proc
    character(len=*), INTENT(IN) :: msg

    PRC%err = .true.
    PRC%proc = proc
    PRC%err_text = msg

END SUBROUTINE SET_ERROR

end MODULE constants
