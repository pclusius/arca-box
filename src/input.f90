! ===================================================================================
! Atmospherically Relevant Chemistry and Aerosol box model
! Copyright (C) 2021  Multi-Scale ing group
! Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
! Contact information arca@helsinki.fi
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ===================================================================================


Module INPUT
! Module to read init file
use second_Monitor, ONLY: SPC_NAMES
USE second_Parameters, ONLY: NREACT
USE second_precision, ONLY: dp
use constants
USE auxillaries

Implicit none

INTEGER :: N_VARS ! This will store the number of variables in NAMES.DAT
INTEGER :: N_XTRS ! This will store the number of variables in AEMS.DAT
INTEGER :: LENV   ! This will store the number of named indices in this code

! INDICES TO PROPERLY COMBINE INPUT TO CORRECT VALUES IN THE MODEL
! ALL VARIABLES THAT ARE TIME DEPENDENT MUST BE HERE
! IF YOU ADD VARIABLES HERE, YOU NEED TO UPDATE:
! - this list
! - subroutine NAME_MODS_SORT_NAMED_INDICES
!------------------------------------------------------------------
INTEGER :: inm_TempK = 0
INTEGER :: inm_pres = 0
INTEGER :: inm_RH = 0
INTEGER :: inm_CS = 0
INTEGER :: inm_CS_NA = 0
INTEGER :: inm_swr = 0
INTEGER :: inm_IPR = 0
INTEGER :: inm_H2SO4 = 0
INTEGER :: inm_NH3 = 0
INTEGER :: inm_DMA = 0
INTEGER :: inm_SO2 = 0
INTEGER :: inm_NO = 0
INTEGER :: inm_NO2 = 0
INTEGER :: inm_CO = 0
INTEGER :: inm_H2 = 0
INTEGER :: inm_O3 = 0
INTEGER :: inm_JIN = 0
INTEGER :: N_CONCENTRATIONS = 0

INTEGER, ALLOCATABLE :: INDRELAY_CH(:)
INTEGER, ALLOCATABLE :: INDRELAY_TIED(:)
INTEGER, ALLOCATABLE :: index_cond(:)
! INTEGER, ALLOCATABLE :: index_inorg(:)

REAL(dp), allocatable, private :: INPUT_ENV(:,:)  ! will be of same shape as the files
REAL(dp), allocatable, private :: INPUT_MCM(:,:)  ! will be of same shape as the files
REAL(dp), allocatable :: timevec(:)     ! Whatever the times were for (currently) ALL measurements
REAL(dp), allocatable :: CONC_MAT(:,:)  ! will be of shape ( len(timevec) : N_VARS )
real(dp), allocatable :: par_data(:,:)
real(dp), allocatable :: GGR(:)
INTEGER               :: H2SO4_ind_in_chemistry = 0 ! Will be checked later and set to correct if H2SO4 is found
INTEGER               :: OH_ind_in_chemistry    = 0 ! Will be checked later and set to correct if H2SO4 is found
REAL(dp)              :: swr_spectrum(84) = 0d0     ! Vector for holding SWR spectral data, read in once from file or the interpolated values from swr_temporal_data
LOGICAL               :: swr_is_time_dependent = .false.
REAL(dp), ALLOCATABLE :: swr_times(:)               ! Vectors for holding SWR spectral data time stamps, used if SWR spectrum time dependent
REAL(dp), ALLOCATABLE :: swr_temporal_data(:,:)     ! Vectors for holding SWR spectral data, used if SWR spectrum time dependent
INTEGER, ALLOCATABLE  :: init_only_these(:)         ! Indices of the variables in MODS that will only be read in the first loop

! variable for storing init file name
CHARACTER(len=256) :: Fname_init ! init file names
CHARACTER(len=5)   :: gui        ! magic word for gui in use
LOGICAL            :: ingui = .False. ! True if program is invoked from gui, from command line
CHARACTER(len=14), parameter  :: namelists(21) = &
['NML_TIME      ','NML_Flag      ','NML_Path      ','NML_MISC      ','NML_VAP       ','NML_PARTICLE  ','NML_ENV       ',&
 'NML_MCM       ','NML_MODS      ','NML_PRECISION ','NML_CUSTOM    ','NML_ACDC      ','NML_NAMES     ','              ',&       ! to help troubleshooting
 '              ','              ','              ','              ','              ','              ','              ']       ! to help troubleshooting

! MAIN PATHS
! CHARACTER(len=256):: WORK_DIR   = ''
CHARACTER(len=256):: INOUT_DIR   = 'INOUT'
CHARACTER(len=256):: CASE_NAME  = 'DEFAULTCASE'
CHARACTER(len=90) :: RUN_NAME   = 'DEFAULTRUN'
NAMELIST /NML_Path/ INOUT_DIR, Case_name, RUN_NAME

! MODULES IN USE OPTIONS
Logical :: Chemistry_flag      = .false.
Logical :: Aerosol_flag        = .false.
Logical :: ACDC_solve_ss       = .false.
! Logical :: NUCLEATION          = .false.
Logical :: ORG_NUCL            = .false.
Logical :: ACDC                = .false.
Logical :: model_H2SO4         = .false.
Logical :: Condensation        = .false.
Logical :: Coagulation         = .false.
Logical :: Deposition          = .false.
Logical :: Chem_Deposition     = .false.
! Logical :: Extra_data          = .false.
Logical :: RESOLVE_BASE        = .false.
Logical :: PRINT_ACDC          = .false.
Logical :: OPTIMIZE_DT         = .false. ! Will be deprecated after 1.2
Logical :: USE_SPEED           = .false. ! Replacing OPTIMIZE_DT
Logical :: AFTER_CHEM_ON       = .false.
Logical :: AFTER_NUCL_ON       = .false.
CHARACTER(len=3) :: FILE_TIME_UNIT  = 'day'
CHARACTER(len=3) :: LOSSFILE_TIME_UNIT  = 'day'

NAMELIST /NML_Flag/ chemistry_flag, Aerosol_flag, ACDC_solve_ss, ACDC, & !NUCLEATION,
         Condensation, Coagulation, Deposition, Chem_Deposition, model_H2SO4, RESOLVE_BASE, &
         PRINT_ACDC, USE_SPEED,OPTIMIZE_DT, ORG_NUCL, AFTER_CHEM_ON, AFTER_NUCL_ON, FILE_TIME_UNIT,LOSSFILE_TIME_UNIT !,INIT_W_MODAL, Extra_data

! TIME OPTIONS
real(dp)  :: runtime = 1d0
real(dp)  :: FSAVE_INTERVAL = 300d0
real(dp)  :: PRINT_INTERVAL = 15*60d0
INTEGER   :: FSAVE_DIVISION = 0
real(dp)  :: dt = -1d0
CHARACTER(len=10)  :: DATE = '1800-01-01', NUMBER = ''
NAMELIST /NML_TIME/ runtime, dt, FSAVE_INTERVAL, PRINT_INTERVAL, FSAVE_DIVISION, DATE, NUMBER

! MODIFIER OPTIONS
! MODS is declared in CONSTANTS.f90, in order to be more widely available
NAMELIST /NML_MODS/ MODS

! ----------------------------------------------------------------------------------------------------------------------
! Particle related variables
INTEGER             :: PSD_MODE = 0
! PSD representation used
! 0 = Fully stationary
! 1 = Moving average
! 2 = Fully moving, not yet implemented
INTEGER             :: n_bins_par = 100    ! number of bins in the particle range
REAL(dp)            :: min_particle_diam = 1d-9 ! lower limit of particle range [m]
REAL(dp)            :: max_particle_diam = 2d-6 ! upper limit of particle range [m]

! DMPS INPUT
REAL(dp)            :: DMPS_read_in_time = 0d0

! 6) PSD input as discussed: real (time, total conc in first two columns, diameter in first row and concentration matrix) (nr_times + 1, nr_channels + 2)
! CHARACTER(len=256)  :: DMPS_dir
! CHARACTER(len=256)  :: extra_p_dir
CHARACTER(len=256)  :: DMPS_file

! 4) name list of species making up the particle phase (we think that this will only be a few species that are measured in the particle phase): integer
! 5) number of nonvolatile species considered: integer
! 7) Composition of the particles: real (time: first column, diameter: first row, species mass fraction matrix [1]) (nr_times + 1, nr_channels)
INTEGER             :: n_xpar_options = 3    ! number of options in the extra_particles file, not including path and name
! The list is needed for every species in the read in particle phase (see "4) name list of species ...") or we make one long list
! (nr_times + 1, nr_bins * particle phase species)
! Note: the number of channels does not have to be the number of size bins
CHARACTER(len=256)  :: extra_particles = '' ! file containing paths to extra particle sumfile
CHARACTER(len=500)  :: mmodal_input    = '' ! String defining the Modal PSD
INTEGER             :: mmodal_input_inuse  = -1 ! Switch toggling the Modal PSD
REAL(dp)            :: dmps_highband_lower_limit = 0d0    !for use_dmps_partial, read all dmps data above this diameter [m]
REAL(dp)            :: dmps_lowband_upper_limit = 0d0  !for use_dmps_partial, read all dmps data below this diameter [m]
logical             :: use_dmps = .false.
logical             :: use_dmps_partial = .false.
REAL(dp)            :: N_MODAL = -1d0   ! number of modes in multimodal intialization
real(dp)            :: dmps_interval = -1d0

NAMELIST /NML_PARTICLE/ PSD_MODE,n_bins_par,min_particle_diam,max_particle_diam, DMPS_file,extra_particles,& !DMPS_dir,extra_p_dir,
                        DMPS_read_in_time,dmps_highband_lower_limit, dmps_lowband_upper_limit,use_dmps,use_dmps_partial, &
                        mmodal_input, N_MODAL, mmodal_input_inuse, dmps_interval

REAL(dp), ALLOCATABLE            :: MMODES(:)       ! Modal PSD parameters
INTEGER,  ALLOCATABLE            :: MMODES_EMS(:)   ! Modal PSD emission parameters
INTEGER                          :: N_MODAL_EMS = 0 ! number of modes in multimodal intialization

type(particle_grid), ALLOCATABLE :: xtras(:)
! BG_PAR is here in case you want to use it Carlton and Lukas, in the end we remove either BG_PAR or par_data mmkay.
type(particle_grid) :: BG_PAR       ! Var to store the particle size distribution. This might become redundant
type(particle_grid) :: PAR_LOSSES   ! Var to store losses file, which could be either one row or a matrix, but has time
                                    ! and size dependency


! ENVIRONMENTAL INPUT
CHARACTER(len=256)  :: ENV_FILE = ''
CHARACTER(len=256)  :: LOSSES_FILE = ''
REAL(dp)            :: CONSTANT_PAR_LOSS_RATE = -1d0
! Chamber properties
REAL(dp)            :: CHAMBER_FLOOR_AREA     = 0d0
REAL(dp)            :: CHAMBER_HEIGHT = 0d0
REAL(dp)            :: CHAMBER_CIRCUMFENCE    = 0d0 ! Deprecated since v1.2, instead calculated assuming square floor
! Aerosol loss parametrisation
REAL(dp)            :: ustar = 5d-2                 ! [m/s] Friction velocity, affects particle wall losses
! vapour loss parametrisation
REAL(dp)            :: EDDYK = 5d-2                 ! [1/s] Coefficient of eddy diffusion - describes turbulence in the chamber
REAL(dp)            :: ALPHAWALL = 5d-5             ! [-] Wall loss accommodation coefficient, wall/component property, assumed constant
REAL(dp)            :: Cw_eqv = 40d-6               ! [mol/m³] Equivalent mass concentration of the wall, equilibrium wall vapour concentration

! Shortwave / actinic flux
CHARACTER(len=256)  :: spectrumfile = ''
Logical             :: SWR_IS_ACTINICFLUX  = .false.
REAL(dp)            :: SWR_IN_LOWER = 300d0  ! [nm] wavelength range of the Global SW Irradiation measurements
REAL(dp)            :: SWR_IN_UPPER = 4000d0 ! [nm] wavelength range of the Global SW Irradiation measurements

NAMELIST /NML_ENV/  ENV_file, LOSSES_FILE, CHAMBER_FLOOR_AREA, CHAMBER_CIRCUMFENCE, CHAMBER_HEIGHT, &
                    EDDYK,ustar,ALPHAWALL,Cw_eqv, spectrumfile,SWR_IS_ACTINICFLUX,&
                    SWR_IN_LOWER,SWR_IN_UPPER

! MCM INPUT
CHARACTER(len=256)  :: MCM_file = ''
NAMELIST /NML_MCM /MCM_file! ,  MCM_path

! MISC OPTIONS
real(dp)  :: lat              ! Latitude for Photochemistry
real(dp)  :: lon              ! Longitude for Photochemistry
real(dp)  :: CH_Albedo = 2d-1 ! Albedo
real(dp)  :: DMA_f = 0
real(dp)  :: resolve_BASE_precision = 1d-2
CHARACTER(3) :: Fill_formation_with = ''
INTEGER   :: JD = -1
Logical   :: skip_acdc = .True.
INTEGER   :: wait_for = 0 ! -1 for no pause, 0 for indefinite and positive value for fixed amount of seconds
CHARACTER(1000)  :: Description = '*'
CHARACTER(len=256)      :: GR_sizes = '6d-9 20d-9 30d-9'
NAMELIST /NML_MISC/ lat, lon, wait_for, Description, CH_Albedo, DMA_f, resolve_BASE_precision, Fill_formation_with, skip_acdc, &
                    GR_sizes

! Logical                 :: VAP_logical = .True. deprecated
Logical                 :: Use_atoms = .False.
CHARACTER(len=256)      :: Vap_names,Vap_file
CHARACTER(len=256)      :: Vap_atoms = ''

NAMELIST /NML_VAP/ Use_atoms, Vap_names, Vap_atoms !, Vap_props

CHARACTER(1000) :: INITIALIZE_WITH = ''
CHARACTER(1000) :: INIT_ONLY = ''
CHARACTER(1000) :: HARD_CORE = 'GENERIC' ! define what compound is used to initialize particles
INTEGER  :: limit_vapours = 999999
INTEGER  :: INITIALIZE_FROM = 0
Logical  :: use_raoult = .True.

! if true, will not save condensible vapour concentration in Particles.nc. They will always be saved also in Chemistry.nc
Logical  :: DONT_SAVE_CONDENSIBLES     = .False. ! DEPRECATED, has no effect
real(dp) :: dmps_tres_min              = 10.
real(dp) :: VP_MULTI                   = 1d0
real(dp) :: start_time_s               = 0d0
real(dp) :: LAST_VBS_BINNING_S         = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: END_DMPS_PARTIAL           = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: FLOAT_CHEMISTRY_AFTER_HRS  = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: FLOAT_CONC_AFTER_HRS       = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: FLOAT_EMIS_AFTER_HRS       = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: dmps_multi                 = 1d6 ! Multiplicator to convert dmps linear concentration to #/m^3
real(dp) :: SURFACE_TENSION            = 0.05 ! Common surface tension /surface energy density N/m^2 or J/m^3
                                              ! Will be overwritten if vapours-file contains 5th column with compound specific surface tension
real(dp) :: ORGANIC_DENSITY            = 1400d0 ! Common density for organic particle (liquid) phase compounds kg/m^3
real(dp) :: HARD_CORE_DENSITY          = 1400d0 ! Common density for organic particle (liquid) phase compounds kg/m^3
real(dp) :: NPF_DIST                   = 1.15d0 ! Maximum relative diameter (w.r.t. to minimum diameter) where new particles are distributed
Logical  :: NO2_IS_NOX                 = .false.
Logical  :: Kelvin_taylor              = .false.
Logical  :: Kelvin_exp                 = .true.
Logical  :: USE_RH_CORRECTION          = .true.
LOGICAL  :: TEMP_DEP_SURFACE_TENSION   = .False.
LOGICAL  :: use_diff_dia_from_diff_vol = .False.
LOGICAL  :: PARAM_AGING                = .False.
real(dp) :: AGING_HL_HRS               = 1d3 ! [hrs] common aging halflife in particle phase organics.
INTEGER  :: Aging_exponent             = 1
REAL(dp), ALLOCATABLE :: GR_bins(:)  ! used for GR calculation [m]
Logical  :: CALC_GR                 = .True.
Logical  :: ENABLE_END_FROM_OUTSIDE = .True.
Logical  :: Use_old_composition     = .false.
Logical  :: NETCDF_OUT              = .true.
Logical  :: FINAL_CHEM_TXT          = .false.
Logical  :: VBS_ONLY                = .false.
Logical  :: VBS_CSAT                = .false.
Logical  :: PP_H2SO4_TO_AMM_SULFATE = .false.
INTEGER  :: BINARY_FILE             = 0

! First one is the Global timestep lower limit, three four are upper limits for individual processes
real(dp) :: DT_UPPER_LIMIT(3)       = [150d0,150d0,150d0]
real(dp) :: Limit_for_Evaporation   = 0_dp ! not in use
real(dp) :: alpha_coa               = 1d0 ! Accomodation coefficient for coagulation
! real(dp) :: DP_large_particle_limit = 1d0 ! smallest diameter vith viscosity effects [m]
! real(dp) :: k_large_particle_limit  = 0d0 ! expoent for viscosity effect: activity = activity * (1 + log10(diam/viscous_diam)**viscous_exp)
real(dp) :: MIN_CONCTOT_CC_FOR_DVAP = 1d3 ! [#/cm3] lower threshold of concentration, which is checked in optimized time step.
                                          ! Here the idea is that we don't worry about changes of gases whic only exist in so
                                          ! small concentrations
! CHARACTER(len=500)  :: mmodal_input_ems  = ''       ! String defining the Modal PSD emissions
! INTEGER             :: mmodal_input_inuse_ems  = -1 ! Switch toggling the Modal PSD
! CHARACTER(4)  :: nonhomspara = 'mohr'
! CHARACTER(4)  :: homspara = 'stol'

! defined in Constants: Logical  :: NO_NEGATIVE_CONCENTRATIONS = .true.
REAL(dp) :: factorsForReactionRates(NREACT) = 1d0   ! factors to modify chemical reaction rates (optionally)
REAL(dp) :: START_CHEM = -9.9d10 ! Start and stop times in seconds for chem and aero
REAL(dp) :: STOP_CHEM  =  9.9d10 ! Start and stop times in seconds for chem and aero
REAL(dp) :: START_AER  = -9.9d10 ! Start and stop times in seconds for chem and aero
REAL(dp) :: STOP_AER   =  9.9d10 ! Start and stop times in seconds for chem and aero

INTEGER,PARAMETER :: NVBS = 6
!                                ULVOC   |  ELVOC   |   LVOC   |  SVOC   |  IVOC   |   REST
!                                      -8.5       -4.5       -0.5       2.5       6.0
REAL(DP) :: VBS_LIMITS(NVBS-1) = [-8.5,-4.5,-0.5,2.5,6.0]
REAL(DP) :: AGING_K(NVBS) = [1.6e-23,1.6e-23,1.6e-23,1.6e-23,1.6e-23,1.6e-23]
CHARACTER(8) :: VBS_NAMES(NVBS) = ['ULVOC   ','ELVOC   ','LVOC    ','SVOC    ','IVOC    ','REST    ']

NAMELIST /NML_CUSTOM/ use_raoult, dmps_tres_min, &
                      start_time_s, dmps_multi, INITIALIZE_WITH,INITIALIZE_FROM, VP_MULTI,&
                      DONT_SAVE_CONDENSIBLES, limit_vapours, END_DMPS_PARTIAL,NO2_IS_NOX,&
                      NO_NEGATIVE_CONCENTRATIONS, FLOAT_CHEMISTRY_AFTER_HRS, USE_RH_CORRECTION, &
                      TEMP_DEP_SURFACE_TENSION, use_diff_dia_from_diff_vol, DT_UPPER_LIMIT, ENABLE_END_FROM_OUTSIDE, &
                      Limit_for_Evaporation,MIN_CONCTOT_CC_FOR_DVAP, Use_old_composition, alpha_coa, Kelvin_taylor,&
                      SURFACE_TENSION, HARD_CORE,ORGANIC_DENSITY,HARD_CORE_DENSITY,FLOAT_CONC_AFTER_HRS,INIT_ONLY, &
                      FLOAT_EMIS_AFTER_HRS,NPF_DIST, PARAM_AGING, AGING_HL_HRS,FINAL_CHEM_TXT,NETCDF_OUT,VBS_LIMITS,VBS_NAMES,&
                      VBS_ONLY,PP_H2SO4_TO_AMM_SULFATE,VBS_CSAT, BINARY_FILE, &
                      LAST_VBS_BINNING_S,Aging_exponent,AGING_K, START_CHEM, STOP_CHEM, START_AER, STOP_AER

! ==================================================================================================================
! Define change range in percentage
REAL(dp), DIMENSION(2) :: Ddiam_range = [5d-1, 3d0] ! -> minimum and maximum relative change in particle diameter
REAL(dp), DIMENSION(2) :: Dpnum_range = [5d-1, 3d0] ! -> minimum and maximum relative change in particle number
REAL(dp), DIMENSION(2) :: Dvapo_range = [5d-1, 3d0] ! -> minimum and maximum relative change in particle concentration
! Defines the minimum/maximum relative change caused by a process within a timestep
NAMELIST /NML_PRECISION/ Ddiam_range,Dpnum_range,Dvapo_range

INTEGER, ALLOCATABLE        :: ACDC_SYSTEMS(:)
CHARACTER(250), ALLOCATABLE :: ACDC_links(:)

NAMELIST /NML_ACDC/ ACDC_SYSTEMS, ACDC_links
! ==================================================================================================================
! ! Options for screen output
! LOGICAL :: clusterfractions,jions,timestep_multipliers,time_efficiency,Jorganic,GR
! NAMELIST /NML_SCREENPRINTS/ clusterfractions,jions,timestep_multipliers,time_efficiency,Jorganic,GR
! Relative path to NAMES.DAT

! ==================================================================================================================
CHARACTER(200) :: NAMESDAT = 'ModelLib/required/NAMES.dat'
CHARACTER(200) :: INORGANIC = 'ModelLib/required/INORGANIC.dat'
CHARACTER(200) :: XTRASDAT = 'ModelLib/required/AEMS.dat'

NAMELIST /NML_NAMES/ NAMESDAT,INORGANIC,XTRASDAT
! ==================================================================================================================

contains

subroutine READ_rates(Chem)
    real(dp) :: RCONST(NREACT) = 1d0 ! This variable is only used in this subroutine; not the same as in the chemistry
    integer :: iii, k
    CHARACTER(len=*),INTENT(IN) :: Chem

    NAMELIST /NML_RATES/  RCONST
    print FMT_MSG, 'Looking for '//'src/chemistry/'//TRIM(ADJUSTL(Chem))//'/RATES.dat'

    OPEN(UNIT=888, FILE='src/chemistry/'//TRIM(ADJUSTL(Chem))//'/RATES.dat', STATUS='OLD', ACTION='READ', iostat=iii)
    if (iii==0) THEN
        do k=1, ROWCOUNT(888)
            READ(888,NML = NML_RATES, IOSTAT=iii)
        end do
        print FMT_MSG, 'Reaction rates are complimented from RATES.dat'
        factorsForReactionRates = RCONST
    else
        print FMT_MSG, 'Using hardcoded reaction rates'
    END IF

end subroutine READ_rates
! ======================================================================================================================
! Subroutine reads all user input from INITFILE and input files and processes the input etc.
! ......................................................................................................................
subroutine READ_INPUT_DATA()
    IMPLICIT NONE
    CHARACTER(len=256)                :: buf
    integer, allocatable              :: Natoms(:,:)
    integer                           :: ioi, ii, ind_core
    integer                           :: i, j, k, jj, path_l(2), N_Xtr = 0
    integer                           :: rows, cols, class
    LOGICAL                           :: elements_missing = .false.
    real(dp)                          :: molar_mass, parameter_A, parameter_B,SURFACE_TENSION_buf=-1d0, ALPHAWALL_buf=-1d0,Henry, rho
    CHARACTER(len=256)                :: species_name
    CHARACTER(len=20), ALLOCATABLE    :: atoms_name(:)
    CHARACTER(len=2)    :: noacd
    INTEGER             :: nr_of_acdc_modules = 0   ! Number of ACDC submodules NOTE this number is (and should be) defined in
                                                    ! makefile variable NMACDC, because the number of ACDC systems is essentially
                                                    ! only limited upon compilation)
    TYPE(vapourfile_input)  :: vap_in,default_vap
    LOGICAL                 :: CSV_INPUT = .false.

    default_vap%Compound = 'NONAME'
    default_vap%Mass     = 0d0
    default_vap%par_A    = 0d0
    default_vap%par_B    = 0d0
    default_vap%st       = SURFACE_TENSION
    default_vap%alpha_w  = ALPHAWALL
    default_vap%density  = ORGANIC_DENSITY

! This block is handled by C preprocessor --------------------------!
#ifdef NMACDC
    noacd = NMACDC
    read(noacd,'(i2)') nr_of_acdc_modules
#endif
! end of C preprocessor  -------------------------------------------!

    if (nr_of_acdc_modules < 2) THEN
        print FMT_FAT0, 'Number of ACDC submodules was not (correctly) defined in makefile (variable -DNMACDC).'
        stop
    ELSE
        ALLOCATE(G_ACDC(nr_of_acdc_modules))
        ALLOCATE(ACDC_SYSTEMS(nr_of_acdc_modules))
        ALLOCATE(ACDC_links(nr_of_acdc_modules))
        ! These are just default values for backwards compatibility and get overrided if the used has defined them.
        ACDC_SYSTEMS        = 0
        ACDC_SYSTEMS(1:2)   = 1
        ACDC_links(1)       = 'A H2SO4 N NH3'
        ACDC_links(2)       = 'A H2SO4 D DMA'
    End IF

    CALL READ_NAMESDAT

    ! CHECK HOW MANY POSSIBLE INPUT VARIABLES (METEOROLOGICAL, MCM ETC.) THERE ARE IN THE MODEL
    OPEN(800, file=TRIM(NAMESDAT), ACTION='READ', status='OLD', iostat=ioi)
    write(*,FMT_HDR) 'READING COMPOUND NAMES AND ORDER FROM '
    write(*,FMT_SUB) TRIM(NAMESDAT)
    write(*,FMT_HDR) ''
    call handle_file_io(ioi, TRIM(NAMESDAT), 'This file is essential.')

    ! Check the number of rows in NAMESDAT and close file
    N_VARS = rowcount(800)
    CLOSE(800)

    OPEN(800, file=TRIM(XTRASDAT), ACTION='READ', status='OLD', iostat=ioi)
    N_XTRS = rowcount(800)
    CLOSE(800)

    N_VARS = N_VARS + N_XTRS

    ! BASED ON N_VARS, ALLOCATE AND INITIALIZE VECTORS
    ALLOCATE(MODS(N_VARS))
    ALLOCATE(INDRELAY_CH(N_VARS))
    ALLOCATE(INDRELAY_TIED(N_VARS))
    INDRELAY_CH = 0
    INDRELAY_TIED = 0

    CALL READ_INIT_FILE
    ! Declare named indices for convenience in deving
    CALL NAME_MODS_SORT_NAMED_INDICES
    ! Save Initfile options to backup file
    CALL CREATE_DIRECTORIES
    ! Time conversions, Julian Day etc.
    CALL PUT_USER_SUPPLIED_TIMEOPTIONS_IN_GTIME
    ! check if some concentratinos are used only for intialization
    CALL PARSE_INIT_ONLY
    ! Report what is being read from where
    CALL REPORT_INPUT_COLUMNS_TO_USER
    ! Determine what to do with shortwave input
    if (spectrumfile == '') spectrumfile = "ModelLib/Photolyse/Spectra/glob_swr_distr.txt"
    CALL SW_PP

    ! Here we turn submodules on or off based on other options
    If (USE_SPEED)         OPTIMIZE_DT  = .true. ! for backward compatibility
    if (Kelvin_taylor)     Kelvin_exp   = .false.
    if (.not.Kelvin_taylor)Kelvin_exp   = .true.
    if (LOSSES_FILE /= '') Deposition   = .true.
    If (.not.Aerosol_flag) Condensation = .false.
    If (.not.Aerosol_flag) Coagulation  = .false.
    If (.not.Aerosol_flag) Deposition   = .false.
    if (.not.Condensation) CALC_GR = .false.
    H2SO4_ind_in_chemistry = IndexFromName( 'H2SO4', SPC_NAMES )
    OH_ind_in_chemistry = IndexFromName( 'OH', SPC_NAMES )

    ! dmps_tres_min is being phazed out and replaced/overridden by dmps_interval
    if (dmps_interval>0d0) dmps_tres_min = dmps_interval

    ! ALLOCATE CONC_MAT Currently both files need to have same time resolution FIX THIS SOME DAY!
    ! The idea here is to count the rows to get time and allocate CONCMAT and TIMEVEC
    IF ((ENV_file /= '') .or. (MCM_file /= '')) THEN
        IF (ENV_file /= '') THEN
            OPEN(unit=801, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
            CALL handle_file_io(ioi, ENV_file, 'stop')
        ELSE
            OPEN(unit=801, File=TRIM(MCM_file), ACTION='READ', STATUS='OLD',iostat=ioi)
            CALL handle_file_io(ioi, MCM_file, 'stop')
        END IF
        ! Now we can allocate CONC_MAT and TIMEVEC and close the input file for now
        ALLOCATE(CONC_MAT(ROWCOUNT(801,'#'),N_VARS))
        ALLOCATE(TIMEVEC(ROWCOUNT(801,'#')))
        CLOSE(801)

    ! Deal with a situation where we have no input. We still need conc_mat and timevec.
    ELSE
        ALLOCATE(CONC_MAT(2,N_VARS))
        ALLOCATE(TIMEVEC(2))
        TIMEVEC = (/0d0, GTIME%SIM_TIME_H/)
    END IF
    ! Initialize concentrations before reading them from input files
    CONC_MAT = 0d0

    ! READ ENV INPUT
    if (ENV_file /= '') THEN
        OPEN(unit=801, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
        CALL handle_file_io(ioi, ENV_file, 'Terminating on CONCMAT allocation')
        rows = ROWCOUNT(801,'#')
        if ('.CSV' == UCASE(ENV_file(len(TRIM(ENV_file))-3:len(TRIM(ENV_file))))) THEN
          cols = COLCOUNT(801, ',')
        ELSE
          cols = COLCOUNT(801)
        END IF
        ALLOCATE(INPUT_ENV(rows,cols))
        INPUT_ENV = 0
        call FILL_INPUT_BUFF(801,cols,INPUT_ENV,ENV_file)
        timevec = INPUT_ENV(:,1)
        CLOSE(801)
    END IF

    ! READ MCM INPUT
    if (MCM_file /= '') THEN
        OPEN(unit=801, File=TRIM(MCM_file), ACTION='READ', STATUS='OLD')
        rows = ROWCOUNT(801,'#')
        if ('.CSV' == UCASE(MCM_file(len(TRIM(MCM_file))-3:len(TRIM(MCM_file))))) THEN
          cols = COLCOUNT(801, ',')
        ELSE
          cols = COLCOUNT(801)
        END IF
        allocate(INPUT_MCM(rows,cols))
        INPUT_MCM = 0
        call FILL_INPUT_BUFF(801,cols,INPUT_MCM,MCM_file)
        timevec = INPUT_MCM(:,1)
        CLOSE(801)
    END IF

    CALL PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)

    ! check IF dmps data is used or not. If no then do nothing.
    IF (USE_DMPS .and. Aerosol_flag) then
        write(*,FMT_MSG) 'Reading DMPS files '// TRIM(DMPS_file)
        CALL PARSE_PARTICLE_GRID(DMPS_file, BG_PAR)
        BG_PAR%name = 'BACKGROUND_CONC'
    END IF

IF ((TRIM(LOSSES_FILE) /= '') .and. Aerosol_flag) THEN
    if (LOSSES_FILE(1:1) == '#') THEN
        READ(LOSSES_FILE(2:), *) CONSTANT_PAR_LOSS_RATE
    ELSE
        CALL PARSE_PARTICLE_GRID(LOSSES_FILE, PAR_LOSSES)
    END IF
END IF

IF (TRIM(extra_particles) /= '') THEN
    ! First we open the extra particle files to count the dimensions needed for the matrix
    OPEN(unit=801, File=TRIM(extra_particles) , STATUS='OLD', iostat=ioi)
    N_Xtr = ROWCOUNT(801)
    allocate(XTRAS(N_Xtr))

    PRINT FMT_SUB, 'Reading XTRAS: '//TRIM(i2chr(N_Xtr))//' lines'

    DO I=1,N_Xtr
        allocate(XTRAS(I)%options(n_xpar_options))

        read(801,'(a)') buf

        ! Get the indexes for slicing, path, name and options from each line of the input file
        path_l = 0
        k = 1
        DO J=2, LEN(TRIM(BUF))
        if (BUF(J:J) == ' ' .and. BUF(J-1:J-1) /= ' ' .and. k<3) THEN
            path_l(k) = J
            k = k+1
        END IF
        END DO

        XTRAS(I)%name = ADJUSTL(TRIM(BUF(path_l(1):path_l(2))))
        read(BUF(path_l(2):), *) XTRAS(I)%options

        CALL PARSE_PARTICLE_GRID(buf(1:path_l(1)), XTRAS(I))
        print FMT_MSG, 'Extra particle input for '//XTRAS(i)%name

    END DO

END IF

print FMT_LEND,

CALL PARSE_ACDC_SYSTEMS

IF (Aerosol_flag) then

    CALL PARSE_MULTIMODAL
    CALL PARSE_MULTIMODAL_EMS

    if (CALC_GR) CALL PARSE_GR_SIZES

    if (PSD_MODE == 1) write(*,FMT_MSG) 'Using fully stationary PSD scheme with '//TRIM(i2chr(n_bins_par))//' bins.'
    if (PSD_MODE == 2) write(*,FMT_MSG) 'Using fixed grid/moving average PSD scheme with '//TRIM(i2chr(n_bins_par))//' bins.'
    print FMT_LEND,
    IF (condensation) THEN
      Vap_file = Vap_names
      if ('.CSV' == UCASE(Vap_file(len(TRIM(Vap_file))-3:len(TRIM(Vap_file))))) CSV_INPUT = .true.
      write(*,FMT_MSG) 'Reading Vapour name file '// TRIM(Vap_file)
    else
      Vap_file = 'ModelLib/required/empty.dat'
      Use_atoms = .false.
    end if
    OPEN(unit=802, File= TRIM(Vap_file) , STATUS='OLD', iostat=ioi)
    call handle_file_io(ioi, Vap_file, &
        'If Condensation is used, "Vapour file" must be properly defined (in tab "Aerosol").')

    rows = ROWCOUNT(802)
    cols = COLCOUNT(802)

    VAPOUR_PROP%n_cond_org = 0
    do j = 1,rows
        read(802,*,iostat=ioi) species_name
        if (j<=limit_vapours .or. j==rows) THEN
            if ((TRIM(species_name) == 'HOA')         & ! HOA is sometimes used in PRAM
                .or.(TRIM(species_name) == 'GENERIC') & ! GENERIC is not active in chemistry
                ! We want to add H2SO4 to the end of vapours, and it should not be in vapour file anyway
                .or.(TRIM(species_name) == 'H2SO4') ) &
                cycle
            k = IndexFromName( TRIM(species_name), SPC_NAMES )
            if (k>0) THEN
                VAPOUR_PROP%n_cond_org = VAPOUR_PROP%n_cond_org + 1
            end if
        end if
    end do
    REWIND(802)
    ! Here we account for the non-volatile pseudocompound GENERIC
    VAPOUR_PROP%n_cond_org = VAPOUR_PROP%n_cond_org + 1

    write(*,FMT_MSG) 'Reading list of particle phase inorganics '// TRIM(INORGANIC)
    OPEN(unit=803, File= TRIM(INORGANIC) , STATUS='OLD', iostat=ioi)
    call handle_file_io(ioi, TRIM(INORGANIC), &
        'Could not access list of inorganics '//TRIM(INORGANIC))
    rows = ROWCOUNT(803)

    VAPOUR_PROP%n_cond_ino = 0
    VAPOUR_PROP%ind_H2SO4 = 0

    do j = 1,rows
        read(803,*,iostat=ioi) species_name
        k = IndexFromName( TRIM(species_name), SPC_NAMES )
        if (k>0) THEN
            VAPOUR_PROP%n_cond_ino = VAPOUR_PROP%n_cond_ino + 1
            if (TRIM(UCASE(species_name)) == 'H2SO4') &
              VAPOUR_PROP%ind_H2SO4 = VAPOUR_PROP%n_cond_ino + VAPOUR_PROP%n_cond_org
        end if
    end do

    REWIND(803)
    ! Here we add place for sulfuric acid, non-organic
    VAPOUR_PROP%n_cond_tot = VAPOUR_PROP%n_cond_org + VAPOUR_PROP%n_cond_ino

    allocate(VAPOUR_PROP%Vapour_names (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%molar_mass   (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%psat_a       (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%psat_b       (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%molec_mass   (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%molec_volume (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%density      (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%surf_tension (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%diff         (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%c_speed      (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%c_sat        (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%cond_type    (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%molec_dia    (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%mfractions   (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%alpha        (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%alphawall    (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%diff_vol     (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%diff_dia     (VAPOUR_PROP%n_cond_tot) )
    allocate(VAPOUR_PROP%VBS_BINS     (VAPOUR_PROP%n_cond_org-1, NVBS) )
    VAPOUR_PROP%VBS_BINS = 0
    VAPOUR_PROP%ind_g_ULVOC = IndexFromName( 'g_ULVOC', SPC_NAMES)
    VAPOUR_PROP%ind_g_ELVOC = IndexFromName( 'g_ELVOC', SPC_NAMES)
    VAPOUR_PROP%ind_g_LVOC = IndexFromName( 'g_LVOC', SPC_NAMES)
    VAPOUR_PROP%ind_g_SVOC = IndexFromName( 'g_SVOC', SPC_NAMES)
    VAPOUR_PROP%ind_g_IVOC = IndexFromName( 'g_IVOC', SPC_NAMES)
    VAPOUR_PROP%ind_g_REST = IndexFromName( 'g_REST', SPC_NAMES)

    ! This is the vector that combines chemistry and condensible vapours.
    ! -1 because GENERIC is not in gas phase. Ever. H2SO4 is picked from chemistry manually
    ALLOCATE(index_cond(VAPOUR_PROP%n_cond_tot))
    ! ALLOCATE(index_inorg(VAPOUR_PROP%n_cond_org+1:VAPOUR_PROP%n_cond_tot))

    index_cond = 0
    ! index_inorg = 0

    if (USE_RH_CORRECTION) THEN
        ! These are only allocated for sulfuric acid, maybe nitric acid in the future.
        ! Uses same index as in VAPOUR_PROP%ind_H2SO4
        allocate(VAPOUR_PROP%wet_dia(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_cond_tot))
        allocate(VAPOUR_PROP%wet_mass(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_cond_tot))
    END IF


    IF (condensation) print FMT_SUB, 'Compounds picked from Vapours file: '//TRIM(i2chr(VAPOUR_PROP%n_cond_org))
    IF (condensation) print FMT_SUB, 'Total number of condensibles      : '//TRIM(i2chr(VAPOUR_PROP%n_cond_tot))

    ! Reading the vap names and vap vapour_properties
    VAPOUR_PROP%ind_GENERIC = VAPOUR_PROP%n_cond_org
    VAPOUR_PROP%Mfractions  = 0.0
    VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0 !

    ! ---------------------------------------------------------------------
    ! ORGANIC VAPOUR PROPERTIES
    ! ---------------------------------------------------------------------
    rows = ROWCOUNT(802)
    ii = 1
    do j = 1, rows
        if (CSV_INPUT) THEN
          read(802,*,iostat=ioi) vap_in
          species_name        = vap_in%Compound
          molar_mass          = vap_in%Mass
          parameter_A         = vap_in%par_A
          parameter_B         = vap_in%par_B
          SURFACE_TENSION_buf = vap_in%st
          ALPHAWALL_buf       = vap_in%alpha_w
          rho                 = vap_in%density
          vap_in              = default_vap
        ELSE
          if (cols==4) read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B
          if (cols==5) read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B, SURFACE_TENSION_buf
          if (cols==6) read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B, SURFACE_TENSION_buf, ALPHAWALL_buf
        end IF
        if (j<=limit_vapours .or. j==rows) THEN ! j==rows is because "GENERIC" is selected even if limitvapours<rows

            ! Check if the compounds exists in Chemistry and only then add to vapours
            k = IndexFromName( species_name, SPC_NAMES )

            if (k>0.and..not.(TRIM(species_name)=='GENERIC')) index_cond(ii) = k

            if ((k>0).or.(TRIM(species_name)=='GENERIC')) THEN
                ! fill the hash table for vapour index -> chemistry index

                VAPOUR_PROP%molar_mass(ii)   = molar_mass *1D-3 ! kg/mol
                VAPOUR_PROP%psat_a(ii)       = parameter_A
                VAPOUR_PROP%psat_b(ii)       = parameter_B
                VAPOUR_PROP%vapour_names(ii) = TRIM(species_name)
                VAPOUR_PROP%molec_mass(ii)   = VAPOUR_PROP%molar_mass(ii)/Na  !kg/#
                if (TRIM(species_name)=='GENERIC') THEN
                    VAPOUR_PROP%cond_type(ii) = 0  ! Generic non-evaporating
                    VAPOUR_PROP%density(ii)   = HARD_CORE_DENSITY  ! kg/m3
                ELSE
                    VAPOUR_PROP%cond_type(ii) = 1  ! Organic compound
                    VAPOUR_PROP%density(ii)   = ORGANIC_DENSITY  ! kg/m3
                END IF
                VAPOUR_PROP%molec_volume(ii) = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
                VAPOUR_PROP%diff_vol(ii)     = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)

                if (cols==4) VAPOUR_PROP%surf_tension(ii) = SURFACE_TENSION
                if (cols>=5.or.CSV_INPUT) THEN
                    if (SURFACE_TENSION_buf<0d0) VAPOUR_PROP%surf_tension(ii) = abs(SURFACE_TENSION_buf)*SURFACE_TENSION
                    if (SURFACE_TENSION_buf>=0d0) VAPOUR_PROP%surf_tension(ii) = SURFACE_TENSION_buf
                end if
                SURFACE_TENSION_buf = -1d0

                if (cols<6) VAPOUR_PROP%alphawall(ii) = ALPHAWALL
                if (cols>=6.or.CSV_INPUT) THEN
                    if (ALPHAWALL_buf<0d0) VAPOUR_PROP%alphawall(ii) = abs(ALPHAWALL_buf)*ALPHAWALL
                    if (ALPHAWALL_buf>=0d0) VAPOUR_PROP%alphawall(ii) = ALPHAWALL_buf
                end if
                ALPHAWALL_buf = -1d0

                VAPOUR_PROP%alpha(ii)         = 1.0
                ! this is just initial value, always gets updated with T
                VAPOUR_PROP%c_sat(ii)         = saturation_conc_m3(VAPOUR_PROP%psat_a(ii),VAPOUR_PROP%psat_b(ii), 293.15d0)

                ii = ii + 1
            END IF
        END IF
    end do
    close(802)

    ! READ INORGANICS
    rows = ROWCOUNT(803)
    ii = VAPOUR_PROP%n_cond_org + 1
    do j = 1, rows
        read(803,*,iostat=ioi)   species_name, molar_mass, Henry, rho, class
        ! Check if the compounds exists in Chemistry and only then add to vapours
        k = IndexFromName( species_name, SPC_NAMES )

        if (k>0) THEN
            index_cond(ii) = k

            ! fill the hash table for vapour index -> chemistry index

            VAPOUR_PROP%molar_mass(ii)   = molar_mass *1D-3 ! kg/mol
            VAPOUR_PROP%psat_a(ii)       = 1d0
            VAPOUR_PROP%psat_b(ii)       = 20000d0
            VAPOUR_PROP%vapour_names(ii) = TRIM(species_name)
            VAPOUR_PROP%molec_mass(ii)   = VAPOUR_PROP%molar_mass(ii)/Na  !kg/#
            VAPOUR_PROP%cond_type(ii)    = class  ! some other compound
            VAPOUR_PROP%density(ii)      = rho  ! kg/m3

            VAPOUR_PROP%molec_volume(ii) = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
            VAPOUR_PROP%diff_vol(ii)     = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
            VAPOUR_PROP%surf_tension(ii) = 0d0
            VAPOUR_PROP%alpha(ii)        = 1.0
            VAPOUR_PROP%c_sat(ii)        = 0d0

            ii = ii + 1
        END IF
    end do
    close(803)


    ! In case GENERIC was not in Vapour file (should not happen if the file was from the GUI), add GENERIC with default values
    if (VAPOUR_PROP%vapour_names(VAPOUR_PROP%n_cond_org) /= 'GENERIC') THEN
        IF (condensation) print FMT_WARN0, 'The vapour file did not contain GENERIC, adding it now. You should update your vapour file.'
        ii = VAPOUR_PROP%n_cond_org
        VAPOUR_PROP%vapour_names(ii)  = 'GENERIC'
        VAPOUR_PROP%cond_type(ii)     = 0  ! Generic non-evaporating
        VAPOUR_PROP%molar_mass(ii)    = 437.0 * 1d-3
        VAPOUR_PROP%psat_a(ii)        = 10
        VAPOUR_PROP%psat_b(ii)        = 1d4
        VAPOUR_PROP%density(ii)       = HARD_CORE_DENSITY
        VAPOUR_PROP%molec_mass(ii)    = VAPOUR_PROP%molar_mass(ii)/Na
        VAPOUR_PROP%molec_volume(ii)  = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
        VAPOUR_PROP%diff_vol(ii)      = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
        VAPOUR_PROP%surf_tension(ii)  = SURFACE_TENSION
        VAPOUR_PROP%alpha(ii)         = 1.0
        ! this is just initial value, always gets updated with T
        VAPOUR_PROP%c_sat(ii)         = saturation_conc_m3(VAPOUR_PROP%psat_a(ii),VAPOUR_PROP%psat_b(ii), 293.15d0)
    END IF

    VAPOUR_PROP%Mfractions  = 0.0
    if (TRIM(HARD_CORE) == '') THEN
        VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0 !
    ELSE
        ind_core = IndexFromName(HARD_CORE,VAPOUR_PROP%vapour_names)
        if (ind_core>0 .and. ind_core<=VAPOUR_PROP%n_cond_tot) THEN
            VAPOUR_PROP%Mfractions(ind_core) = 1d0
        ELSE
            VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0
        END IF
    END IF

    if (Use_atoms) THEN
        OPEN(unit=804, File=TRIM(Vap_atoms) , STATUS='OLD', iostat=ioi)
        call handle_file_io(ioi, Vap_atoms, 'Terminating the program.')
        write(*,FMT_MSG) 'Reading the list of elemental composition: '// TRIM(Vap_atoms)

        allocate(Natoms(5,ROWCOUNT(804))) ! C,O,N,H
        allocate(atoms_name(ROWCOUNT(804)))

        Natoms = 0

        DO j=1,ROWCOUNT(804)
            READ(804,*, iostat=ioi) atoms_name(j), molar_mass, Natoms(:,j)
        END DO

        CLOSE(804)

        DO j=1,VAPOUR_PROP%n_cond_org
            jj = IndexFromName(VAPOUR_PROP%vapour_names(j), atoms_name)
            if (jj>0) THEN
                vapour_prop%diff_vol(j) = (Natoms(1,jj)*15.9D0 + Natoms(2,jj)*6.11D0 &
                                          + Natoms(4,jj)*2.31D0 + Natoms(3,jj)*4.54D0) + Natoms(5,jj) * 22.9D0 ![Å^3]
            ELSE
                elements_missing = .true.
            END IF
        END DO
        if (elements_missing) print FMT_WARN0, 'Not all organics had atom content, using generic diameter'

        deallocate(Natoms)
        deallocate(atoms_name)

    end if

    ! ii = VAPOUR_PROP%n_cond_tot
    ! VAPOUR_PROP%vapour_names(ii)  = 'H2SO4'
    ! VAPOUR_PROP%cond_type(ii)     = 2  ! Acid
    ! VAPOUR_PROP%molar_mass(ii)    = 98.0785 * 1d-3
    ! VAPOUR_PROP%psat_a(ii)        = 0
    ! VAPOUR_PROP%psat_b(ii)        = 20000d0
    ! VAPOUR_PROP%density(ii)       = 1830.5 ! kg/m3
    ! VAPOUR_PROP%molec_mass(ii)    = VAPOUR_PROP%molar_mass(ii)/Na
    ! VAPOUR_PROP%molec_volume(ii)  = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
    ! VAPOUR_PROP%diff_vol(ii)      = 4*6.11D0 + 2*2.31D0 + 22.9D0 ! O=4, H=2, S=1
    ! VAPOUR_PROP%surf_tension(ii)  = 0.07
    ! VAPOUR_PROP%molec_dia(ii)     = (6D0 * VAPOUR_PROP%molec_volume(ii) / pi )**(1D0/3D0)  ! molecular diameter [m]
    ! VAPOUR_PROP%alpha(ii)         = 1.0
    ! VAPOUR_PROP%c_sat(ii)         = 0.0 ! Sulfuric acid stays put
    ! VAPOUR_PROP%diff_dia(ii)      = VAPOUR_PROP%molec_dia(ii)

    ! Now the volumes are updated, the diameter can be calculated
    VAPOUR_PROP%molec_dia(:ii-1) = (6D0 * VAPOUR_PROP%molec_volume(:ii-1) / pi )**(1D0/3D0)  ! molecular diameter [m]
    if (use_diff_dia_from_diff_vol) THEN
        VAPOUR_PROP%diff_dia(:ii-1) = (6D0 * 1d-30 * VAPOUR_PROP%diff_vol(:ii-1) / pi )**(1D0/3D0)  ! molecular diameter [m]
    ELSE
        VAPOUR_PROP%diff_dia(:ii-1) = VAPOUR_PROP%molec_dia(:ii-1)  ! molecular diameter [m]
    END IF

end if

  CALL CHECK_MODIFIERS ! Print out which modifiers differ from default values
end subroutine READ_INPUT_DATA


!-------------------------------------------------------------------------------
! Reads the path to NAMESDAT
!-------------------------------------------------------------------------------
subroutine READ_NAMESDAT
  implicit none
  integer :: IOS(20)=0, i=1, k
  CALL GETARG(1,Fname_init)
  ! Advice user if no INITFILE was provided
  IF (Fname_init == '') THEN
    write(*,FMT_LEND)
    write(*,FMT_MSG) 'No INITFILE Defined. Send the relative path to proper INITFILE as command line option.'
    CALL GETARG(0,Fname_init)
    write(*,FMT_MSG) 'For example: $ '//TRIM(Fname_init)//' model_init'
    write(*,FMT_LEND)
    STOP
  END IF

  CALL GETARG(2,gui)
  if (gui == '--gui') ingui = .true.

  OPEN(UNIT=888, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=IOS(1))
  ! Handle file not found error
  IF (IOS(1) /= 0) THEN
    write(*,FMT_FAT0) 'There is no INITFILE '//TRIM(ADJUSTL(Fname_init))//', exiting. Good bye.'
    write(*,FMT_LEND)
    STOP
  END IF

  ! if INITFILE was found, we read it. In case there is a problem in namelist filling, give en error.
  write(*,FMT_HDR) 'READING USER DEFINED INTIAL VALUES FROM:'
  write(*,FMT_SUB) TRIM(ADJUSTL(Fname_init))

  ! if (.not. ingui) &
  !   CALL EXECUTE_COMMAND_LINE('sh ModelLib/gui/modules/compatibility_layer.sh '//TRIM(ADJUSTL(Fname_init)))


  do k=1, ROWCOUNT(888); READ(888,NML = NML_NAMES, IOSTAT=IOS(i)) ! #1
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1
  CLOSE(888)

end subroutine READ_NAMESDAT

subroutine READ_INIT_FILE
  implicit none
  integer :: IOS(20)=0, i=1, k
  ! CALL GETARG(1,Fname_init)
  ! ! Advice user if no INITFILE was provided
  ! IF (Fname_init == '') THEN
  !   write(*,FMT_LEND)
  !   write(*,FMT_MSG) 'No INITFILE Defined. Send the relative path to proper INITFILE as command line option.'
  !   CALL GETARG(0,Fname_init)
  !   write(*,FMT_MSG) 'For example: $ '//TRIM(Fname_init)//' model_init'
  !   write(*,FMT_LEND)
  !   STOP
  ! END IF
  !
  ! CALL GETARG(2,gui)
  ! if (gui == '--gui') ingui = .true.
  !
  OPEN(UNIT=888, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=IOS(1))
  ! ! Handle file not found error
  ! IF (IOS(1) /= 0) THEN
  !   write(*,FMT_FAT0) 'There is no INITFILE '//TRIM(ADJUSTL(Fname_init))//', exiting. Good bye.'
  !   write(*,FMT_LEND)
  !   STOP
  ! END IF

  ! ! if INITFILE was found, we read it. In case there is a problem in namelist filling, give en error.
  ! write(*,FMT_HDR) 'READING USER DEFINED INTIAL VALUES FROM:'
  ! write(*,FMT_HDR) TRIM(ADJUSTL(Fname_init))
  !
  ! if (.not. ingui) &
  !   CALL EXECUTE_COMMAND_LINE('sh ModelLib/gui/modules/compatibility_layer.sh '//TRIM(ADJUSTL(Fname_init)))
  !
  !
  ! do k=1, ROWCOUNT(888); READ(888,NML = NML_NAMES, IOSTAT=IOS(i)) ! #1
  ! ! IF (IOS(i) == 0) EXIT;
  ! end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_TIME, IOSTAT=IOS(i)) ! #2
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_Flag, IOSTAT=IOS(i)) ! #3
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_Path, IOSTAT=IOS(i)) ! #4
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MISC, IOSTAT=IOS(i)) ! #5
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_VAP, IOSTAT=IOS(i)) ! #6
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_PARTICLE, IOSTAT=IOS(i)) ! #7
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_ENV, IOSTAT=IOS(i)) ! #8
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MCM, IOSTAT=IOS(i)) ! #9
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MODS, IOSTAT=IOS(i)) ! #10
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_PRECISION, IOSTAT=IOS(i)) ! #11
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_CUSTOM, IOSTAT=IOS(i)) ! #12
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_ACDC, IOSTAT=IOS(i)) ! #13
  ! IF (IOS(i) == 0) EXIT;
  end do; REWIND(888); i=i+1
  REWIND(888); i=i+1

  CLOSE(888)

  IF (INOUT_DIR(len(TRIM(INOUT_DIR)):len(TRIM(INOUT_DIR))) == '/') INOUT_DIR(len(TRIM(INOUT_DIR)):len(TRIM(INOUT_DIR))) = ' '

end subroutine READ_INIT_FILE

subroutine CREATE_DIRECTORIES
  implicit none
  integer :: i
  ! Also save all settings to initfile. Use this file to rerun if necessary
  open(889, file=TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(NUMBER)//'/'//TRIM(RUN_NAME)//'/NMLS.conf', action='WRITE',iostat=i)
  if (i/=0) THEN
    ! Create output directory
    print FMT_MSG, 'Output directory did not exist -> creating...'
    call system('mkdir -p '//TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(NUMBER)//'/'//TRIM(RUN_NAME), I)
    call handle_file_io(I, TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(NUMBER)//'/'//TRIM(RUN_NAME),&
                     'Could not create output directories')
    open(889, file=TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(NUMBER)//'/'//TRIM(RUN_NAME)//'/NMLS.conf', action='WRITE',iostat=i)
    call handle_file_io(I, TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(NUMBER)//'/'//TRIM(RUN_NAME)//'/NMLS.conf',&
    'Could not write files, do you have permissions?')
  END if

  write(889,NML = NML_TIME       ) ! directories and test cases
  write(889,NML = NML_Flag       ) ! flags
  write(889,NML = NML_Path       ) ! time related stuff
  write(889,NML = NML_MISC       ) ! dmps_file information
  write(889,NML = NML_VAP        ) ! environmental information
  write(889,NML = NML_PARTICLE   ) ! MCM_file information
  write(889,NML = NML_ENV        ) ! modification parameters
  write(889,NML = NML_MCM        ) ! misc input
  write(889,NML = NML_PRECISION  ) ! custom input
  write(889,NML = NML_CUSTOM     ) ! custom input
  write(889,NML = NML_ACDC       ) ! ACDC systems
  write(889,NML = NML_MODS       ) ! vapour input
  write(889,NML = NML_NAMES      ) ! vapour input
  close(889)

end subroutine CREATE_DIRECTORIES

subroutine PUT_USER_SUPPLIED_TIMEOPTIONS_IN_GTIME
    implicit none
    INTEGER :: days(12), non_leap_year(12) = ([31,28,31,30,31,30,31,31,30,31,30,31])
    INTEGER :: leap_year(12) = ([31,29,31,30,31,30,31,31,30,31,30,31])
    INTEGER :: y,m,d,ioi

    GTIME%SIM_TIME_H = runtime
    GTIME%SIM_TIME_S = runtime*3600d0
    GTIME%dt = DT
    ! figure out the correct save interval
    IF (FSAVE_DIVISION > 0) THEN
        GTIME%FSAVE_INTERVAL = max(GTIME%dt, INT(GTIME%SIM_TIME_S/GTIME%dt) / FSAVE_DIVISION * GTIME%dt)
    ELSEIF (FSAVE_INTERVAL > 0) THEN
        GTIME%FSAVE_INTERVAL = FSAVE_INTERVAL
    END IF
    GTIME%PRINT_INTERVAL = PRINT_INTERVAL
    GTIME%PRINTACDC = PRINT_ACDC

    if (TRIM(DATE) /= '1800-01-01') THEN
    read(date(1:4),*,iostat=ioi) y
    if (ioi == 0) read(date(6:7),*,iostat=ioi) m
    if (ioi == 0) read(date(9:) ,*,iostat=ioi) d
    if (ioi /= 0) THEN
        if (NUMBER /= '') THEN
            print FMT_HDR,
            print FMT_HDR, 'USING NUMBER INSTEAD OF DATE -> ASSUMING LIGHT DIRECTION IS FROM DIRECTLY UP'
            print FMT_HDR,
        ELSE
            print FMT_WARN0, 'Date provided in the INITFILE is not a proper date'
        END IF
    ELSE
        if (MODULO(y,4) == 0) THEN
            if ((MODULO(y,100) == 0) .and. (MODULO(y,400) /= 0)) THEN
                days = non_leap_year
            ELSE
                days = leap_year
            END IF
        ELSE
            days = non_leap_year
        END IF
        GTIME%JD = sum(days(:m-1)) + d
        print FMT_MSG, 'Date: '//TRIM(date)//' -> Julian Day: '//TRIM(i2chr(GTIME%JD))
    END IF
END if


  ! In case user wants to start at later point, move clock, but only if the movement skips one filesave
    if (start_time_s > 0) THEN
        GTIME = ADD(GTIME, start_time_s)
    end if
end subroutine PUT_USER_SUPPLIED_TIMEOPTIONS_IN_GTIME

subroutine NAME_MODS_SORT_NAMED_INDICES
  implicit none
  INTEGER :: i, j

  DO i=1, N_VARS
    IF (TRIM(MODS(i)%UNIT) == 'X') THEN
        MODS(i)%UNIT = '#'
    ELSE
        MODS(i)%ISPROVIDED = .true.
    END IF
  END DO

  OPEN(800, file=TRIM(NAMESDAT), ACTION='READ', status='OLD')
  DO i = 1,N_VARS - N_XTRS
    READ(800, *) MODS(I)%NAME
    IF (TRIM(MODS(I)%NAME) == 'TEMPK'        ) inm_TempK = i
    IF (TRIM(MODS(I)%NAME) == 'PRESSURE'     ) inm_pres = i
    IF (TRIM(MODS(I)%NAME) == 'REL_HUMIDITY' ) inm_RH = i
    IF (TRIM(MODS(I)%NAME) == 'CONDENS_SINK' ) inm_CS = i
    IF (TRIM(MODS(I)%NAME) == 'CON_SIN_NITR' ) inm_CS_NA = i
    IF (TRIM(MODS(I)%NAME) == 'SW_RADIATION' ) inm_swr = i
    IF (TRIM(MODS(I)%NAME) == 'ION_PROD_RATE') inm_IPR = i
    IF (TRIM(MODS(I)%NAME) == 'H2SO4'        ) inm_H2SO4 = i
    IF (TRIM(MODS(I)%NAME) == 'NH3'          ) inm_NH3 = i
    IF (TRIM(MODS(I)%NAME) == 'DMA'          ) inm_DMA = i
    IF (TRIM(MODS(I)%NAME) == 'SO2'          ) inm_SO2 = i
    IF (TRIM(MODS(I)%NAME) == 'NO'           ) inm_NO = i
    IF (TRIM(MODS(I)%NAME) == 'NO2'          ) inm_NO2 = i
    IF (TRIM(MODS(I)%NAME) == 'CO'           ) inm_CO = i
    IF (TRIM(MODS(I)%NAME) == 'H2'           ) inm_H2 = i
    IF (TRIM(MODS(I)%NAME) == 'O3'           ) inm_O3 = i
    IF (TRIM(MODS(I)%NAME) == 'NUC_RATE_IN ' ) inm_JIN = i
    IF (TRIM(MODS(I)%NAME) == '#'            ) LENV = i
    IF (.not.TRIM(UCASE(MODS(I)%NAME(1:3))) == 'EMI') N_CONCENTRATIONS = i
  END DO
  close(800)

  OPEN(800, file=TRIM(XTRASDAT), ACTION='READ', status='OLD')
  DO j = 1, N_XTRS
    READ(800, *) MODS(i+j-1)%NAME
  END DO
  close(800)

end subroutine NAME_MODS_SORT_NAMED_INDICES


subroutine REPORT_INPUT_COLUMNS_TO_USER
    implicit none
    integer :: i

    DO i=1,N_VARS
        IF (I==1) print FMT_MSG, 'ENV values from '//TRIM(ENV_file)//':'
        IF ((I==LENV) .and. (maxval(MODS(LENV:)%col)>-1)) print FMT_MSG, 'MCM values from '//TRIM(MCM_file)//':'

        IF (MODS(I)%col > -1) THEN
            IF ((TRIM(ENV_file) == '') .and. (I<LENV) .and. MODS(I)%mode==0) THEN
                print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(i2chr(MODS(I)%col))//' but no ENV_FILE'
                print FMT_SUB, 'To avoid accidental error, if referring to file it must exist.'
                STOP  'Change column number to -1, link to another variable, or use parametric input'
            ELSE IF ((TRIM(MCM_file) == '') .and. (I>LENV) .and. MODS(I)%mode==0) THEN
                print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(i2chr(MODS(I)%col))//' but no MCM_FILE'
                print FMT_SUB, 'To avoid accidental error, if referring to file it must exist.'
                STOP  'Change column number to -1, link to another variable, or use parametric input'
            ELSE
                print FMT_SUB, TRIM(MODS(I)%NAME)//' will be read from column: '//TRIM(i2chr(MODS(I)%col))
            END IF
        END IF
    END DO
end subroutine REPORT_INPUT_COLUMNS_TO_USER

SUBROUTINE PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)
    implicit none
    REAL(DP), intent(in)    :: INPUT_ENV(:,:)
    REAL(DP), intent(in)    :: INPUT_MCM(:,:)
    REAL(DP), intent(inout) :: CONC_MAT(:,:)
    integer                 :: i

    DO i=1,N_VARS
        IF ((I < lenv) .and. (ENV_file /= '') .and. (MODS(I)%col > -1)) THEN
            CONC_MAT(:,I) = input_ENV(:,MODS(I)%col)
        END IF

        IF ((I>lenv) .and. (MCM_file /= '') .and. (MODS(I)%col > -1)) THEN
            CONC_MAT(:,I) = input_MCM(:,MODS(I)%col)
        END IF
    END DO

    if (NO2_IS_NOX) CONC_MAT(:,inm_NO2) =  CONC_MAT(:,inm_NO2)- CONC_MAT(:,inm_NO)

END SUBROUTINE PUT_INPUT_IN_THEIR_PLACES



subroutine FILL_INPUT_BUFF(unit,cols,INPUT_BF,Input_file)
    implicit none
    real(dp), intent(inout) :: INPUT_BF(:,:)
    integer, intent(in)     :: cols
    CHARACTER(*)            :: Input_file
    integer                 :: i,j,k,ioi,unit
    CHARACTER(len=6000)     :: dump

    i = 1
    print FMT_MSG, 'Filling input matrices...'
    DO k = 1, ROWCOUNT(unit)
        READ(unit,*, iostat=ioi) (INPUT_BF(i,j),j=1,cols)
        IF ((ioi /= 0) .and. (k==1)) THEN
            REWIND(unit)
            READ(unit,*, iostat=ioi) dump
            print FMT_SUB, 'First row omitted from file "'// TRIM(Input_file) //'".'
        ELSE IF ((ioi /= 0) .and. (k>1)) THEN
            print FMT_WARN1, 'Bad value in file '// TRIM(Input_file) //'". Maybe a non-numeric on line '//i2chr(k)
        ELSE
            i=i+1
        END IF
    END DO
    print FMT_SUB, 'Done filling input matrices...'
end subroutine FILL_INPUT_BUFF


! ================================================================================================
! Print out which modifiers differ from default values. If SHIFTER or MULTIPLYER differ from their
! null value, this subroutine will pick them and print them for the user before the main loop starts
! ================================================================================================
SUBROUTINE CHECK_MODIFIERS()
    IMPLICIT NONE
    type(input_mod) :: test
    integer         :: i,j=0
    CHARACTER(4)    :: cprf

    print FMT_HDR, 'Checking input validity'

    CALL CONVERT_TEMPS_TO_KELVINS
    CALL CONVERT_PRESSURE_AND_VALIDATE_UNITS

    do i=1,size(MODS)
        IF (MODS(i)%MODE > 0) THEN
            print FMT_MSG, '- Replacing input for '//TRIM(MODS(i)%name)//' with parametrized function.'
            j=1
        ELSE
            IF (ABS(MODS(i)%multi - test%multi) > 1d-9) THEN
                print FMT_MSG, '- Multiplying '//TRIM(MODS(i)%name)//' with: '//TRIM(f2chr(MODS(i)%multi))
                j=1
            END IF
            if (ABS(MODS(i)%shift - test%shift) > 1d-9) THEN
                if (TRIM(MODS(i)%UNIT) == '#') THEN
                    cprf = '/cm3'
                else
                    cprf = ''
                end if
                print FMT_MSG, '- Adding a constant to '//TRIM(MODS(i)%name)//', [in '//TRIM(MODS(i)%UNIT)//TRIM(cprf)//']: '//TRIM(f2chr(MODS(i)%shift))
                j=1
            END IF
        END IF
    END DO
    if (j == 1) print FMT_LEND

END SUBROUTINE CHECK_MODIFIERS


! ================================================================================================
! Subroutine converts temperature to Kelvins based on the user input. Since there is more than one
! way to define unit for temperature, this routine tries to tackle with all of them.
! ================================================================================================
SUBROUTINE CONVERT_TEMPS_TO_KELVINS
    IMPLICIT NONE
    CHARACTER(len=1) :: TempUnit = 'X'
    !
    ! if (TRIM(UCASE(MODS(inm_TempK)%UNIT)) == '#') THEN
    !     print FMT_WARN0, "No unit for temperature. Use either 'K' or 'C'. Now assuming Kelvins. This may lead to SIGFPE."
    !     TempUnit = 'K'
    ! ELSEIF (TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'K' .or. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'C' )THEN
    !     TempUnit = TRIM(UCASE(MODS(inm_TempK)%UNIT))
    ! END IF

    TempUnit = TRIM(UCASE(MODS(inm_TempK)%UNIT))

    IF (UCASE(TempUnit) == 'K') THEN
        print FMT_MSG, '- Temperature input in Kelvins.'
    ELSEIF (UCASE(TempUnit) == 'C') THEN
        print FMT_MSG, '- Converting temperature from degrees C -> K.'
        MODS(inm_TempK)%min = MODS(inm_TempK)%min + 273.15d0
        MODS(inm_TempK)%max = MODS(inm_TempK)%max + 273.15d0
        CONC_MAT(:,inm_TempK) = CONC_MAT(:,inm_TempK) + 273.15d0
    ELSE
        print FMT_WARN0, "Could not recognize temperature unit. Use either 'K' or 'C'. Now assuming Kelvins."
    END IF

    ! Check if accidental double conversion is attempted
    IF ((TempUnit == 'C') .and. (  ABS(MODS(inm_TempK)%shift - 273.15)<1d0  )) THEN
        print FMT_WARN1, 'Temperature will be converted to Kelvins, but an additional constant is added: ',MODS(inm_TempK)%shift
    END IF

    MODS(inm_TempK)%UNIT = 'K'

END SUBROUTINE CONVERT_TEMPS_TO_KELVINS


! ================================================================================================
! This subroutine converts pressure from all possible input units to Pa
! ================================================================================================
SUBROUTINE CONVERT_PRESSURE_AND_VALIDATE_UNITS
  IMPLICIT NONE
  INTEGER      :: i
  CHARACTER(5) :: buf

  buf = UCASE(TRIM(MODS(inm_pres)%UNIT))
  if (TRIM(buf) == 'HPA' .or. TRIM(buf) == 'MBAR') THEN
      print FMT_MSG, '- Converting pressure from hPa (mbar) to Pascals.'
    ELSEIF (TRIM(buf) == 'KPA') THEN
      print FMT_MSG, '- Converting pressure from kPa to Pascals.'
    ELSEIF (TRIM(buf) == 'ATM') THEN
      print FMT_MSG, '- Converting pressure from atm to Pascals.'
    ELSEIF (TRIM(buf) == 'BAR') THEN
      print FMT_MSG, '- Converting pressure from bar to Pascals.'
    ELSEIF (TRIM(buf) == 'PA') THEN
      print FMT_MSG, '- Pressure is given in Pascals.'
    else
      if ((MODS(inm_pres)%MODE > 0) .or. (MODS(inm_pres)%col > 1)  .or. (ABS(MODS(inm_pres)%multi - 1d0)>1d-9) .or. (ABS(MODS(inm_pres)%shift)>1d-16)) THEN
          if (TRIM(buf) == '#') THEN
              print FMT_MSG, '- Assuming Pascals for pressure.'
          else
              print FMT_FAT0, 'Cannot recognize given unit "'//TRIM(MODS(inm_pres)%UNIT)//'" for pressure. Exiting. '
              stop
          end if
      end if
  end if

  do i=1,N_VARS
      buf = UCASE(TRIM(MODS(i)%UNIT))
      if (TRIM(buf) /= '#' .and. i /= inm_pres .and. i /= inm_tempK) THEN
          IF (TRIM(buf) /= 'PPM' .and. TRIM(buf) /= 'PPB' .and. TRIM(buf) /= 'PPT' .and. TRIM(buf) /= 'PPQ') THEN
              print FMT_FAT0, 'Cannot recognize unit "'//TRIM(MODS(i)%UNIT)//'" for '//TRIM(MODS(i)%name)//'. Exiting. '
              stop
          ELSE
              if ((MODS(i)%MODE > 0) .or. (MODS(i)%col > 1)  .or. (ABS(MODS(i)%multi - 1d0)>1d-9) .or. (ABS(MODS(i)%shift)>1d-16)) THEN
                  print FMT_MSG, '- Converting '//TRIM(MODS(i)%name)//' from '//TRIM(MODS(i)%UNIT)
              end if
          END IF
      end if
  end do

END SUBROUTINE CONVERT_PRESSURE_AND_VALIDATE_UNITS


SUBROUTINE PARSE_PARTICLE_GRID(file, parvar)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN)        :: file   ! 1. row is diameter, 1. column is time
    TYPE(particle_grid), INTENT(INOUT)  :: parvar ! The particle grid to be allocated and filled with data in file
    integer                             :: rows, cols, ioi, I
    real(dp)                            :: fl_buff(2)
    LOGICAL                             :: extracolumn

    ALLOCATE(parvar%conc_modelbins(n_bins_par))

    open(8889, file=TRIM(file), IOSTAT=ioi)
    call handle_file_io(ioi, file, 'Exiting the program')

    rows = rowcount(8889)
    IF ('.CSV' == UCASE(file(len(TRIM(file))-3:len(TRIM(file))))) THEN
      cols = colcount(8889, separator=',')
    ELSE
      cols = colcount(8889)
    END IF

    read(8889,*) fl_buff

    REWIND(8889)
    IF (fl_buff(2) > 0d0) THEN
        allocate(parvar%conc_matrix(rows-1,cols-1))
        ALLOCATE(parvar%sections(cols-1))
        extracolumn = .false.
    ELSE
        allocate(parvar%conc_matrix(rows-1,cols-2))
        ALLOCATE(parvar%sections(cols-2))
        extracolumn = .true.
    END IF

    ALLOCATE(parvar%time(rows-1))
    IF (extracolumn) THEN
        read(8889,*) fl_buff(1), fl_buff(2), parvar%sections(:)
    ELSE
        read(8889,*) fl_buff(1), parvar%sections(:)
    END IF
    DO I=1,rows-1
        IF (extracolumn) THEN
            read(8889,*) parvar%time(I), fl_buff(1), parvar%conc_matrix(I,:)
        ELSE
            read(8889,*) parvar%time(I), parvar%conc_matrix(I,:)
        END IF
    END DO
    CLOSE(8889)

END SUBROUTINE PARSE_PARTICLE_GRID


SUBROUTINE PARSE_MULTIMODAL()
    IMPLICIT NONE
    real(dp) :: buffer(100)
    INTEGER  :: ioi, sze

    buffer = -9999d0
    read(mmodal_input, *, IOSTAT=ioi) buffer
    sze =  SIZE(PACK(buffer,buffer>0))
    if (sze>=3) THEN
        ALLOCATE(MMODES(sze))
        MMODES = PACK(buffer,buffer>0)
    ELSE
        N_MODAL = -1d0
    ENDIF
    IF (mmodal_input_inuse == 0) N_MODAL = -1d0

END SUBROUTINE PARSE_MULTIMODAL


SUBROUTINE PARSE_MULTIMODAL_EMS()
    IMPLICIT NONE
    integer  :: temp(297) = 0, i,jjj1,jjj2,jjj3, M

    do i=N_VARS-N_XTRS, N_VARS
      if (MODS(i)%ISPROVIDED) THEN
        jjj1 = index(MODS(i)%NAME, '_GMD') - 1
        jjj2 = index(MODS(i)%NAME, '_ln(s)') - 1
        jjj3 = index(MODS(i)%NAME, '_emission') - 1

        if (jjj1>0) then
          read(MODS(i)%NAME(jjj1-1:jjj1),*) M
          temp(M * 3 - 2) = i
          N_MODAL_EMS = N_MODAL_EMS + 1

        else if (jjj2>0) then
          read(MODS(i)%NAME(jjj2-1:jjj2),*) M
          temp(M * 3 - 1) = i
          N_MODAL_EMS = N_MODAL_EMS + 1

        else if (jjj3>0) then
          read(MODS(i)%NAME(jjj3-1:jjj3),*) M
          temp(M * 3) = i
          N_MODAL_EMS = N_MODAL_EMS + 1

        END IF
      END IF
    end do

    if (MOD(N_MODAL_EMS,3) /= 0 ) &
      STOP 'AEROSOL EMISSION INPUT IS INCOMPLETE, PROGRAM WILL EXIT.'
    N_MODAL_EMS = N_MODAL_EMS/3
    ALLOCATE(MMODES_EMS(SIZE(PACK(temp,temp>0))))
    MMODES_EMS = PACK(temp,temp>0)
END SUBROUTINE PARSE_MULTIMODAL_EMS


SUBROUTINE PARSE_GR_SIZES()
    IMPLICIT NONE
    real(dp) :: buffer(100), bin_ratio, bin2
    INTEGER  :: ioi, sze

    bin_ratio = (log(max_particle_diam/min_particle_diam))/n_bins_par
    ! Define second bin diameter
    bin2 = exp(bin_ratio + log(min_particle_diam))

    buffer = -9999d0
    read(GR_sizes, *, IOSTAT=ioi) buffer ! split user supplied string into string array
    ! see how many of the user supplied diameters fall within the PSD
    sze =  SIZE(    PACK(PACK(buffer,buffer>0),PACK(buffer,buffer>0)>bin2)    )
    if (sze>0) THEN
        ALLOCATE(GR_bins(sze))
        ALLOCATE(GGR(sze))
        GR_bins = PACK(PACK(buffer,buffer>0),PACK(buffer,buffer>0)>bin2)
    ELSE
        CALC_GR = .false.
    ENDIF

END SUBROUTINE PARSE_GR_SIZES


SUBROUTINE SW_PP()
    IMPLICIT NONE
    INTEGER :: cols,rows,ii,allrows,aa,bb,skiphdr
    REAL(dp) :: dummy, DL
    REAL(dp), ALLOCATABLE :: WL(:), Weight(:)
    character(1) :: bufr

    print FMT_HDR, ''

    ! Refuse to continue with default spectral file if it should represent actinic flux
    if (TRIM(spectrumfile) == "ModelLib/Photolyse/Spectra/glob_swr_distr.txt" .and. SWR_IS_ACTINICFLUX) &
        STOP 'Default spectrum is not valid for Actinic Flux'

    open(UNIT=685,FILE=TRIM(spectrumfile), STATUS='OLD', ACTION='READ', iostat=ii)
    CALL handle_file_io(ii, spectrumfile, 'stop')
    cols = COLCOUNT(685)

    allrows = ROWCOUNT(685,'%')
    rows = ROWCOUNT(685,'#')
    skiphdr = allrows-rows

    REWIND(685)

    if (rows==1.and.cols>2) THEN
        print FMT_MSG, 'Spectral data is time independent. Assuming that data starts at 280 nm, with dL = 5 nm.'
        print FMT_SUB, ' Single row vector -> FIRST COLUMN IS TIME AND OMITTED!'
        do ii=1,allrows
            if (ii<=skiphdr) THEN
                READ(685,*) bufr
            else
                READ(685,*) dummy, swr_spectrum(1:84)
            end if
        end do

    else if (rows>1.and.cols==1) THEN
        print FMT_MSG, 'Spectral data is time independent. Assuming that data starts at 280 nm, with dL = 5 nm.'
        do ii=1,min(allrows, 84+skiphdr)
            if (ii<=skiphdr) THEN
                READ(685,*) bufr
            else
                READ(685,*) swr_spectrum(ii-skiphdr)
            end if
        end do

    else if (rows>1.and.cols==2) THEN
        print FMT_MSG, 'Spectral data is time independent. Assuming that data starts at 280 nm, with dL = 5 nm.'
        print FMT_SUB, ' IGNORING WAVELENGTHS: FIRST COLUMN IS OMITTED!'
        do ii=1,min(allrows, 84+skiphdr)
            if (ii<=skiphdr) THEN
                READ(685,*) bufr
            else
                READ(685,*) dummy, swr_spectrum(ii-skiphdr)
            end if
        end do

    else if (rows>1.and.cols>2) THEN
        print FMT_MSG, 'Spectral data is time dependent. Assuming that data starts at 280 nm, with dL = 5 nm.'
        print FMT_SUB, ' MATRIX DATA: FIRST COLUMN IS TIME.'
        ALLOCATE(swr_temporal_data(rows,84))
        ALLOCATE(swr_times(rows))
        swr_is_time_dependent = .true.

        swr_times = 0d0
        swr_temporal_data = 0d0

        do ii=1,allrows
            if (ii<=skiphdr) THEN
                READ(685,*) bufr
            else
                READ(685,*) swr_times(ii-skiphdr), swr_temporal_data(ii-skiphdr,1:84)
            end if
        end do

    END IF


    if ( (TRIM(spectrumfile) == "ModelLib/Photolyse/Spectra/glob_swr_distr.txt") .and. &
        ((.not. equal(SWR_IN_LOWER, 300d0)).or.(.not. equal(SWR_IN_UPPER, 4000d0))) ) &
        THEN ! normalize
        REWIND(685)

        ALLOCATE(WL(rows))
        ALLOCATE(Weight(rows))
        do ii=1, rows
            READ(685, *)  WL(ii), Weight(ii)
            ! print* ,WL(ii), Weight(ii)
            ! if (WL(ii)<SWR_IN_LOWER) Weight(ii) = 0d0
            ! if (WL(ii)>SWR_IN_UPPER) Weight(ii) = 0d0
            ! print* ,WL(ii), Weight(ii)
        end do
        aa = minloc(abs(WL-SWR_IN_LOWER),1)
        bb = minloc(abs(WL-SWR_IN_upper),1)
        DL = WL(2) - WL(1)

        swr_spectrum = swr_spectrum / sum(Weight(aa:bb)*DL)
        print FMT_MSG, 'Normalized the spectrum with '//TRIM(f2chr( 1/sum(Weight(aa:bb)*DL)))

        DEALLOCATE(WL)
        DEALLOCATE(Weight)

    end if

    close(685)

    If (SWR_IS_ACTINICFLUX) print FMT_WARN0, 'TREATING SHORTWAVE DATA AS IT WERE ACTINIC FLUX!'

    print FMT_HDR, ''

END SUBROUTINE SW_PP


SUBROUTINE PARSE_ACDC_SYSTEMS

    IMPLICIT NONE
    INTEGER :: ii,jj,ioi=0, counter,mm
    CHARACTER(len=16)   :: name(24)
    CHARACTER(len=256)  :: System,Energies,Dipoles,path,Nickname
    NAMELIST /ACDC_RECORD_NML/ System,Energies,Dipoles,Nickname
    print *, ''
    print FMT_HDR, 'Allocating ACDC systems to the selected input'


    do jj=1,size(G_ACDC)

        name = '----------------'
        read(ACDC_links(jj),*, iostat=ioi) name(:)

        counter = 0
        do ii=0,11
            if (name((ii*2)+1)/= '----------------') THEN
                counter = counter + 1
            end if
        end do

        if (ACDC_SYSTEMS(jj)==1) THEN
            print FMT_MSG, 'ACDC submodule #'//i2chr(jj)//' is initialized with input definitions.'

            WRITE(path,'(a,i0.2,a)') 'src/ACDC/ACDC_',jj,'/ACDC_RECORD_NML'
            OPEN(UNIT=889, FILE=TRIM(path), STATUS='OLD', ACTION='READ', iostat=ii)
            if (ii==0) THEN
                Nickname = '-not defined-'
                do while (ii==0)
                    READ(889,NML = ACDC_RECORD_NML, IOSTAT=ii)
                end do
                close(889)
                G_ACDC(jj)%SYSTEM_FILE=System
                G_ACDC(jj)%ENERGY_FILE=Energies
                G_ACDC(jj)%DIPOLE_FILE=Dipoles
                G_ACDC(jj)%NICKNAME=Nickname
                print FMT_SUB, 'System name: '//TRIM(G_ACDC(jj)%NICKNAME)
                print FMT_SUB, 'System based on '//TRIM(G_ACDC(jj)%SYSTEM_FILE)
                print FMT_SUB, 'Energies from '//TRIM(G_ACDC(jj)%ENERGY_FILE)
                print FMT_SUB, 'Dipoles from '//TRIM(G_ACDC(jj)%DIPOLE_FILE)
            END IF


            G_ACDC(jj)%inuse = .true.

            ALLOCATE(G_ACDC(jj)%ACDC_LINK_IND(counter))
            ALLOCATE(G_ACDC(jj)%ACDC_monConc(counter))
            G_ACDC(jj)%ACDC_monConc(counter) = 0d0
            ALLOCATE(G_ACDC(jj)%ACDC_MONOMER_NAMES(counter))

            do ii=1,counter
                if (name(((ii-1)*2)+1)/= '----------------') THEN
                    G_ACDC(jj)%ACDC_MONOMER_NAMES(ii)(:) = TRIM(name(((ii-1)*2)+1))
                    G_ACDC(jj)%ACDC_LINK_IND(ii) = IndexFromName(TRIM(name(((ii-1)*2)+2)), [(MODS(mm)%NAME, mm=1,size(MODS))] )
                    IF (G_ACDC(jj)%ACDC_LINK_IND(ii) == 0) THEN
                        print FMT_WARN0, 'In ACDC '//i2chr(jj)//': Compound '//TRIM(name(((ii-1)*2)+2))//' does not exist in input'
                        print FMT_SUB, 'Searching compound '//TRIM(name(((ii-1)*2)+2))//' from chemistry...'
                        G_ACDC(jj)%ACDC_LINK_IND(ii) = -1 * IndexFromName(TRIM(name(((ii-1)*2)+2)), SPC_NAMES )
                        if (G_ACDC(jj)%ACDC_LINK_IND(ii)<0) print FMT_SUB, 'Found '//TRIM(name(((ii-1)*2)+2))//' from chemistry. Using as input '//i2chr(ii)//' for ACDC '//i2chr(jj)//'.'
                        if (G_ACDC(jj)%ACDC_LINK_IND(ii)==0) THEN
                            print FMT_FAT0, 'In ACDC '//i2chr(jj)//': Could not find '//TRIM(name(((ii-1)*2)+2))
                            stop
                        END IF
                    END IF
                end if
            end do
        ELSE
            print FMT_MSG, 'ACDC submodule #'//i2chr(jj)//' is not used.'
        END IF

    end do

    print FMT_LEND,


END SUBROUTINE PARSE_ACDC_SYSTEMS

SUBROUTINE PARSE_INIT_ONLY()
  implicit none
  CHARACTER(25) :: buffer = ''
  CHARACTER(25) :: bufferlist(100) = '#'
  INTEGER :: ii

  iF (INIT_ONLY == '') RETURN

  ALLOCATE(init_only_these(0))
  read(INIT_ONLY, *,IOSTAT=ii) bufferlist(:)
  do ii=1,100
    if (TRIM(bufferlist(ii)) == '#') EXIT
    if (IndexFromName(TRIM(bufferlist(ii)),SPC_NAMES) > 0) &
      init_only_these = [init_only_these,IndexFromName(TRIM(bufferlist(ii)),SPC_NAMES)]
  end do

END SUBROUTINE PARSE_INIT_ONLY

end module INPUT
