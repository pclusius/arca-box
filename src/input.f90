Module INPUT
! Module to read init file
use second_Monitor
USE second_precision, ONLY: dp
use constants
USE auxillaries

Implicit none

! Relative path to NAMES.DAT
CHARACTER(27), PARAMETER :: NAMESDAT = 'ModelLib/required/NAMES.dat'

INTEGER :: N_VARS ! This will store the number of variables in NAMES.DAT
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

INTEGER, ALLOCATABLE :: INDRELAY_CH(:)
INTEGER, ALLOCATABLE :: INDRELAY_TIED(:)
INTEGER, ALLOCATABLE :: index_cond(:)

REAL(dp), allocatable, private :: INPUT_ENV(:,:)  ! will be of same shape as the files
REAL(dp), allocatable, private :: INPUT_MCM(:,:)  ! will be of same shape as the files
REAL(dp), allocatable :: timevec(:)     ! Whatever the times were for (currently) ALL measurements
REAL(dp), allocatable :: CONC_MAT(:,:)  ! will be of shape ( len(timevec) : N_VARS )
real(dp), allocatable :: par_data(:,:)
real(dp), allocatable :: GGR(:)

! variable for storing init file name
CHARACTER(len=256) :: Fname_init ! init file names
CHARACTER(len=5)   :: gui        ! magic word for gui in use
LOGICAL            :: ingui = .False. ! True if program is invoked from gui, from command line
CHARACTER(len=14), parameter  :: namelists(21) = &
['NML_TIME      ','NML_Flag      ','NML_Path      ','NML_MISC      ','NML_VAP       ','NML_PARTICLE  ','NML_ENV       ',&
 'NML_MCM       ','NML_MODS      ','NML_PRECISION ','NML_CUSTOM    ','              ','              ','              ',&       ! to help troubleshooting
 '              ','              ','              ','              ','              ','              ','              ']       ! to help troubleshooting

! MAIN PATHS
! CHARACTER(len=256):: WORK_DIR   = ''
CHARACTER(len=256):: INOUT_DIR   = 'INOUT'
CHARACTER(len=256):: CASE_NAME  = 'DEFAULTCASE'
CHARACTER(len=30) :: RUN_NAME   = 'DEFAULTRUN'
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
Logical :: Use_speed           = .false.
Logical :: AFTER_CHEM_ON       = .false.
Logical :: AFTER_NUCL_ON       = .false.
! Logical :: INIT_W_MODAL        = .true.

NAMELIST /NML_Flag/ chemistry_flag, Aerosol_flag, ACDC_solve_ss, ACDC, & !NUCLEATION,
         Condensation, Coagulation, Deposition, Chem_Deposition, model_H2SO4, RESOLVE_BASE, &
         PRINT_ACDC, Use_speed, ORG_NUCL, AFTER_CHEM_ON, AFTER_NUCL_ON !,INIT_W_MODAL, Extra_data

! Logical :: USE_OPENMP   = .false.
! NAMELIST /NML_PARALLEL/ USE_OPENMP

! TIME OPTIONS
real(dp)  :: runtime = 1d0
real(dp)  :: FSAVE_INTERVAL = 300d0
real(dp)  :: PRINT_INTERVAL = 15*60d0
INTEGER   :: FSAVE_DIVISION = 0
real(dp)  :: dt = -1d0
CHARACTER(len=10)  :: DATE = '1800-01-01', INDEX = ''
NAMELIST /NML_TIME/ runtime, dt, FSAVE_INTERVAL, PRINT_INTERVAL, FSAVE_DIVISION, DATE, INDEX

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
REAL(dp)            :: N_MODAL = -1d0   !for use_dmps_partial, read all dmps data below this diameter [m]
NAMELIST /NML_PARTICLE/ PSD_MODE,n_bins_par,min_particle_diam,max_particle_diam, DMPS_file,extra_particles,& !DMPS_dir,extra_p_dir,
                        DMPS_read_in_time,dmps_highband_lower_limit, dmps_lowband_upper_limit,use_dmps,use_dmps_partial, &
                        mmodal_input, N_MODAL, mmodal_input_inuse

REAL(dp), ALLOCATABLE            :: MMODES(:)  !for use_dmps_partial, read all dmps data below this diameter [m]
type(particle_grid), ALLOCATABLE :: xtras(:)
! BG_PAR is here in case you want to use it Carlton and Lukas, in the end we remove either BG_PAR or par_data mmkay.
type(particle_grid) :: BG_PAR       ! Var to store the particle size distribution. This might become redundant
type(particle_grid) :: PAR_LOSSES   ! Var to store losses file, which could be either one row or a matrix, but has time
                                    ! and size dependency


! ENVIRONMENTAL INPUT
CHARACTER(len=256)  :: ENV_FILE = ''
CHARACTER(len=256)  :: LOSSES_FILE = ''
REAL(dp)            :: CHAMBER_FLOOR_AREA = 0d0
REAL(dp)            :: CHAMBER_CIRCUMFENCE = 0d0
REAL(dp)            :: CHAMBER_HEIGHT = 0d0
REAL(dp)            :: EDDYK = 5d-2
REAL(dp)            :: ustar = 5d-2
REAL(dp)            :: ALPHAWALL = 5d-5
NAMELIST /NML_ENV/ ENV_file, LOSSES_FILE, CHAMBER_FLOOR_AREA, CHAMBER_CIRCUMFENCE, CHAMBER_HEIGHT,EDDYK, ustar, ALPHAWALL

! MCM INPUT
! CHARACTER(len=256)  :: MCM_path = ''
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
CHARACTER(1000)  :: Description='*'
CHARACTER(len=256)      :: GR_sizes = '6d-9 20d-9 30d-9'
NAMELIST /NML_MISC/ lat, lon, wait_for, Description, CH_Albedo, DMA_f, resolve_BASE_precision, Fill_formation_with, skip_acdc, &
                    GR_sizes

Logical                 :: VAP_logical = .True.
Logical                 :: Use_atoms = .False.
CHARACTER(len=256)      :: Vap_names
CHARACTER(len=256)      :: Vap_atoms = ''

NAMELIST /NML_VAP/ VAP_logical, Use_atoms, Vap_names, Vap_atoms !, Vap_props

CHARACTER(1000) :: INITIALIZE_WITH = ''
INTEGER  :: limit_vapours = 999999
INTEGER  :: acdc_iterations = 1
INTEGER  :: INITIALIZE_FROM = 0
Logical  :: use_raoult = .True.
Logical  :: variable_density = .False.
! if true, will not save condensible vapour concentration in Particles.nc. They will always be saved also in Chemistry.nc
Logical  :: DONT_SAVE_CONDENSIBLES = .False.
! If True, skips ACDC with very low concentrations and negligible formation rates
real(dp) :: dmps_tres_min = 10.
real(dp) :: VP_MULTI = 1d0
real(dp) :: start_time_s = 0d0
real(dp) :: END_DMPS_SPECIAL = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: FLOAT_CHEMISTRY_AFTER_HRS = 1d100 ! Any number larger than runtime will do as defaults
real(dp) :: dmps_multi       = 1d6 ! Multiplicator to convert dmps linear concentration to #/m^3
Logical  :: NO2_IS_NOX        = .false.
Logical  :: USE_RH_CORRECTION = .true.
LOGICAL  :: TEMP_DEP_SURFACE_TENSION = .False.
LOGICAL  :: use_diff_dia_from_diff_vol = .False.
REAL(dp), ALLOCATABLE   :: GR_bins(:)  ! used for GR calculation [m]
Logical                 :: CALC_GR = .True.
real(dp)                :: speed_dt_limit(8) =  [1d-5,60d0,&
                                                1d-5,60d0,&
                                                1d-5,60d0,&
                                                1d-5,60d0]

! defined in Constants: Logical  :: NO_NEGATIVE_CONCENTRATIONS = .false.

NAMELIST /NML_CUSTOM/ use_raoult, acdc_iterations,variable_density,dmps_tres_min, &
                      start_time_s, dmps_multi, INITIALIZE_WITH,INITIALIZE_FROM, VP_MULTI, &
                      DONT_SAVE_CONDENSIBLES, limit_vapours, END_DMPS_SPECIAL,NO2_IS_NOX,&
                      NO_NEGATIVE_CONCENTRATIONS, FLOAT_CHEMISTRY_AFTER_HRS, USE_RH_CORRECTION, &
                      TEMP_DEP_SURFACE_TENSION, use_diff_dia_from_diff_vol, speed_dt_limit

! ==================================================================================================================
! Define change range in percentage
REAL(dp), DIMENSION(2) :: diameter_prec_def = [1.d-1, 2d0] ! -> minimum and maximum relative change in particle diameter
REAL(dp), DIMENSION(2) :: pnumber_prec_def  = [1.d-1, 2d0] ! -> minimum and maximum relative change in particle number
REAL(dp), DIMENSION(2) :: vapour_prec_def   = [1.d-1, 2d0] ! -> minimum and maximum relative change in particle concentration
CHARACTER(22)          :: range_names(3)    = ['particle diameter     ',&
                                               'particle concentration',&
                                               'vapour concentration  ']
! Defines the minimum/maximum relative change caused by a process within a timestep
REAL(dp), DIMENSION(3,2) :: change_range = 0
NAMELIST /NML_PRECISION/ change_range, diameter_prec_def,pnumber_prec_def,vapour_prec_def

! Options for screen output
LOGICAL :: clusterfractions,jions,timestep_multipliers,time_efficiency,Jorganic,GR
NAMELIST /NML_SCREENPRINTS/ clusterfractions,jions,timestep_multipliers,time_efficiency,Jorganic,GR

contains

! ======================================================================================================================
! Subroutine reads all user input from INITFILE and input files and processes the input etc.
! ......................................................................................................................
subroutine READ_INPUT_DATA()
    IMPLICIT NONE
    CHARACTER(len=256)                :: buf
    integer, allocatable              :: Natoms(:,:)
    integer                           :: ioi, ii
    integer                           :: i, j, k, jj, path_l(2), N_Xtr = 0
    integer                           :: rows, cols
    LOGICAL                           :: elements_missing = .false.
    real(dp)                          :: molar_mass, parameter_A, parameter_B
    CHARACTER(len=256)                :: species_name
    CHARACTER(len=20), ALLOCATABLE    :: atoms_name(:)

    ! CHECK HOW MANY POSSIBLE INPUT VARIABLES (METEOROLOGICAL, MCM ETC.) THERE ARE IN THE MODEL
    OPEN(800, file=NAMESDAT, ACTION='READ', status='OLD', iostat=ioi)
    call handle_file_io(ioi, NAMESDAT, 'This file is essential and should not be changed.')
    N_VARS = rowcount(800)
    CLOSE(800)

    ! BASED ON N_VARS, ALLOCATE VECTORS
    ALLOCATE(MODS(N_VARS))
    ALLOCATE(INDRELAY_CH(N_VARS))
    ALLOCATE(INDRELAY_TIED(N_VARS))
    INDRELAY_CH = 0
    INDRELAY_TIED = 0

    CALL NAME_MODS_SORT_NAMED_INDICES
    CALL READ_INIT_FILE
    CALL PUT_USER_SUPPLIED_TIMEOPTIONS_IN_GTIME
    CALL REPORT_INPUT_COLUMNS_TO_USER

    ! set up number of ACDC iterations.
    if (ACDC_solve_ss) acdc_iterations = 1

    ! ALLOCATE CONC_MAT Currently both files need to have same time resolution FIX THIS SOON!
    ! The idea here is to count the rows to get time and allocate CONCMAT and TIMEVEC
    ! This is very much under construction
    IF ((ENV_file /= '') .or. (MCM_file /= '' .and. Chemistry_flag)) THEN
        IF (ENV_file /= '') THEN
            OPEN(unit=801, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
            CALL handle_file_io(ioi, ENV_file, 'stop')

        ELSE
            OPEN(unit=801, File=TRIM(MCM_file), ACTION='READ', STATUS='OLD',iostat=ioi)
            CALL handle_file_io(ioi, MCM_file, 'stop')

        END IF
        ALLOCATE(CONC_MAT(ROWCOUNT(801,'#'),N_VARS))
        ALLOCATE(TIMEVEC(ROWCOUNT(801,'#')))
        CLOSE(801)
    ! Deal with a situation where we have no input. We still need conc_mat and timevec.
    ELSE
        ALLOCATE(CONC_MAT(2,N_VARS))
        ALLOCATE(TIMEVEC(2))
        TIMEVEC = (/0d0, GTIME%SIM_TIME_H/)
    END IF

    CONC_MAT = 0d0

    ! READ ENV INPUT
    if (ENV_file /= '') THEN
        OPEN(unit=801, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
        CALL handle_file_io(ioi, ENV_file, 'Terminating on CONCMAT allocation')

        rows = ROWCOUNT(801,'#')
        cols = COLCOUNT(801)

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
        cols = COLCOUNT(801)
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

IF ((TRIM(LOSSES_FILE) /= '') .and. Aerosol_flag) CALL PARSE_PARTICLE_GRID(LOSSES_FILE, PAR_LOSSES)

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

IF (Aerosol_flag) then

    CALL PARSE_MULTIMODAL
    CALL PARSE_GR_SIZES

    write(*,FMT_MSG) 'Reading Vapour name file '// TRIM(Vap_names)
    OPEN(unit=802, File= TRIM(Vap_names) , STATUS='OLD', iostat=ioi)
    call handle_file_io(ioi, Vap_names, &
        'If Condensation is used, "Vapour file" must be defined (in tab "Advanced").')

    rows = ROWCOUNT(802, '#')
    cols = COLCOUNT(802)

    VAPOUR_PROP%n_condorg = 0
    do j = 1,rows
        read(802,*,iostat=ioi) species_name
        if (j<=limit_vapours .or. j==rows) THEN
            k = IndexFromName( TRIM(species_name), SPC_NAMES )
            if (k>0) THEN
                VAPOUR_PROP%n_condorg = VAPOUR_PROP%n_condorg + 1
            end if
        end if
    end do
    REWIND(802)

    VAPOUR_PROP%n_condtot = VAPOUR_PROP%n_condorg + 1
    allocate(VAPOUR_PROP%Vapour_names (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%molar_mass   (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%psat_a       (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%psat_b       (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%molec_mass   (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%molec_volume (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%density      (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%surf_tension (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%diff         (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%c_speed      (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%c_sat        (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%cond_type    (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%molec_dia    (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%mfractions   (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%alpha        (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%diff_vol     (VAPOUR_PROP%n_condtot) )
    allocate(VAPOUR_PROP%diff_dia     (VAPOUR_PROP%n_condtot) )
    ALLOCATE(index_cond(VAPOUR_PROP%n_condorg))
    index_cond = 0

    if (USE_RH_CORRECTION) THEN
        ! These are only allocated for sulfuric acid, maybe nitric acid in the future.
        ! Uses same index as in VAPOUR_PROP%ind_H2SO4
        allocate(VAPOUR_PROP%wet_dia(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_condtot))
        allocate(VAPOUR_PROP%wet_mass(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_condtot))
    END IF


    print FMT_SUB, 'Compounds picked from Vapours file: '//TRIM(i2chr(VAPOUR_PROP%n_condorg))
    print FMT_SUB, 'Total number of condensibes       : '//TRIM(i2chr(VAPOUR_PROP%n_condtot))

    !!! reading the vap names and vap vapour_properties
    VAPOUR_PROP%ind_GENERIC = VAPOUR_PROP%n_condorg
    VAPOUR_PROP%ind_H2SO4   = VAPOUR_PROP%n_condtot
    VAPOUR_PROP%Mfractions  = 0.0
    VAPOUR_PROP%Mfractions(VAPOUR_PROP%n_condorg) = 1d0 !

    ! ---------------------------------------------------------------------
    ! ORGANIC VAPOUR PROPERTIES
    ! ---------------------------------------------------------------------
    ii = 1
    do j = 1, rows
        read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B
        if (j<=limit_vapours .or. j==rows) THEN ! j==rows is because "GENERIC" is selected even if limitvapours<rows

            ! Check if the compounds exists in Chemistry and only then add to vapours
            k = IndexFromName( species_name, SPC_NAMES )

            if (k>0) THEN
                ! fill the hash table for vapour index -> chemistry index
                index_cond(ii) = k

                VAPOUR_PROP%molar_mass(ii)   = molar_mass *1D-3 ! kg/mol
                VAPOUR_PROP%psat_a(ii)       = parameter_A
                VAPOUR_PROP%psat_b(ii)       = parameter_B
                VAPOUR_PROP%vapour_names(ii) = TRIM(species_name)
                VAPOUR_PROP%molec_mass(ii)   = VAPOUR_PROP%molar_mass(ii)/Na  !kg/#

                ! Option for simple parametrisation of organic vapour liquid density. Use with caution, not yet thoroughly implemented
                IF (variable_density) THEN
                    VAPOUR_PROP%density(ii)  = -30d0 * (parameter_A - parameter_B/293.15d0) + 1029d0  ! kg/m3
                ELSE
                    VAPOUR_PROP%density(ii)  = 1400.0  ! kg/m3
                END IF
                VAPOUR_PROP%molec_volume(ii) = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
                VAPOUR_PROP%diff_vol(ii)     = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
                VAPOUR_PROP%surf_tension(ii) = 0.05
                VAPOUR_PROP%cond_type(ii)    = 1  ! not an acid (H2SO4 or HCL)
                VAPOUR_PROP%alpha(ii)        = 1.0
                VAPOUR_PROP%c_sat(ii)        = saturation_conc_m3(VAPOUR_PROP%psat_a(ii),VAPOUR_PROP%psat_b(ii), 293.15d0)

                ii = ii + 1
            END IF
        END IF
    end do

    close(802)

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

        DO j=1,VAPOUR_PROP%n_condorg
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

    ! Now the volumes are updated, the diameter can be calculated
    VAPOUR_PROP%molec_dia = (6D0 * VAPOUR_PROP%molec_volume / pi )**(1D0/3D0)  ! molecular diameter [m]
    if (use_diff_dia_from_diff_vol) THEN
        VAPOUR_PROP%diff_dia = (6D0 * 1d-30 * VAPOUR_PROP%diff_vol / pi )**(1D0/3D0)  ! molecular diameter [m]
    ELSE
        VAPOUR_PROP%diff_dia = VAPOUR_PROP%molec_dia  ! molecular diameter [m]
    END IF
    ! ---------------------------------------------------------------------
    ! Sulfuric acid treated separately
    ! ---------------------------------------------------------------------
    ii = VAPOUR_PROP%n_condtot
    VAPOUR_PROP%vapour_names(ii)  = 'H2SO4'
    VAPOUR_PROP%cond_type(ii)     = 2  ! Acid
    VAPOUR_PROP%molar_mass(ii)    = 98.0785 * 1d-3
    VAPOUR_PROP%psat_a(ii)        = 0
    VAPOUR_PROP%psat_b(ii)        = 20000d0
    VAPOUR_PROP%density(ii)       = 1819.3946 ! kg/m3
    VAPOUR_PROP%molec_mass(ii)    = VAPOUR_PROP%molar_mass(ii)/Na
    VAPOUR_PROP%molec_volume(ii)  = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
    VAPOUR_PROP%diff_vol(ii)      = 4*6.11D0 + 2*2.31D0 + 22.9D0 ! O=4, H=2, S=1
    VAPOUR_PROP%diff_dia(ii)      = VAPOUR_PROP%molec_dia(ii)
    VAPOUR_PROP%surf_tension(ii)  = 0.07
    VAPOUR_PROP%molec_dia(ii)     = (6D0 * VAPOUR_PROP%molec_volume(ii) / pi )**(1D0/3D0)  ! molecular diameter [m]
    VAPOUR_PROP%alpha(ii)         = 1.0
    VAPOUR_PROP%c_sat(ii)         = 0.0 ! Sulfuric acid stays put

  end if

  CALL CHECK_MODIFIERS ! Print out which modifiers differ from default values
end subroutine READ_INPUT_DATA


subroutine READ_INIT_FILE
!-------------------------------------------------------------------------------
! Reads the init file and fills user-provided variables
!-------------------------------------------------------------------------------
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
  write(*,FMT_HDR) 'READING USER DEFINED INTIAL VALUES FROM: '//TRIM(ADJUSTL(Fname_init))

  do k=1, ROWCOUNT(888); READ(888,NML = NML_TIME, IOSTAT=IOS(i)) ! #1
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_Flag, IOSTAT=IOS(i)) ! #2
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_Path, IOSTAT=IOS(i)) ! #3
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MISC, IOSTAT=IOS(i)) ! #4
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_VAP, IOSTAT=IOS(i)) ! #5
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_PARTICLE, IOSTAT=IOS(i)) ! #6
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_ENV, IOSTAT=IOS(i)) ! #7
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MCM, IOSTAT=IOS(i)) ! #8
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_MODS, IOSTAT=IOS(i)) ! #9
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_PRECISION, IOSTAT=IOS(i)) ! #10
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  do k=1, ROWCOUNT(888); READ(888,NML = NML_CUSTOM, IOSTAT=IOS(i)) ! #11
  IF (IOS(i) == 0) EXIT;end do; REWIND(888); i=i+1

  CLOSE(888)


  IF (SUM(ABS(IOS)) /= 0) then
    write(*,FMT_MSG) 'There was a problem with INITFILE. Check the file and refer to manual. Exiting now.'
    DO i=1,size(IOS,1)
        write(*,FMT_SUB) TRIM(namelists(i))//' ('//TRIM(i2chr(I))//') returned '//TRIM(i2chr(IOS(i)))
    END DO
    write(*,FMT_LEND)
    STOP
  end if
  IF (INOUT_DIR(len(TRIM(INOUT_DIR)):len(TRIM(INOUT_DIR))) == '/') INOUT_DIR(len(TRIM(INOUT_DIR)):len(TRIM(INOUT_DIR))) = ' '
  ! Also save all settings to initfile. Use this file to rerun if necessary
  open(889, file=TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(INDEX)//'/'//TRIM(RUN_NAME)//'/NMLS.conf', action='WRITE')
  write(889,NML = NML_TIME       ) ! directories and test cases
  write(889,NML = NML_Flag       ) ! flags
  write(889,NML = NML_Path       ) ! time related stuff
  write(889,NML = NML_MISC       ) ! dmps_file information
  write(889,NML = NML_VAP        ) ! environmental information
  write(889,NML = NML_PARTICLE   ) ! MCM_file information
  write(889,NML = NML_ENV        ) ! modification parameters
  write(889,NML = NML_MCM        ) ! misc input
  write(889,NML = NML_CUSTOM     ) ! custom input
  write(889,NML = NML_MODS       ) ! vapour input
  close(889)

  DO i=1, N_VARS
    IF (TRIM(MODS(i)%UNIT) == '-') THEN
        MODS(i)%UNIT = '#'
    ELSE
        MODS(i)%ISPROVIDED = .true.
    END IF
  END DO

end subroutine READ_INIT_FILE

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
        GTIME%FSAVE_INTERVAL = INT(GTIME%SIM_TIME_S/GTIME%dt) / FSAVE_DIVISION * GTIME%dt
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
        if (INDEX /= '') THEN
            print FMT_HDR,
            print FMT_HDR, 'USING INDEX INSTEAD OF DATE -> ASSUMING LIGHT DIRECTION IS FROM DIRECTLY UP'
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
    INTEGER :: i

    OPEN(800, file=NAMESDAT, ACTION='READ', status='OLD')
    DO i = 1,N_VARS
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
            IF ((TRIM(ENV_file) == '') .and. (I<LENV)) THEN
                print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(i2chr(MODS(I)%col))//' but no ENV_FILE. Using parametrisation?'
            ELSE IF ((TRIM(MCM_file) == '') .and. (I>LENV)) THEN
                print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(i2chr(MODS(I)%col))//' but no MCM_FILE. Using parametrisation?'
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
        IF ((ioi /= 0) .and. (i==1)) THEN
            REWIND(unit)
            READ(unit,*, iostat=ioi) dump
            print FMT_SUB, 'First row omitted from file "'// TRIM(Input_file) //'".'
        ELSE IF ((ioi /= 0) .and. (i>1)) THEN
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

    if (TRIM(UCASE(MODS(inm_TempK)%UNIT)) == '#') THEN
        print FMT_WARN0, "No unit for temperature. Use either 'K' or 'C'. Now assuming Kelvins. This may lead to SIGFPE."
        TempUnit = 'K'
    ELSEIF (TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'K' .or. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'C' )THEN
        TempUnit = TRIM(UCASE(MODS(inm_TempK)%UNIT))
    END IF

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
    cols = colcount(8889)

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

end module INPUT
