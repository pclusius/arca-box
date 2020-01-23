Module INPUT
! Module to read init file

USE second_precision, ONLY: dp
use constants
USE auxillaries
USE Aerosol_auxillaries

Implicit none

! public :: vapours, ambient
INTEGER :: N_VARS ! This will store the number of variables in NAMES.DAT
INTEGER :: LENV   ! This will store the number of named indices in this code

! INDICES TO PROPERLY COMBINE INPUT TO CORRECT VALUES IN THE MODEL
! ALL VARIABLES THAT ARE TIME DEPENDENT MUST BE HERE
! IF YOU ADD VARIABLES HERE, YOU NEED TO UPDATE:
! - this list
! - subroutine NAME_MODS_SORT_NAMED_INDICES
!------------------------------------------------------------------
INTEGER :: inm_TempK
INTEGER :: inm_pres
INTEGER :: inm_RH
INTEGER :: inm_CS
INTEGER :: inm_CS_NA
INTEGER :: inm_swr
INTEGER :: inm_IPR
INTEGER :: inm_H2SO4
INTEGER :: inm_NH3
INTEGER :: inm_DMA
INTEGER :: inm_SO2
INTEGER :: inm_NO
INTEGER :: inm_NO2
INTEGER :: inm_CO
INTEGER :: inm_H2
INTEGER :: inm_O3
INTEGER :: inm_JIN
INTEGER :: inm_JDMA

INTEGER, ALLOCATABLE :: INDRELAY_CH(:)

REAL(dp), allocatable, private :: INPUT_ENV(:,:)  ! will be of same shape as the files
REAL(dp), allocatable, private :: INPUT_MCM(:,:)  ! will be of same shape as the files
REAL(dp), allocatable :: timevec(:)     ! Whatever the times were for (currently) ALL measurements
REAL(dp), allocatable :: CONC_MAT(:,:)  ! will be of shape ( len(timevec) : N_VARS )
real(dp), allocatable :: par_data(:,:)

! variable for storing init file name
character(len=256), private :: Fname_init ! init file names

! MAIN PATHS
character(len=256):: WORK_DIR   = ''
character(len=256):: INOUT_DIR   = 'INOUT'
character(len=256):: CASE_NAME  = 'DEFAULTCASE'
character(len=30) :: RUN_NAME   = 'DEFAULTRUN'
NAMELIST /NML_Path/ Work_dir, INOUT_DIR, Case_name, RUN_NAME

! MODULES IN USE OPTIONS
Logical :: Chemistry_flag      = .false.
Logical :: Aerosol_flag        = .false.
Logical :: ACDC_solve_ss       = .false.
Logical :: NUCLEATION          = .false.
Logical :: ACDC                = .true.
Logical :: model_H2SO4         = .false.
Logical :: Condensation        = .false.
Logical :: Coagulation         = .false.
Logical :: Extra_data          = .false.
Logical :: Current_case        = .false.
Logical :: RESOLVE_BASE        = .false.
Logical :: PRINT_ACDC          = .false.
NAMELIST /NML_Flag/ chemistry_flag, Aerosol_flag, ACDC_solve_ss, NUCLEATION, ACDC, &
         Extra_data, Current_case, Condensation, Coagulation, model_H2SO4, RESOLVE_BASE, PRINT_ACDC

! TIME OPTIONS
real(dp)  :: runtime = 1d0
real(dp)  :: FSAVE_INTERVAL = 300d0
real(dp)  :: PRINT_INTERVAL = 15*60d0
INTEGER   :: FSAVE_DIVISION = 0
INTEGER   :: dt = -1
character(len=10)  :: DATE = '1800-01-01', INDEX = ''
NAMELIST /NML_TIME/ runtime, dt, FSAVE_INTERVAL, PRINT_INTERVAL, FSAVE_DIVISION, DATE, INDEX

! MODIFIER OPTIONS
! MODS is declared in CONSTANTS.f90, in order to be more widely available
NAMELIST /NML_MODS/ MODS

! ----------------------------------------------------------------------------------------------------------------------
! Particle related variables
INTEGER             :: PSD_MODE = 1
! PSD representation used
! 0 = Lukas basic
! 1 = Lukas advanced
! 2 = Lukas premium
INTEGER             :: n_bins_particle = 100    ! number of bins in the particle range
REAL(dp)            :: min_particle_diam = 1d-9 ! lower limit of particle range [m]
REAL(dp)            :: max_particle_diam = 1d-6 ! upper limit of particle range [m]

! DMPS INPUT
REAL(dp)            :: DMPS_read_in_time = 0d0

! 6) PSD input as discussed: real (time, total conc in first two columns, diameter in first row and concentration matrix) (nr_times + 1, nr_channels + 2)
! character(len=256)  :: DMPS_dir
! character(len=256)  :: extra_p_dir
character(len=256)  :: DMPS_file

! 4) name list of species making up the particle phase (we think that this will only be a few species that are measured in the particle phase): integer
! 5) number of nonvolatile species considered: integer
! 7) Composition of the particles: real (time: first column, diameter: first row, species mass fraction matrix [1]) (nr_times + 1, nr_channels)
INTEGER             :: n_xpar_options = 3    ! number of options in the extra_particles file, not including path and name
! The list is needed for every species in the read in particle phase (see "4) name list of species ...") or we make one long list
! (nr_times + 1, nr_bins * particle phase species)
! Note: the number of channels does not have to be the number of size bins
character(len=256)  :: extra_particles = '' ! file containing paths to extra particle sumfile
REAL(dp)            :: dmps_highband_lower_limit = 15d-9    !for use_dmps_special, read all dmps data above this diameter [m]
REAL(dp)            :: dmps_lowband_upper_limit = 6.*1d-10  !for use_dmps_special, read all dmps data below this diameter [m]
logical             :: use_dmps = .false.
logical             :: use_dmps_special = .false.
NAMELIST /NML_PARTICLE/ PSD_MODE,n_bins_particle,min_particle_diam,max_particle_diam, DMPS_file,extra_particles,& !DMPS_dir,extra_p_dir,
                        DMPS_read_in_time,dmps_highband_lower_limit, dmps_lowband_upper_limit,use_dmps,use_dmps_special

type(inert_particles), ALLOCATABLE :: xtras(:)
! BG_PAR is here in case you want to use it Carlton and Lukas, in the end we remove either BG_PAR or par_data mmkay.
type(inert_particles) :: BG_PAR


! ENVIRONMENTAL INPUT
character(len=256)  :: ENV_path = ''
character(len=256)  :: ENV_file = ''
character(len=1)  :: TempUnit = '' ! K or C
NAMELIST /NML_ENV/ ENV_path, ENV_file, TempUnit

! MCM INPUT
character(len=256)  :: MCM_path = ''
character(len=256)  :: MCM_file = ''
NAMELIST /NML_MCM / MCM_path, MCM_file

! MISC OPTIONS
real(dp)  :: lat              ! Latitude for Photochemistry
real(dp)  :: lon              ! Longitude for Photochemistry
real(dp)  :: CH_Albedo = 2d-1 ! Albedo
real(dp)  :: DMA_f = 0
real(dp)  :: resolve_BASE_precision = 1d-2
CHARACTER(3) :: Fill_formation_with = ''
INTEGER   :: JD = -1
INTEGER   :: wait_for = 0 ! -1 for no pause, 0 for indefinite and positive value for fixed amount of seconds
LOGICAL   :: python   = .false. ! 1 = the program will offer to plot the output file GENERAL
CHARACTER(1000)  :: Description='*'
CHARACTER(100)   :: Solver = ''
NAMELIST /NML_MISC/ JD, lat, lon, wait_for,python, Description,Solver, CH_Albedo, DMA_f, resolve_BASE_precision, Fill_formation_with

Logical  :: VAP_logical = .False.
Logical  :: Use_atoms = .False.
character(len=256)  :: Vap_names
character(len=256)  :: Vap_props
character(len=256)  :: Vap_atoms
NAMELIST /NML_VAP/ VAP_logical, Use_atoms, Vap_names, Vap_props, Vap_atoms

type(vapour_ambient)  :: vapours
type(ambient_properties)  :: ambient
type(atoms):: Natoms  ! atoms of hydrogen, oxygen, nitrogen and carbon. Used for calculating diffusion

real(dp),allocatable :: Diff_org(:), Vol_org(:)

contains

subroutine READ_INPUT_DATA()
  IMPLICIT NONE
  character(len=256)  :: buf
  type(nrowcol)       :: rowcol_count

  integer             :: ioi,ioi2, ii, iosp, ioprop, ioi3
  integer             :: i, j, k, z, xp, yp, path_l(2)
  !!! for vapour FILES
  character(len=256)  :: species_name
  real(dp)            :: molar_mass, parameter_A, parameter_B, fl_buff(2)

  ! Welcoming message
  print'(a,t35,a)', achar(10),  '--~:| HLS-BOX v.0.1 |:~--'//achar(10)

  ! CHECK HOW MANY POSSIBLE INPUT VARIABLES (METEOROLOGICAL, MCM ETC.) THERE ARE IN THE MODEL
  OPEN(2151, file='src/NAMES.dat', ACTION='READ', status='OLD', iostat=ioi)
  IF (ioi /= 0) THEN
    print FMT_FAT0, 'Could not open NAMES.DAT. This is a mandatory file and should be in src/ directory.'
    STOP
  END IF
  N_VARS = rowcount(2151)
  CLOSE(2151)

  ! BASED ON N_VARS, ALLOCATE VECTORS
  ALLOCATE(MODS(N_VARS))

  ALLOCATE(INDRELAY_CH(N_VARS))
  INDRELAY_CH = 0
  CALL NAME_MODS_SORT_NAMED_INDICES
  CALL READ_INIT_FILE
  CALL PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME
  CALL REPORT_INPUT_COLUMNS_TO_USER

  ! ALLOCATE CONC_MAT Currently both files need to have same time resolution FIX THIS SOON!
  ! The idea here is to count the rows to get time and allocate CONCMAT and TIMEVEC
  ! This is very much under construction
  IF ((ENV_file /= '') .or. (MCM_file /= '')) THEN
    IF (ENV_file /= '') THEN
      OPEN(unit=51, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
      CALL handle_file_io(ioi, ENV_file, 'stop')

    ELSE
      OPEN(unit=51, File=TRIM(MCM_file), ACTION='READ', STATUS='OLD',iostat=ioi)
      CALL handle_file_io(ioi, MCM_file, 'stop')

    END IF
    ALLOCATE(CONC_MAT(ROWCOUNT(51,'#'),N_VARS))
    ALLOCATE(TIMEVEC(ROWCOUNT(51,'#')))
    CLOSE(51)
    ! Deal with a situation where we have no input. We still need conc_mat and timevec.
  ELSE
    ALLOCATE(CONC_MAT(2,N_VARS))
    ALLOCATE(TIMEVEC(2))
    TIMEVEC = (/0d0, MODELTIME%SIM_TIME_H/)
  END IF
  CONC_MAT = 0d0

  ! READ ENV INPUT
  if (ENV_file /= '') THEN
    OPEN(unit=51, File=TRIM(ENV_file), ACTION='READ', STATUS='OLD', iostat=ioi)
    CALL handle_file_io(ioi, ENV_file, 'Terminating on CONCMAT allocation')

    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)

    ALLOCATE(INPUT_ENV(rowcol_count%rows,rowcol_count%cols))
    INPUT_ENV = 0
    call FILL_INPUT_BUFF(51,rowcol_count,INPUT_ENV,ENV_file)
    timevec = INPUT_ENV(:,1)
    CLOSE(51)
  END IF

  ! READ MCM INPUT
  if (MCM_file /= '') THEN
    OPEN(unit=51, File=TRIM(MCM_file), ACTION='READ', STATUS='OLD')
    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)
    allocate(INPUT_MCM(rowcol_count%rows,rowcol_count%cols))
    INPUT_MCM = 0
    call FILL_INPUT_BUFF(51,rowcol_count,INPUT_MCM,MCM_file)
    timevec = INPUT_MCM(:,1)
    CLOSE(51)
  END IF

  CALL PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)

  ! check IF dmps data is used or not. If no then do nothing. Use dmps only
  !if aerosol_flag is true
  IF (USE_DMPS .and. Aerosol_flag) then

    write(*,FMT_SUB) 'Reading DMPS files '// TRIM(DMPS_file)
    OPEN(unit=51, File=TRIM(DMPS_file), STATUS='OLD', iostat=ioi)
    CALL handle_file_io(ioi, DMPS_file, 'Terminating on DMPS pardata')

    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)

    ! Allocate par_data, for this we find out if the particle file has the total particle number as second column in
    ! which case it gets omitted
    read(51,*) fl_buff
    IF (fl_buff(2) > 0d0) THEN
      allocate(par_data(rowcol_count%rows,rowcol_count%cols))
    ELSE
      allocate(par_data(rowcol_count%rows,rowcol_count%cols-1))
    END IF
    REWIND(51)

    DO I=1,rowcol_count%rows
      IF (fl_buff(2) > 0d0) then
        read(51,*) par_data(I,:)
      END IF
      IF (fl_buff(2) < 1d-200) THEN
        read(51,*) par_data(I,1), fl_buff(1), par_data(I,2:)
      END IF
    END DO
    CLOSE(51)

    ! BG_PAR is here in case you want to use it Carlton and Lukas, in the end we remove either BG_PAR or par_data mmkay.
    allocate(BG_PAR%binseries(size(par_data(:,1))-1, size(par_data(1,:))-1))
    allocate(BG_PAR%options(1))
    allocate(BG_PAR%time(size(par_data(:,1))-1))
    allocate(BG_PAR%sections(size(par_data(1,:))-1))
    BG_PAR%name = 'BACKGROUND_CONC'
    BG_PAR%time = par_data(2:, 1)
    BG_PAR%sections = par_data(1, 2:)
    BG_PAR%binseries = par_data(2:, 2:)

  END IF



  IF (extra_particles /= '') THEN
    ! First we open the extra particle files to count the dimensions needed for the matrix
    OPEN(unit=51, File=TRIM(extra_particles) , STATUS='OLD', iostat=ioi)
    Z = ROWCOUNT(51)
    write(buf, '(a, i0, a)') 'reading XTRAS: ', Z,' lines'
    PRINT FMT_SUB, TRIM(buf)
    ! Now we can allocate it
    allocate(XTRAS(Z))

    DO I=1,Z
      ! PRINT*,'Z', Z
      allocate(XTRAS(I)%options(n_xpar_options))

      read(51,'(a)') buf

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

      OPEN(unit=51+I, File=TRIM(buf(1:path_l(1))) , STATUS='OLD', iostat=ioi)
      CALL handle_file_io(ioi, buf(1:path_l(1)), 'Terminating when trying to open extra particles, check the paths in extra_particles')

      yp = ROWCOUNT(51+I)
      xp = COLCOUNT(51+I)

      allocate(XTRAS(I)%time(yp-1))

      read(51+I,*) fl_buff
      IF (fl_buff(2) > 0d0) THEN
        allocate(XTRAS(I)%binseries(yp-1,xp-1))
        allocate(XTRAS(I)%sections(xp-1))
      ELSE
        allocate(XTRAS(I)%binseries(yp-1,xp-2))
        allocate(XTRAS(I)%sections(xp-2))
      END IF
      REWIND(51+I)

      DO j=1,yp
        if (j==1) THEN
          read(51+I,*) fl_buff(1:xp-size(XTRAS(I)%sections)), XTRAS(I)%sections
        ELSE
          read(51+I,*) fl_buff(1:xp-size(XTRAS(I)%sections)), XTRAS(I)%binseries(J-1,:)
          XTRAS(I)%time(j-1) = fl_buff(1)
        END IF
      END DO
      CLOSE(51+I)

    END DO
    CLOSE(51)
    do i=1,Z
      print FMT_MSG, 'Extra particle input for '//XTRAS(i)%name
    end do
  END IF





  print FMT_LEND,

  IF (VAP_logical) then
   write(*,FMT_MSG) 'Reading Vapour name file '// TRIM(Vap_names)
   write(*,FMT_MSG) 'Reading Vapour prop file '// TRIM(Vap_props)
   OPEN(unit=52, File= TRIM(Vap_names) , STATUS='OLD', iostat=ioi)
   OPEN(unit=53, File= TRIM(Vap_props) , STATUS='OLD', iostat=ioi2)

   IF (ioi /= 0) THEN
     print FMT_FAT0, 'Vap_names file was defined but not readable, exiting. Check NML_vap in INIT file'
     STOP
   ELSEIF (ioi2 /= 0) THEN
       print FMT_FAT0, 'Vap_props was defined but not readable, exiting. Check NML_vap in INIT file'
       STOP
   END IF

   rowcol_count%rows = ROWCOUNT(52)
   rowcol_count%cols = COLCOUNT(53)

   allocate(vapours%Vapour_names(rowcol_count%rows + 1 ))
   allocate(vapours%molar_mass(rowcol_count%rows + 1))
   allocate(vapours%parameter_A(rowcol_count%rows + 1))
   allocate(vapours%parameter_B(rowcol_count%rows + 1))
   allocate(vapours%molec_mass(rowcol_count%rows + 1))
   allocate(vapours%molec_volume(rowcol_count%rows + 1))
   allocate(vapours%density(rowcol_count%rows + 1))
   allocate(vapours%surf_tension(rowcol_count%rows + 1))
   allocate(vapours%c_sat(rowcol_count%rows + 1))
   allocate(vapours%vap_conc(rowcol_count%rows + 1))
   allocate(vapours%cond_type(rowcol_count%rows + 1))
   allocate(vapours%molec_dia(rowcol_count%rows + 1))
   allocate(vapours%mfractions(rowcol_count%rows + 1))
   allocate(vapours%alpha(rowcol_count%rows + 1))

   vapours%vapour_number = rowcol_count%rows
   vapours%vbs_bins      = rowcol_count%rows + 1

   write(buf,'(i0)') vapours%vapour_number
   print FMT_SUB, 'Vapours available = '//TRIM(buf)
   write(buf,'(i0)') vapours%vbs_bins
   print FMT_SUB, 'Vapour bins + H2SO4 = '//TRIM(buf)
   !!! reading the vap names and vap vapour_properties

   vapours%Mfractions(vapours%vapour_number) = 1 !

 do ii = 1, vapours%vapour_number+1
   if (ii <= vapours%vapour_number) then !!! all compounds
     read(52,*,iostat=iosp)   species_name
     read(53,*,iostat=ioprop) molar_mass, parameter_A, parameter_B
     vapours%molar_mass(ii)    = molar_mass *1D-3 ! kg/mol
     vapours%parameter_A(ii)   = parameter_A
     vapours%parameter_B(ii)   = parameter_B
     vapours%vapour_names(ii)  = TRIM(species_name)
     vapours%molec_mass(ii)    = calculate_molecular_mass(vapours%molar_mass(ii))  !kg/#
     vapours%density(ii)       = 1400.0  !!! kg/m3
     vapours%molec_volume(ii)  = calculate_molecular_volume(vapours%density(ii),vapours%molec_mass(ii))
     vapours%surf_tension(ii)  = 0.05
     vapours%cond_type(ii)     = 1  ! not an acid (H2SO4 or HCL)
     vapours%molec_dia(ii)     = (6D0 * vapours%molec_volume(ii) / pi )**(1D0/3D0)  ! molecular diameter [m]
     ! vapours%mfractions(ii)    = 0.0
     vapours%alpha(ii)         = 1.0
     vapours%c_sat(ii)         = calculate_saturation_vp(vapours%parameter_A(ii),vapours%parameter_B(ii), ambient%temperature)
   else !! h2so4
     vapours%molar_mass(ii)    = 98.0785 *1d-3
     vapours%parameter_A(ii)   = 3.869717803774
     vapours%parameter_B(ii)   = 313.607405085
     vapours%vapour_names(ii)  = 'H2S04'
     vapours%molec_mass(ii)    = calculate_molecular_mass(vapours%molar_mass(ii))
     vapours%density(ii)       = 1819.3946 !!! kg/m3
     vapours%molec_volume(ii)  = calculate_molecular_volume(vapours%density(ii),vapours%molec_mass(ii))
     vapours%surf_tension(ii)  = 0.07
     vapours%cond_type(ii)     = 2  ! Acid
     vapours%molec_dia(ii)     = (6D0 * vapours%molec_volume(ii) / pi )**(1D0/3D0)  ! molecular diameter [m]
     vapours%mfractions(ii)    = 0.0
     vapours%alpha(ii)         = 1.0
     vapours%c_sat(ii)         = 0.0 !calculate_saturation_vp(vapours%parameter_A(ii),vapours%parameter_B(ii), ambient%temperature)
   end if
 end do


   close(52)
   close(53)


  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Calculating the diffusivity of organics based on fullers method !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (Use_atoms) THEN

    write(*,FMT_MSG) 'Reading O_C file '// TRIM(Vap_atoms)
    ! OPEN(unit=61, File=TRIM(ADJUSTL(CASE_DIR)) // '/'//TRIM(Vap_atoms) , STATUS='OLD', iostat=ioi3)
    OPEN(unit=61, File=TRIM(Vap_atoms) , STATUS='OLD', iostat=ioi3)

    rowcol_count%rows = ROWCOUNT(61)
    rowcol_count%cols = COLCOUNT(61)

    allocate(Natoms%N_Carbon(rowcol_count%rows))
    allocate(Natoms%N_Hydrogen(rowcol_count%rows))
    allocate(Natoms%N_Oxygen(rowcol_count%rows))
    allocate(Natoms%N_Nitrogen(rowcol_count%rows))
    allocate(Natoms%comp_prop(rowcol_count%cols, rowcol_count%rows))
    allocate(Vol_org(rowcol_count%rows))
    allocate(Diff_org(rowcol_count%rows))

    READ(61,*) Natoms%comp_prop

    Natoms%N_Carbon   = Natoms%comp_prop(2,:)
    Natoms%N_Oxygen   = Natoms%comp_prop(3,:)
    Natoms%N_Nitrogen = Natoms%comp_prop(4,:)
    Natoms%N_Hydrogen = Natoms%comp_prop(8,:)

   !Diffusivity
   ! dimensionless diffusion volumes
   ! Atomic diffusion values from (Tang et al., ACP 2015) , Reid et al., 1987

   Vol_org = Natoms%N_Carbon*15.9D0 + Natoms%N_Oxygen*6.11D0+&
        Natoms%N_Hydrogen*2.31D0 + Natoms%N_Nitrogen*4.54D0

  ! Fuller's method /! Diffusivity organic compounds (m^2 s^-1)
  ! http://umich.edu/~elements/course/lectures/eleven/exam4.htm

   Diff_org = 1D-7*ambient%temperature**1.75D0*SQRT(1/(Mair*1D3) + 1/Natoms%comp_prop(5,:)) / &
   (ambient%pressure/1.01325D5*(Vol_org**(1D0/3D0) + 20.1**(1D0/3D0))**2D0)


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

  OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=IOS(1))
  ! Handle file not found error
  IF (IOS(1) /= 0) THEN
    write(*,FMT_FAT0) 'No INITFILE called '//TRIM(ADJUSTL(Fname_init))//', exiting. Good bye.'
    write(*,FMT_LEND)
    STOP
  END IF

  ! if INITFILE was found, we read it. In case there is a problem in namelist filling, give en error.
  write(*,FMT_HDR) 'READING USER DEFINED INTIAL VALUES FROM: '//TRIM(ADJUSTL(Fname_init))

  do k=1, ROWCOUNT(50); READ(50,NML = NML_TIME, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_Flag, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_Path, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_MISC, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_VAP, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_PARTICLE, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_ENV, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_MCM, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  do k=1, ROWCOUNT(50); READ(50,NML = NML_MODS, IOSTAT=IOS(i))
  IF (IOS(i) == 0) EXIT;end do; REWIND(50); i=i+1

  CLOSE(50)

  IF (SUM(ABS(IOS)) /= 0) then
    write(*,FMT_MSG) 'There was a problem with INITFILE. Check the file and refer to manual. Exiting now.'
    write(*,*) IOS
    write(*,FMT_LEND)
    STOP
  end if

end subroutine READ_INIT_FILE

subroutine PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME
  implicit none
  INTEGER :: days(12), non_leap_year(12) = ([31,28,31,30,31,30,31,31,30,31,30,31])
  INTEGER :: leap_year(12) = ([31,29,31,30,31,30,31,31,30,31,30,31])
  INTEGER :: y,m,d, ioi
  CHARACTER(100) :: buf
  MODELTIME%SIM_TIME_H = runtime
  MODELTIME%SIM_TIME_S = runtime*3600d0
  ! figure out the correct save interval
  IF (FSAVE_DIVISION > 0) THEN
    MODELTIME%FSAVE_INTERVAL = INT(MODELTIME%SIM_TIME_S/MODELTIME%dt) / FSAVE_DIVISION * MODELTIME%dt
  ELSEIF (FSAVE_INTERVAL > 0) THEN
    MODELTIME%FSAVE_INTERVAL = FSAVE_INTERVAL
  END IF
  MODELTIME%PRINT_INTERVAL = PRINT_INTERVAL
  ! IF julian day was provided, use it
  IF (JD > 0) MODELTIME%JD = JD
  if (TRIM(DATE) /= '1800-01-01') THEN

    read(date(1:4),*,iostat=ioi) y
    if (ioi == 0) read(date(6:7),*,iostat=ioi) m
    if (ioi == 0) read(date(9:) ,*,iostat=ioi) d
    if (ioi /= 0) THEN
      print FMT_WARN0, 'Proper date not provided in INITFILE'
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
      MODELTIME%PRINTACDC = PRINT_ACDC
      MODELTIME%JD = sum(days(:m-1)) + d
      write(buf, '(i0)') MODELTIME%JD
      print FMT_MSG, 'Date: '//TRIM(date)//' -> Julian Day: '//TRIM(buf)
    END IF
  END if
end subroutine PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME

subroutine NAME_MODS_SORT_NAMED_INDICES
  implicit none
  INTEGER :: i
  OPEN(2151, file='src/NAMES.dat', ACTION='READ', status='OLD')
  DO i = 1,N_VARS
    READ(2151, *) MODS(I)%NAME
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
  close(2151)

end subroutine NAME_MODS_SORT_NAMED_INDICES

subroutine REPORT_INPUT_COLUMNS_TO_USER
  implicit none
  integer::i
  character(4) :: buffer
  DO i=1,N_VARS
    IF (I==1) print FMT_MSG, 'ENV values from '//TRIM(ENV_file)//':'
    IF ((I==LENV) .and. (maxval(MODS(LENV:)%col)>-1)) print FMT_MSG, 'MCM values from '//TRIM(ENV_file)//':'

    IF (MODS(I)%col > -1) THEN
      write(buffer, '(i0)') MODS(I)%col
      IF ((TRIM(ENV_file) == '') .and. (I<LENV)) THEN
        print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(buffer)//' but no ENV_FILE. Using parametrisation?'
      ELSE IF ((TRIM(MCM_file) == '') .and. (I>LENV)) THEN
        print FMT_SUB, TRIM(MODS(I)%NAME)//' should be read from column: '//TRIM(buffer)//' but no MCM_FILE. Using parametrisation?'
      ELSE
        print FMT_SUB, TRIM(MODS(I)%NAME)//' will be read from column: '//TRIM(buffer)
      END IF
    END IF
  END DO

end subroutine REPORT_INPUT_COLUMNS_TO_USER

SUBROUTINE PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)
  implicit none
  REAL(DP), intent(in)    :: INPUT_ENV(:,:)
  REAL(DP), intent(in)    :: INPUT_MCM(:,:)
  REAL(DP), intent(inout) :: CONC_MAT(:,:)
  integer :: i

  DO i=1,N_VARS
    IF ((I < lenv) .and. (ENV_file /= '') .and. (MODS(I)%col > -1)) THEN
      CONC_MAT(:,I) = input_ENV(:,MODS(I)%col)
    END IF

    IF ((I>lenv) .and. (MCM_file /= '') .and. (MODS(I)%col > -1)) THEN
      CONC_MAT(:,I) = input_MCM(:,MODS(I)%col)
    END IF

  END DO

END SUBROUTINE PUT_INPUT_IN_THEIR_PLACES

subroutine FILL_INPUT_BUFF(unit,rowcol_count, INPUT_BF,Input_file)
  implicit none
  type(nrowcol), intent(in) :: rowcol_count
  real(dp), intent(inout) :: INPUT_BF(:,:)
  character(*) :: Input_file
  integer :: i,j,k, ioi, unit

  ! Reading data into the variable
  i = 1
  print FMT_MSG, 'Filling input matrices...'
  DO k = 1, ROWCOUNT(unit)
    READ(unit,*, iostat=ioi) (INPUT_BF(i,j),j=1,rowcol_count%cols)
    IF ((ioi /= 0) .and. (i==1)) THEN
      print FMT_SUB, 'First row omitted from file "'// TRIM(Input_file) //'".'
    ELSE IF ((ioi /= 0) .and. (i>1)) THEN
      print FMT_WARN1, 'Bad value in file '// TRIM(Input_file) //'". Maybe an illegal value in line number', real(k)
    ELSE
      i=i+1
    END IF
  END DO
  print FMT_MSG, 'Done filling input matrices...'
end subroutine FILL_INPUT_BUFF


! ================================================================================================
! Print out which modifiers differ from default values. If SHIFTER or MULTIPLYER differ from their
! null value, this subroutine will pick them and print them for the user before the main loop starts
! ================================================================================================
SUBROUTINE CHECK_MODIFIERS()
  IMPLICIT NONE
  type(input_mod) :: test
  integer         :: i,j=0
  character(4)    :: cprf
  character(20)    :: buf

  print FMT_HDR, 'Check input validity'

  ! Save input temperature and pressure for archive purpose
  ORIGINAL_TEMP(1)  = MODS(inm_TempK)
  ORIGINAL_press(1) = MODS(inm_pres)

  CALL CONVERT_TEMPS_TO_KELVINS
  CALL CONVERT_PRESSURE_AND_VALIDATE_UNITS

  do i=1,size(MODS)
      IF (MODS(i)%MODE > 0) THEN
          print FMT_MSG, '- Replacing input for '//TRIM(MODS(i)%name)//' with parametrized function.'
          j=1
      ELSE
          IF (ABS(MODS(i)%multi - test%multi) > 1d-9) THEN
              write(buf, '(es9.3)') MODS(i)%multi
              print FMT_MSG, '- Multiplying '//TRIM(MODS(i)%name)//' with: '//TRIM(buf)
              j=1
          END IF
          if (ABS(MODS(i)%shift - test%shift) > 1d-9) THEN
              if (TRIM(MODS(i)%UNIT) == '#') THEN
                  cprf = '/cm3'
              else
                  cprf = ''
              end if
              write(buf, '(es9.3)') MODS(i)%shift
              print FMT_MSG, '- Adding a constant to '//TRIM(MODS(i)%name)//', [in '//TRIM(MODS(i)%UNIT)//TRIM(cprf)//']: '//TRIM(buf)
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
    !use constants, ONLY: UCASE
    IMPLICIT NONE

    if ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == '#') THEN
        print FMT_WARN0, "No unit for temperature. Use either 'K' or 'C'. Now assuming Kelvins."
        TempUnit = 'K'
    elseif ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'K') THEN
        TempUnit = 'K'
    elseif ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'C') THEN
        TempUnit = 'C'
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
    ! Check if a double conversion is attempted
    IF ((TempUnit == 'C') .and. (  ABS(MODS(inm_TempK)%shift - 273.15)<1d0  )) THEN
        print FMT_WARN1, 'Temperature will be converted to Kelvins, but an additional constant is added: ',MODS(inm_TempK)%shift
    END IF
    TempUnit = 'K'
    MODS(inm_TempK)%UNIT = TempUnit
  END SUBROUTINE CONVERT_TEMPS_TO_KELVINS


! ================================================================================================
! This subroutine converts pressure from all possible input units to Pa
! ================================================================================================
SUBROUTINE CONVERT_PRESSURE_AND_VALIDATE_UNITS
  !use constants, ONLY: UCASE
  IMPLICIT NONE
  INTEGER      :: i
  character(5) :: buf

  buf = UCASE(TRIM(MODS(inm_pres)%UNIT))
  if (TRIM(buf) == 'HPA' .or. TRIM(buf) == 'MBAR') THEN
      ! CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 100d0
      ! MODS(inm_pres)%shift = MODS(inm_pres)%shift *100d0
      print FMT_MSG, '- Converting pressure from hPa (mbar) to Pascals.'
    elseif (TRIM(buf) == 'KPA') THEN
      ! CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 1000d0
      ! MODS(inm_pres)%shift = MODS(inm_pres)%shift *1000d0
      print FMT_MSG, '- Converting pressure from kPa to Pascals.'
    elseif (TRIM(buf) == 'ATM') THEN
      ! CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 1.01325d5
      ! MODS(inm_pres)%shift = MODS(inm_pres)%shift * 1.01325d5
      print FMT_MSG, '- Converting pressure from atm to Pascals.'
    elseif (TRIM(buf) == 'PA') THEN
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




end module INPUT
