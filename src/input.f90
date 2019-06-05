Module INPUT
! Module to read init file

USE second_precision
use constants
USE auxillaries
USE Aerosol_auxillaries

Implicit none

INTEGER :: N_VARS ! This will store the number of variables in NAMES.DAT
INTEGER :: LENV   ! This will store the number of named indices in this code

! INDICES TO PROPERLY COMBINE INPUT TO CORRECT VALUES IN THE MODEL
! ALL VARIABLES THAT ARE TIME DEPENDENT MUST BE HERE
! IF YOU ADD VARIABLES HERE, YOU NEED TO UPDATE:
! - this list
!------------------------------------------------------------------
INTEGER :: inm_tempK
INTEGER :: inm_pres
INTEGER :: inm_RH
INTEGER :: inm_CS
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
character(len=256):: CASE_DIR   = ''
character(len=256):: CASE_NAME  = ''
character(len=30) :: RUN_NAME   = ''
NAMELIST /NML_Path/ Work_dir, Case_Dir, Case_name, RUN_NAME

! MODULES IN USE OPTIONS
Logical :: Aerosol_flag   = .false.
Logical :: Chemistry_flag = .false.
Logical :: Particle_flag  = .false.
Logical :: ACDC_solve_ss  = .false.
Logical :: NUCLEATION     = .true.
Logical :: ACDC           = .true.
Logical :: Extra_data     = .false.
Logical :: Current_case   = .false.
NAMELIST /NML_Flag/ Aerosol_flag, chemistry_flag, particle_flag,ACDC_solve_ss,NUCLEATION,ACDC,Extra_data, Current_case

! TIME OPTIONS
real(dp)  :: runtime = 1d0
real(dp)  :: FSAVE_INTERVAL = 300d0
real(dp)  :: PRINT_INTERVAL = 15*60d0
INTEGER   :: FSAVE_DIVISION = 0
NAMELIST /NML_TIME/ runtime, FSAVE_INTERVAL, PRINT_INTERVAL, FSAVE_DIVISION

! MODIFIER OPTIONS
! MODS and UNITS are declared in CONSTANTS.f90, in order to be more widely available
NAMELIST /NML_MODS/ MODS
NAMELIST /NML_UNITS/ UNITS

! DMPS INPUT
character(len=256)  :: DMPS_dir
character(len=256)  :: DMPS_file
REAL(dp)            :: read_in_time = 0d0 ![seconds] !for use_dmps_special, read dmps data above this cut_off_diameter(m)
REAL(dp)            :: dmps_upper_band_limit = 18.*1d-9 !for use_dmps_special, read dmps data above this cut_off_diameter(m)
REAL(dp)            :: dmps_lower_band_limit = 6.*1d-10 !for use_dmps_special, read dmps data below this take_in_diameter(m)
logical             :: use_dmps = .false.
logical             :: use_dmps_special = .false.
NAMELIST /NML_DMPS/ DMPS_dir, DMPS_file,read_in_time,dmps_upper_band_limit, dmps_lower_band_limit,&
use_dmps,use_dmps_special

! ENVIRONMENTAL INPUT
character(len=256)  :: ENV_path = ''
character(len=256)  :: ENV_file = ''
character(len=1)  :: TempUnit = 'K' ! or C
NAMELIST /NML_ENV/ ENV_path, ENV_file, TempUnit

! MCM INPUT
character(len=256)  :: MCM_path = ''
character(len=256)  :: MCM_file = ''
NAMELIST /NML_MCM / MCM_path, MCM_file

! MISC OPTIONS
real(dp)  :: lat
real(dp)  :: lon
INTEGER   :: JD = -1
CHARACTER(1000)  :: Description
NAMELIST /NML_MISC/ JD, lat, lon, Description

Logical  :: VAP_logical = .False.
character(len=256)  :: Vap_names
character(len=256)  :: Vap_props
NAMELIST /NML_VAP/ VAP_logical, Vap_names, Vap_props

type(vapour_ambient) :: vapours

contains

subroutine READ_INPUT_DATA()
  IMPLICIT NONE
  character(len=256)  :: data_dir
  type(nrowcol)       :: rowcol_count
  integer             :: ioi,ioi2, ii, iosp, ioprop
  !!! for vapour FILES
  character(len=256)  :: species_name
  real(dp)            :: molar_mass, parameter_A, parameter_B


  ! Welcoming message
  print'(a,t23,a)', achar(10),  '--~:| Gas and Aerosol Box Model - GAB v.0.1 |:~--'//achar(10)

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
  ALLOCATE(UNITS(N_VARS))
  UNITS = '#'
  ALLOCATE(INDRELAY_CH(N_VARS))
  INDRELAY_CH = 0
  CALL NAME_MODS_SORT_NAMED_INDICES
  CALL READ_INIT_FILE
  CALL FILL_INDRELAY_CH_WITH_INDICES
  CALL PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME

  ! ALLOCATE CONC_MAT
  IF ((ENV_file /= '') .or. (MCM_file /= '')) THEN
    IF (ENV_file /= '') THEN
      OPEN(unit=51, File=TRIM(ENV_path) //'/'//TRIM(ENV_file), ACTION='READ', STATUS='OLD')
    ELSE
      OPEN(unit=51, File=TRIM(MCM_path) //'/'//TRIM(MCM_file), ACTION='READ', STATUS='OLD')
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
    OPEN(unit=51, File=TRIM(ENV_path) //'/'//TRIM(ENV_file), ACTION='READ', STATUS='OLD')
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
    OPEN(unit=51, File=TRIM(MCM_path) //'/'//TRIM(MCM_file), ACTION='READ', STATUS='OLD')
    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)
    allocate(INPUT_MCM(rowcol_count%rows,rowcol_count%cols))
    INPUT_MCM = 0
    call FILL_INPUT_BUFF(51,rowcol_count,INPUT_MCM,MCM_file)
    timevec = INPUT_MCM(:,1)
    CLOSE(51)
  END IF

  CALL PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)

  ! check IF dmps data is used or not. If no then do nothing
  IF (USE_DMPS) then
   write(*,FMT_SUB),'Reading DMPS file '// TRIM(DMPS_file)
   OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' //TRIM(DMPS_dir)// '/'//TRIM(DMPS_file) , STATUS='OLD', iostat=ioi)
   IF (ioi /= 0) THEN
     print FMT_FAT0, 'DMPS file was defined but not readable, exiting. Check NML_DMPS in INIT file'
     STOP
   END IF
   rowcol_count%rows = ROWCOUNT(51,'#')
   rowcol_count%cols = COLCOUNT(51)
   allocate(par_data(rowcol_count%rows,rowcol_count%cols))
  end if

  CLOSE(51)
  print FMT_LEND,

  IF (VAP_logical) then
   write(*,FMT_SUB),'Reading Vapour name file '// TRIM(Vap_names)
   write(*,FMT_SUB),'Reading Vapour prop file '// TRIM(Vap_props)
   OPEN(unit=52, File=TRIM(ADJUSTL(CASE_DIR)) // '/'//TRIM(Vap_names) , STATUS='OLD', iostat=ioi)
   OPEN(unit=53, File=TRIM(ADJUSTL(CASE_DIR)) // '/'//TRIM(Vap_props) , STATUS='OLD', iostat=ioi2)

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

   vapours%vapour_number = rowcol_count%rows
   vapours%vbs_bins      = rowcol_count%rows + 1
   print*, 'Vapours available= ',  vapours%vapour_number
   print*, 'Vapour bins + H2SO4', vapours%vbs_bins
   !!! reading the vap names and vap vapour_properties

 do ii = 1, vapours%vapour_number+1
   if (ii <= vapours%vapour_number) then !!! all compounds
     read(52,*,iostat=iosp)   species_name
     read(53,*,iostat=ioprop) molar_mass, parameter_A, parameter_B
     vapours%molar_mass(ii)    = molar_mass
     vapours%parameter_A(ii)   = parameter_A
     vapours%parameter_B(ii)   = parameter_B
     vapours%vapour_names(ii)  = TRIM(species_name)
     vapours%molec_mass(ii)    = calculate_molecular_mass(molar_mass)
     vapours%molec_volume(ii)  = calculate_molecular_volume(vapours%molec_mass(ii), vapours%density)
   else !! h2so4
     vapours%molar_mass(ii)    = 98.0785
     vapours%parameter_A(ii)   = 3.869717803774
     vapours%parameter_B(ii)   = 313.607405085
     vapours%vapour_names(ii)  = 'H2S04'
     vapours%molec_mass(ii)    = calculate_molecular_mass(vapours%molar_mass(ii))
     vapours%molec_volume(ii)  = calculate_molecular_volume(vapours%molec_mass(ii), vapours%density)
   end if
 end do

   close(52)
   close(53)

  end if


end subroutine READ_INPUT_DATA


subroutine READ_INIT_FILE
!-------------------------------------------------------------------------------
! Reads the init file and fills user-provided variables
!-------------------------------------------------------------------------------
  implicit none
  integer :: IOS(20)
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

  ! Handle file not found error
  IOS = 0
  OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=IOS(1))
  IF (IOS(1) /= 0) THEN
    write(*,FMT_FAT0) 'No INITFILE called '//TRIM(ADJUSTL(Fname_init))//', exiting. Good bye.'
    write(*,FMT_LEND)
    STOP
  END IF
  ! INITFILE was found, reading it. In case there is a problem in namelist filling, give en error.
  write(*,FMT_HDR), 'READING USER DEFINED INTIAL VALUES FROM: '//TRIM(ADJUSTL(Fname_init))
  READ(50,NML=NML_Path,  IOSTAT= IOS(2)) ! directories and test cases
    IF (IOS(2) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Path, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_Flag,  IOSTAT= IOS(3)) ! flags
    IF (IOS(3) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Flag, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_TIME,  IOSTAT= IOS(4)) ! time related stuff
    IF (IOS(4) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_TIME, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_DMPS,  IOSTAT= IOS(5)) ! dmps_file information
    IF (IOS(5) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_DMPS, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_ENV,  IOSTAT= IOS(6)) ! environmental information
    IF (IOS(6) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Temp, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_MCM,   IOSTAT= IOS(7)) ! MCM_file information
    IF (IOS(7) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MCM, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_MODS,  IOSTAT= IOS(8)) ! modification parameters
    IF (IOS(8) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MODS, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_UNITS,  IOSTAT= IOS(9)) ! units
    IF (IOS(9) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_UNITS, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_MISC,  IOSTAT= IOS(10)) ! misc input
    IF (IOS(10) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MISC, maybe some undefinded INITFILE input?'
  READ(50,NML=NML_VAP,  IOSTAT= IOS(11)) ! vapour input
      IF (IOS(11) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_VAP, maybe some undefinded INITFILE input?'
  CLOSE(50)

  IF (SUM(IOS) /= 0) then
    write(*,FMT_MSG) 'There was a problem with INITFILE. Check the file and refer to manual. Exiting now.'
    write(*,FMT_LEND)
    STOP
  end if

end subroutine READ_INIT_FILE

subroutine PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME
  implicit none
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
    IF (TRIM(MODS(I)%NAME) == '#'            ) LENV = i

  END DO
  close(2151)

end subroutine NAME_MODS_SORT_NAMED_INDICES

subroutine FILL_INDRELAY_CH_WITH_INDICES
  implicit none
  integer::i
  character(4) :: buffer
  DO i=1,N_VARS
    IF (I==1) print FMT_MSG, 'Values from '//TRIM(ENV_file)//':'
    IF (I==LENV) print FMT_MSG, 'Values from '//TRIM(ENV_file)//':'
    IF (MODS(I)%col > -1) THEN
      write(buffer, '(i0)') MODS(I)%col
      print FMT_SUB, 'Values for '//TRIM(MODS(I)%NAME)//' will be read from column: '//TRIM(buffer)
    END IF
  END DO

end subroutine FILL_INDRELAY_CH_WITH_INDICES

SUBROUTINE PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)
  implicit none
  REAL(DP), intent(in)    :: INPUT_ENV(:,:)
  REAL(DP), intent(in)    :: INPUT_MCM(:,:)
  REAL(DP), intent(inout) :: CONC_MAT(:,:)
  integer :: i

  DO i=1,N_VARS
    IF (MODS(I)%NAME == '#') lenv = i
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




end module INPUT
