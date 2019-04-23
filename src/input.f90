Module INPUT
! Module to read init file

USE second_precision
use constants
USE auxillaries

Implicit none

! INDICES TO PROPERLY COMBINE INPUT TO CORRECT VALUES IN THE MODEL
! ALL VARIABLES THAT ARE TIME DEPENDENT MUST BE HERE
! IF YOU ADD VARIABLES HERE, YOU NEED TO UPDATE:
! - this list
! - parameter inm_LAST
! - NAMELIST NML_ICOLS
! - SUBROUTINE PUT_INPUT_IN_THEIR_PLACES
! - subroutine NAME_MODS
!------------------------------------------------------------------

! The following will be provided in ENV_file
INTEGER :: inf_H2SO4 = -1 ; INTEGER, parameter :: inm_H2SO4 = 1
INTEGER :: inf_NH3   = -1 ; INTEGER, parameter :: inm_NH3   = 2
INTEGER :: inf_DMA   = -1 ; INTEGER, parameter :: inm_DMA   = 3
INTEGER :: inf_CS    = -1 ; INTEGER, parameter :: inm_CS    = 4
INTEGER :: inf_swr   = -1 ; INTEGER, parameter :: inm_swr   = 5
INTEGER :: inf_RH    = -1 ; INTEGER, parameter :: inm_RH    = 6
INTEGER :: inf_pres  = -1 ; INTEGER, parameter :: inm_pres  = 7
INTEGER :: inf_temp  = -1 ; INTEGER, parameter :: inm_temp  = 8
INTEGER :: inf_SO2   = -1 ; INTEGER, parameter :: inm_SO2   = 9
INTEGER :: inf_NO    = -1 ; INTEGER, parameter :: inm_NO    = 10
INTEGER :: inf_NO2   = -1 ; INTEGER, parameter :: inm_NO2   = 11
INTEGER :: inf_CO    = -1 ; INTEGER, parameter :: inm_CO    = 12
INTEGER :: inf_H2    = -1 ; INTEGER, parameter :: inm_H2    = 13
INTEGER :: inf_O3    = -1 ; INTEGER, parameter :: inm_O3    = 14
INTEGER :: inf_IPR   = -1 ; INTEGER, parameter :: inm_IPR   = 15
! The following will be provided in VOC_file
INTEGER :: inf_CH3O  = -1 ; INTEGER, parameter :: inm_CH3O  = 16
INTEGER :: inf_CH3C  = -1 ; INTEGER, parameter :: inm_CH3C  = 17
INTEGER :: inf_C2H5  = -1 ; INTEGER, parameter :: inm_C2H5  = 18
INTEGER :: inf_C5H8  = -1 ; INTEGER, parameter :: inm_C5H8  = 19
INTEGER :: inf_MVK   = -1 ; INTEGER, parameter :: inm_MVK   = 20
INTEGER :: inf_MEK   = -1 ; INTEGER, parameter :: inm_MEK   = 21
INTEGER :: inf_BENZ  = -1 ; INTEGER, parameter :: inm_BENZ  = 22
INTEGER :: inf_APIN  = -1 ; INTEGER, parameter :: inm_APIN  = 23
INTEGER :: inf_BPIN  = -1 ; INTEGER, parameter :: inm_BPIN  = 24
INTEGER :: inf_LIMO  = -1 ; INTEGER, parameter :: inm_LIMO  = 25
INTEGER :: inf_Care  = -1 ; INTEGER, parameter :: inm_Care  = 26
INTEGER :: inf_TOLU  = -1 ; INTEGER, parameter :: inm_TOLU  = 27
! Some reserves for VOC, create more if these run out
INTEGER :: inf_RES20 = -1 ; INTEGER, parameter :: inm_RES20 = 28
INTEGER :: inf_RES19 = -1 ; INTEGER, parameter :: inm_RES19 = 29
INTEGER :: inf_RES18 = -1 ; INTEGER, parameter :: inm_RES18 = 30
INTEGER :: inf_RES17 = -1 ; INTEGER, parameter :: inm_RES17 = 31
INTEGER :: inf_RES16 = -1 ; INTEGER, parameter :: inm_RES16 = 32
INTEGER :: inf_RES15 = -1 ; INTEGER, parameter :: inm_RES15 = 33
INTEGER :: inf_RES14 = -1 ; INTEGER, parameter :: inm_RES14 = 34
INTEGER :: inf_RES13 = -1 ; INTEGER, parameter :: inm_RES13 = 35
INTEGER :: inf_RES12 = -1 ; INTEGER, parameter :: inm_RES12 = 36
INTEGER :: inf_RES11 = -1 ; INTEGER, parameter :: inm_RES11 = 37
INTEGER :: inf_RES10 = -1 ; INTEGER, parameter :: inm_RES10 = 38
INTEGER :: inf_RES9  = -1 ; INTEGER, parameter :: inm_RES9  = 39
INTEGER :: inf_RES8  = -1 ; INTEGER, parameter :: inm_RES8  = 40
INTEGER :: inf_RES7  = -1 ; INTEGER, parameter :: inm_RES7  = 41
INTEGER :: inf_RES6  = -1 ; INTEGER, parameter :: inm_RES6  = 42
INTEGER :: inf_RES5  = -1 ; INTEGER, parameter :: inm_RES5  = 43
INTEGER :: inf_RES4  = -1 ; INTEGER, parameter :: inm_RES4  = 44
INTEGER :: inf_RES3  = -1 ; INTEGER, parameter :: inm_RES3  = 45
INTEGER :: inf_RES2  = -1 ; INTEGER, parameter :: inm_RES2  = 46
INTEGER :: inf_RES1  = -1 ; INTEGER, parameter :: inm_RES1  = 47

INTEGER, parameter :: inm_LAST  = 47 ! <-whatever the last index on previous line is

! INDICES
NAMELIST /NML_ICOLS/ inf_H2SO4,inf_NH3,inf_DMA,inf_CS,inf_swr,&
inf_RH,inf_pres,inf_temp,inf_SO2,inf_NO,inf_NO2,inf_CO,inf_H2,inf_O3,&
inf_IPR,inf_CH3O,inf_CH3C,inf_C2H5,inf_C5H8,inf_MVK,inf_MEK,inf_BENZ,&
inf_APIN,inf_BPIN,inf_LIMO,inf_Care,inf_TOLU, inf_RES20,inf_RES19,inf_RES18,&
inf_RES17,inf_RES16,inf_RES15,inf_RES14,inf_RES13,inf_RES12,inf_RES11,inf_RES10,&
inf_RES9 ,inf_RES8 ,inf_RES7 ,inf_RES6 ,inf_RES5 ,inf_RES4 ,inf_RES3 ,inf_RES2 ,inf_RES1

INTEGER :: INDRELAY(inm_LAST,2)

REAL(dp), allocatable, private :: INPUT_ENV(:,:)  ! will be of shape ( len(timevec) : inm_last+1 )
REAL(dp), allocatable, private :: INPUT_VOC(:,:)  ! will be of shape ( len(timevec) : inm_last+1 )
REAL(dp), allocatable :: timevec(:)     ! Whatever the times were for ALL measurements
REAL(dp), allocatable :: CONC_MAT(:,:)  ! will be of shape ( len(timevec) : inm_last )
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
INTEGER   :: JD = -1
NAMELIST /NML_TIME/ runtime, FSAVE_INTERVAL, PRINT_INTERVAL, FSAVE_DIVISION, JD

! MODIFIER OPTIONS

type(input_mod)     :: MODS(inm_LAST) ! THIS VECTOR HOLDS ALL MODIFICATION PARAMETERS
NAMELIST /NML_MODS/ MODS

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
NAMELIST /NML_ENV/ ENV_path, ENV_file

! VOC INPUT
character(len=256)  :: VOC_path = ''
character(len=256)  :: VOC_file = ''
NAMELIST /NML_VOC / VOC_path, VOC_file

! MISC OPTIONS
real(dp)  :: lat
real(dp)  :: lon
CHARACTER(1000)  :: Description
NAMELIST /NML_MISC/ lat, lon, Description

contains

subroutine read_input_data()
  IMPLICIT NONE

  character(len=256)  :: data_dir
  type(nrowcol)       :: rowcol_count
  integer             :: ioi, N_indices
  character(6000)     :: buffer
  print'(a,t23,a)', achar(10),  '--~:| Gas and Aerosol Box Model - GABO v.0.1 |:~--'//achar(10)
  write(buffer,NML_ICOLS)
  N_INDICES = CNTNONTYPESNML(TRIM(buffer))
  IF (N_INDICES /= INM_LAST) THEN
    print FMT_FAT0, 'Number of possible input variables seems to differ from what is assumed in INM_LAST.'
    print FMT_MSG, 'Check that: in input.f90, NML_ICOLS has all indices and that INM_LAST is correct, then recompile.'
    print FMT_LEND
    print*,
    STOP
  END IF

  CALL NAME_MODS
  CALL READ_INIT_FILE
  CALL FILL_INDRELAY_WITH_INDICES
  CALL PUT_USER_SUPPLIED_TIMEOPTIONS_IN_MODELTIME


  ! ALLOCATE CONC_MAT
  IF ((ENV_file /= '') .or. (VOC_file /= '')) THEN
    IF (ENV_file /= '') THEN
      OPEN(unit=51, File=TRIM(ENV_path) //'/'//TRIM(ENV_file), STATUS='OLD')
    ELSE
      OPEN(unit=51, File=TRIM(VOC_path) //'/'//TRIM(VOC_file), STATUS='OLD')
    END IF
    ALLOCATE(CONC_MAT(ROWCOUNT(51,'#'),INM_LAST))
    ALLOCATE(TIMEVEC(ROWCOUNT(51,'#')))
    CLOSE(51)
    ! Deal with a situation where we have no input. We still need conc_mat and timevec.
  ELSE
    ALLOCATE(CONC_MAT(2,INM_LAST))
    ALLOCATE(TIMEVEC(2))
    TIMEVEC = (/0d0, MODELTIME%SIM_TIME_H/)
  END IF
  CONC_MAT = 0d0

  ! READ ENV INPUT
  if (ENV_file /= '') THEN
    OPEN(unit=51, File=TRIM(ENV_path) //'/'//TRIM(ENV_file), STATUS='OLD')
    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)
    ALLOCATE(INPUT_ENV(rowcol_count%rows,rowcol_count%cols))
    INPUT_ENV = 0
    call fill_input_buff(51,rowcol_count,INPUT_ENV,ENV_file)
    timevec = INPUT_ENV(:,1)
    CLOSE(51)
  END IF

  if (VOC_file /= '') THEN
    ! READ VOC INPUT
    OPEN(unit=51, File=TRIM(VOC_path) //'/'//TRIM(VOC_file), STATUS='OLD')
    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)
    allocate(INPUT_VOC(rowcol_count%rows,rowcol_count%cols))
    INPUT_VOC = 0
    call fill_input_buff(51,rowcol_count,INPUT_VOC,VOC_file)
    timevec = INPUT_VOC(:,1)
    CLOSE(51)
  END IF

  CALL PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_VOC,CONC_MAT)

  ! check IF dmps data is used or not. If no then do nothing
  IF (USE_DMPS) then
   write(*,FMT_SUB),'Reading DMPS file '// TRIM(DMPS_file)
   OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' //TRIM(DMPS_dir)// '/'//TRIM(DMPS_file) ,STATUS='OLD', iostat=ioi)
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
end subroutine read_input_data


subroutine READ_INIT_FILE
!-------------------------------------------------------------------------------
! Reads the init file and fills user-provided variables
!-------------------------------------------------------------------------------
  implicit none
  integer :: IOS(20)
  CALL GETARG(1,Fname_init)
  write(*,FMT_HDR), 'READING USER DEFINED INTIAL VALUES FROM: '//TRIM(ADJUSTL(Fname_init))

  IOS = 0
  OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', iostat=IOS(1))
  IF (IOS(1) /= 0) THEN
    write(*,FMT_FAT0) 'No such INIT file, exiting. Good bye.'
    write(*,FMT_LEND)
    STOP
  END IF
  READ(50,NML=NML_Path,  IOSTAT= IOS(2)) ! directories and test cases
    IF (IOS(2) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_Path, maybe some undefinded input?'
  READ(50,NML=NML_Flag,  IOSTAT= IOS(3)) ! flags
    IF (IOS(3) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_Flag, maybe some undefinded input?'
  READ(50,NML=NML_TIME,  IOSTAT= IOS(4)) ! time related stuff
    IF (IOS(4) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_TIME, maybe some undefinded input?'
  READ(50,NML=NML_DMPS,  IOSTAT= IOS(5)) ! dmps_file information
    IF (IOS(5) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_DMPS, maybe some undefinded input?'
  READ(50,NML=NML_ENV,  IOSTAT= IOS(6)) ! environmental information
    IF (IOS(6) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_Temp, maybe some undefinded input?'
  READ(50,NML=NML_VOC,   IOSTAT= IOS(7)) ! voc_file information
    IF (IOS(7) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_VOC, maybe some undefinded input?'
  READ(50,NML=NML_MODS,  IOSTAT= IOS(8)) ! modification parameters
    IF (IOS(8) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_MODS, maybe some undefinded input?'
  READ(50,NML=NML_MISC,  IOSTAT= IOS(9)) ! misc input
    IF (IOS(9) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_MISC, maybe some undefinded input?'
  READ(50,NML=NML_ICOLS,  IOSTAT= IOS(10)) ! indices input
    IF (IOS(10) /= 0) write(*,FMT_FAT0) 'Problem in init file; NML_ICOLS, maybe some undefinded input?'
  CLOSE(50)
  IF (SUM(IOS) /= 0) then
    write(*,FMT_MSG) 'Problems with the init file, exiting now. Good bye.'
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

subroutine NAME_MODS
  implicit none
  INTEGER :: i = 0
  ! Maximum length is 16 characters - for no good reason.
  i=i+1; MODS(inm_H2SO4)%NAME = 'H2SO4'
  i=i+1; MODS(inm_NH3  )%NAME = 'NH3'
  i=i+1; MODS(inm_DMA  )%NAME = 'DMA'
  i=i+1; MODS(inm_CS   )%NAME = 'condens_sink'
  i=i+1; MODS(inm_swr  )%NAME = 'shortwave_rad'
  i=i+1; MODS(inm_RH   )%NAME = 'relative_humid'
  i=i+1; MODS(inm_pres )%NAME = 'pressure'
  i=i+1; MODS(inm_temp )%NAME = 'temperature'
  i=i+1; MODS(inm_SO2  )%NAME = 'SO2'
  i=i+1; MODS(inm_NO   )%NAME = 'NO'
  i=i+1; MODS(inm_NO2  )%NAME = 'NO2'
  i=i+1; MODS(inm_CO   )%NAME = 'CO'
  i=i+1; MODS(inm_H2   )%NAME = 'H2'
  i=i+1; MODS(inm_O3   )%NAME = 'O3'
  i=i+1; MODS(inm_IPR  )%NAME = 'Ion_Prod_Rate'
  i=i+1; MODS(inm_CH3O )%NAME = 'CH3O'
  i=i+1; MODS(inm_CH3C )%NAME = 'CH3C'
  i=i+1; MODS(inm_C2H5 )%NAME = 'C2H5'
  i=i+1; MODS(inm_C5H8 )%NAME = 'C5H8'
  i=i+1; MODS(inm_MVK  )%NAME = 'MVK'
  i=i+1; MODS(inm_MEK  )%NAME = 'MEK'
  i=i+1; MODS(inm_BENZ )%NAME = 'BENZENE'
  i=i+1; MODS(inm_APIN )%NAME = 'ALPHAPINENE'
  i=i+1; MODS(inm_BPIN )%NAME = 'BETAPINENE'
  i=i+1; MODS(inm_LIMO )%NAME = 'LIMONENE'
  i=i+1; MODS(inm_Care )%NAME = 'CARENE'
  i=i+1; MODS(inm_TOLU )%NAME = 'TOLUENE'

  i=i+1; MODS(inm_RES20)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES19)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES18)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES17)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES16)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES15)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES14)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES13)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES12)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES11)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES10)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES9 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES8 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES7 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES6 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES5 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES4 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES3 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES2 )%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES1 )%NAME = 'RESERVE'
  ! As safety feature, check that here nothing is left out
  IF (i /= INM_LAST) THEN
    print FMT_FAT0, 'Trouble in NAME_MODS (input.f90). Number of input variables differs from what is assumed'
    print FMT_MSG, 'in INM_LAST. Fix and recompile.'
    print FMT_LEND
    print*,
    STOP
  END IF
end subroutine NAME_MODS

subroutine FILL_INDRELAY_WITH_INDICES
  implicit none
  integer::i
  i = 0
  INDRELAY(inm_H2SO4,1:2) = (/inm_H2SO4 , inf_H2SO4/) ;i=i+1
  INDRELAY(inm_NH3  ,1:2) = (/inm_NH3   , inf_NH3  /) ;i=i+1
  INDRELAY(inm_DMA  ,1:2) = (/inm_DMA   , inf_DMA  /) ;i=i+1
  INDRELAY(inm_CS   ,1:2) = (/inm_CS    , inf_CS   /) ;i=i+1
  INDRELAY(inm_swr  ,1:2) = (/inm_swr   , inf_swr  /) ;i=i+1
  INDRELAY(inm_RH   ,1:2) = (/inm_RH    , inf_RH   /) ;i=i+1
  INDRELAY(inm_pres ,1:2) = (/inm_pres  , inf_pres /) ;i=i+1
  INDRELAY(inm_temp ,1:2) = (/inm_temp  , inf_temp /) ;i=i+1
  INDRELAY(inm_SO2  ,1:2) = (/inm_SO2   , inf_SO2  /) ;i=i+1
  INDRELAY(inm_NO   ,1:2) = (/inm_NO    , inf_NO   /) ;i=i+1
  INDRELAY(inm_NO2  ,1:2) = (/inm_NO2   , inf_NO2  /) ;i=i+1
  INDRELAY(inm_CO   ,1:2) = (/inm_CO    , inf_CO   /) ;i=i+1
  INDRELAY(inm_H2   ,1:2) = (/inm_H2    , inf_H2   /) ;i=i+1
  INDRELAY(inm_O3   ,1:2) = (/inm_O3    , inf_O3   /) ;i=i+1
  INDRELAY(inm_IPR  ,1:2) = (/inm_IPR   , inf_IPR  /) ;i=i+1
  INDRELAY(inm_CH3O ,1:2) = (/inm_CH3O  , inf_CH3O /) ;i=i+1
  INDRELAY(inm_CH3C ,1:2) = (/inm_CH3C  , inf_CH3C /) ;i=i+1
  INDRELAY(inm_C2H5 ,1:2) = (/inm_C2H5  , inf_C2H5 /) ;i=i+1
  INDRELAY(inm_C5H8 ,1:2) = (/inm_C5H8  , inf_C5H8 /) ;i=i+1
  INDRELAY(inm_MVK  ,1:2) = (/inm_MVK   , inf_MVK  /) ;i=i+1
  INDRELAY(inm_MEK  ,1:2) = (/inm_MEK   , inf_MEK  /) ;i=i+1
  INDRELAY(inm_BENZ ,1:2) = (/inm_BENZ  , inf_BENZ /) ;i=i+1
  INDRELAY(inm_APIN ,1:2) = (/inm_APIN  , inf_APIN /) ;i=i+1
  INDRELAY(inm_BPIN ,1:2) = (/inm_BPIN  , inf_BPIN /) ;i=i+1
  INDRELAY(inm_LIMO ,1:2) = (/inm_LIMO  , inf_LIMO /) ;i=i+1
  INDRELAY(inm_Care ,1:2) = (/inm_Care  , inf_Care /) ;i=i+1
  INDRELAY(inm_TOLU ,1:2) = (/inm_TOLU  , inf_TOLU /) ;i=i+1

  INDRELAY(inm_RES20 ,1:2) = (/inm_RES20, inf_RES20/) ;i=i+1
  INDRELAY(inm_RES19 ,1:2) = (/inm_RES19, inf_RES19/) ;i=i+1
  INDRELAY(inm_RES18 ,1:2) = (/inm_RES18, inf_RES18/) ;i=i+1
  INDRELAY(inm_RES17 ,1:2) = (/inm_RES17, inf_RES17/) ;i=i+1
  INDRELAY(inm_RES16 ,1:2) = (/inm_RES16, inf_RES16/) ;i=i+1
  INDRELAY(inm_RES15 ,1:2) = (/inm_RES15, inf_RES15/) ;i=i+1
  INDRELAY(inm_RES14 ,1:2) = (/inm_RES14, inf_RES14/) ;i=i+1
  INDRELAY(inm_RES13 ,1:2) = (/inm_RES13, inf_RES13/) ;i=i+1
  INDRELAY(inm_RES12 ,1:2) = (/inm_RES12, inf_RES12/) ;i=i+1
  INDRELAY(inm_RES11 ,1:2) = (/inm_RES11, inf_RES11/) ;i=i+1
  INDRELAY(inm_RES10 ,1:2) = (/inm_RES10, inf_RES10/) ;i=i+1
  INDRELAY(inm_RES9  ,1:2) = (/inm_RES9 , inf_RES9 /) ;i=i+1
  INDRELAY(inm_RES8  ,1:2) = (/inm_RES8 , inf_RES8 /) ;i=i+1
  INDRELAY(inm_RES7  ,1:2) = (/inm_RES7 , inf_RES7 /) ;i=i+1
  INDRELAY(inm_RES6  ,1:2) = (/inm_RES6 , inf_RES6 /) ;i=i+1
  INDRELAY(inm_RES5  ,1:2) = (/inm_RES5 , inf_RES5 /) ;i=i+1
  INDRELAY(inm_RES4  ,1:2) = (/inm_RES4 , inf_RES4 /) ;i=i+1
  INDRELAY(inm_RES3  ,1:2) = (/inm_RES3 , inf_RES3 /) ;i=i+1
  INDRELAY(inm_RES2  ,1:2) = (/inm_RES2 , inf_RES2 /) ;i=i+1
  INDRELAY(inm_RES1  ,1:2) = (/inm_RES1 , inf_RES1 /) ;i=i+1

  ! As safety feature, check that here nothing is left out
  IF (i /= INM_LAST) THEN
    print FMT_FAT0, 'Trouble in FILL_INDRELAY_WITH_INDICES (input.f90). Number of input variables differs'
    print FMT_MSG, 'from what is assumed in INM_LAST. Fix and recompile.'
    print FMT_LEND
    print*,
    STOP
  END IF
end subroutine FILL_INDRELAY_WITH_INDICES

SUBROUTINE PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_VOC,CONC_MAT)
  implicit none
  REAL(DP), intent(inout) :: CONC_MAT(:,:)
  REAL(DP), intent(in)    :: INPUT_ENV(:,:)
  REAL(DP), intent(in)    :: INPUT_VOC(:,:)
  integer::i
  i = 0

  IF (inf_H2SO4 > 0) CONC_MAT(:,inm_H2SO4) = input_ENV(:,inf_H2SO4) ;i=i+1
  IF (inf_NH3   > 0) CONC_MAT(:,inm_NH3  ) = input_ENV(:,inf_NH3  ) ;i=i+1
  IF (inf_DMA   > 0) CONC_MAT(:,inm_DMA  ) = input_ENV(:,inf_DMA  ) ;i=i+1
  IF (inf_CS    > 0) CONC_MAT(:,inm_CS   ) = input_ENV(:,inf_CS   ) ;i=i+1
  IF (inf_swr   > 0) CONC_MAT(:,inm_swr  ) = input_ENV(:,inf_swr  ) ;i=i+1
  IF (inf_RH    > 0) CONC_MAT(:,inm_RH   ) = input_ENV(:,inf_RH   ) ;i=i+1
  IF (inf_pres  > 0) CONC_MAT(:,inm_pres ) = input_ENV(:,inf_pres ) ;i=i+1
  IF (inf_temp  > 0) CONC_MAT(:,inm_temp ) = input_ENV(:,inf_temp ) ;i=i+1
  IF (inf_SO2   > 0) CONC_MAT(:,inm_SO2  ) = input_ENV(:,inf_SO2  ) ;i=i+1
  IF (inf_NO    > 0) CONC_MAT(:,inm_NO   ) = input_ENV(:,inf_NO   ) ;i=i+1
  IF (inf_NO2   > 0) CONC_MAT(:,inm_NO2  ) = input_ENV(:,inf_NO2  ) ;i=i+1
  IF (inf_CO    > 0) CONC_MAT(:,inm_CO   ) = input_ENV(:,inf_CO   ) ;i=i+1
  IF (inf_H2    > 0) CONC_MAT(:,inm_H2   ) = input_ENV(:,inf_H2   ) ;i=i+1
  IF (inf_O3    > 0) CONC_MAT(:,inm_O3   ) = input_ENV(:,inf_O3   ) ;i=i+1
  IF (inf_IPR   > 0) CONC_MAT(:,inm_IPR  ) = input_ENV(:,inf_IPR  ) ;i=i+1

  IF (inf_CH3O  > 0) CONC_MAT(:,inm_CH3O ) = input_VOC(:,inf_CH3O ) ;i=i+1
  IF (inf_CH3C  > 0) CONC_MAT(:,inm_CH3C ) = input_VOC(:,inf_CH3C ) ;i=i+1
  IF (inf_C2H5  > 0) CONC_MAT(:,inm_C2H5 ) = input_VOC(:,inf_C2H5 ) ;i=i+1
  IF (inf_C5H8  > 0) CONC_MAT(:,inm_C5H8 ) = input_VOC(:,inf_C5H8 ) ;i=i+1
  IF (inf_MVK   > 0) CONC_MAT(:,inm_MVK  ) = input_VOC(:,inf_MVK  ) ;i=i+1
  IF (inf_MEK   > 0) CONC_MAT(:,inm_MEK  ) = input_VOC(:,inf_MEK  ) ;i=i+1
  IF (inf_BENZ  > 0) CONC_MAT(:,inm_BENZ ) = input_VOC(:,inf_BENZ ) ;i=i+1
  IF (inf_APIN  > 0) CONC_MAT(:,inm_APIN ) = input_VOC(:,inf_APIN ) ;i=i+1
  IF (inf_BPIN  > 0) CONC_MAT(:,inm_BPIN ) = input_VOC(:,inf_BPIN ) ;i=i+1
  IF (inf_LIMO  > 0) CONC_MAT(:,inm_LIMO ) = input_VOC(:,inf_LIMO ) ;i=i+1
  IF (inf_Care  > 0) CONC_MAT(:,inm_Care ) = input_VOC(:,inf_Care ) ;i=i+1
  IF (inf_TOLU  > 0) CONC_MAT(:,inm_TOLU ) = input_VOC(:,inf_TOLU ) ;i=i+1

  IF (inf_RES20 > 0) CONC_MAT(:,inm_RES20) = input_VOC(:,inf_RES20) ;i=i+1
  IF (inf_RES19 > 0) CONC_MAT(:,inm_RES19) = input_VOC(:,inf_RES19) ;i=i+1
  IF (inf_RES18 > 0) CONC_MAT(:,inm_RES18) = input_VOC(:,inf_RES18) ;i=i+1
  IF (inf_RES17 > 0) CONC_MAT(:,inm_RES17) = input_VOC(:,inf_RES17) ;i=i+1
  IF (inf_RES16 > 0) CONC_MAT(:,inm_RES16) = input_VOC(:,inf_RES16) ;i=i+1
  IF (inf_RES15 > 0) CONC_MAT(:,inm_RES15) = input_VOC(:,inf_RES15) ;i=i+1
  IF (inf_RES14 > 0) CONC_MAT(:,inm_RES14) = input_VOC(:,inf_RES14) ;i=i+1
  IF (inf_RES13 > 0) CONC_MAT(:,inm_RES13) = input_VOC(:,inf_RES13) ;i=i+1
  IF (inf_RES12 > 0) CONC_MAT(:,inm_RES12) = input_VOC(:,inf_RES12) ;i=i+1
  IF (inf_RES11 > 0) CONC_MAT(:,inm_RES11) = input_VOC(:,inf_RES11) ;i=i+1
  IF (inf_RES10 > 0) CONC_MAT(:,inm_RES10) = input_VOC(:,inf_RES10) ;i=i+1
  IF (inf_RES9  > 0) CONC_MAT(:,inm_RES9 ) = input_VOC(:,inf_RES9 ) ;i=i+1
  IF (inf_RES8  > 0) CONC_MAT(:,inm_RES8 ) = input_VOC(:,inf_RES8 ) ;i=i+1
  IF (inf_RES7  > 0) CONC_MAT(:,inm_RES7 ) = input_VOC(:,inf_RES7 ) ;i=i+1
  IF (inf_RES6  > 0) CONC_MAT(:,inm_RES6 ) = input_VOC(:,inf_RES6 ) ;i=i+1
  IF (inf_RES5  > 0) CONC_MAT(:,inm_RES5 ) = input_VOC(:,inf_RES5 ) ;i=i+1
  IF (inf_RES4  > 0) CONC_MAT(:,inm_RES4 ) = input_VOC(:,inf_RES4 ) ;i=i+1
  IF (inf_RES3  > 0) CONC_MAT(:,inm_RES3 ) = input_VOC(:,inf_RES3 ) ;i=i+1
  IF (inf_RES2  > 0) CONC_MAT(:,inm_RES2 ) = input_VOC(:,inf_RES2 ) ;i=i+1
  IF (inf_RES1  > 0) CONC_MAT(:,inm_RES1 ) = input_VOC(:,inf_RES1 ) ;i=i+1


  ! As safety feature, check that here nothing is left out
  IF (i /= INM_LAST) THEN
    print FMT_FAT0, 'Trouble in PUT_INPUT_IN_THEIR_PLACES (input.f90). Number of input variables differs'
    print FMT_MSG, 'from what is assumed in INM_LAST. Fix and recompile.'
    print FMT_LEND
    print*,
    STOP
  END IF

END SUBROUTINE PUT_INPUT_IN_THEIR_PLACES

subroutine fill_input_buff(unit,rowcol_count, INPUT_BF,Input_file)
  implicit none
  type(nrowcol), intent(in) :: rowcol_count
  real(dp), intent(inout) :: INPUT_BF(:,:)
  character(*) :: Input_file
  integer :: i,j,k, ioi, unit
  ! Reading data into the variable
  i = 1
  DO k = 1, rowcol_count%rows
    READ(unit,*, iostat=ioi) (INPUT_BF(i,j),j=1,rowcol_count%cols)
    IF ((ioi /= 0) .and. (i==1)) THEN
      print FMT_SUB, 'Header row omitted from file "'// TRIM(Input_file) //'".'
    ELSE IF ((ioi /= 0) .and. (i>1)) THEN
      print FMT_WARN1, 'Bad value in file '// TRIM(Input_file) //'". Maybe an illegal value in line #', real(k)
    ELSE
      i=i+1
    END IF
  END DO
end subroutine fill_input_buff

end module INPUT
