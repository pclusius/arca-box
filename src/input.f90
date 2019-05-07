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
INTEGER :: inf_temp  = -1 ; INTEGER, parameter :: inm_temp  = 1
INTEGER :: inf_pres  = -1 ; INTEGER, parameter :: inm_pres  = 2
INTEGER :: inf_RH    = -1 ; INTEGER, parameter :: inm_RH    = 3
INTEGER :: inf_CS    = -1 ; INTEGER, parameter :: inm_CS    = 4
INTEGER :: inf_swr   = -1 ; INTEGER, parameter :: inm_swr   = 5
INTEGER :: inf_IPR   = -1 ; INTEGER, parameter :: inm_IPR   = 6
! Inorganics
INTEGER :: inf_H2SO4 = -1 ; INTEGER, parameter :: inm_H2SO4 = 7
INTEGER :: inf_NH3   = -1 ; INTEGER, parameter :: inm_NH3   = 8
INTEGER :: inf_DMA   = -1 ; INTEGER, parameter :: inm_DMA   = 9
INTEGER :: inf_SO2   = -1 ; INTEGER, parameter :: inm_SO2   = 10
INTEGER :: inf_NO    = -1 ; INTEGER, parameter :: inm_NO    = 11
INTEGER :: inf_NO2   = -1 ; INTEGER, parameter :: inm_NO2   = 12
INTEGER :: inf_CO    = -1 ; INTEGER, parameter :: inm_CO    = 13
INTEGER :: inf_H2    = -1 ; INTEGER, parameter :: inm_H2    = 14
INTEGER :: inf_O3    = -1 ; INTEGER, parameter :: inm_O3    = 15
INTEGER :: inf_RES16 = -1 ; INTEGER, parameter :: inm_RES16 = 16
INTEGER :: inf_RES17 = -1 ; INTEGER, parameter :: inm_RES17 = 17
INTEGER :: inf_RES18 = -1 ; INTEGER, parameter :: inm_RES18 = 18
INTEGER :: inf_RES19 = -1 ; INTEGER, parameter :: inm_RES19 = 19

! The following will be provided in MCM_file
INTEGER :: INF_CH3OH        = -1 ; INTEGER, parameter :: INM_CH3OH      = 20
INTEGER :: INF_C2H5OH       = -1 ; INTEGER, parameter :: INM_C2H5OH     = 21
INTEGER :: INF_NPROPOL      = -1 ; INTEGER, parameter :: INM_NPROPOL    = 22
INTEGER :: INF_IPROPOL      = -1 ; INTEGER, parameter :: INM_IPROPOL    = 23
INTEGER :: INF_NBUTOL       = -1 ; INTEGER, parameter :: INM_NBUTOL     = 24
INTEGER :: INF_BUT2OL       = -1 ; INTEGER, parameter :: INM_BUT2OL     = 25
INTEGER :: INF_IBUTOL       = -1 ; INTEGER, parameter :: INM_IBUTOL     = 26
INTEGER :: INF_TBUTOL       = -1 ; INTEGER, parameter :: INM_TBUTOL     = 27
INTEGER :: INF_PECOH        = -1 ; INTEGER, parameter :: INM_PECOH      = 28
INTEGER :: INF_IPEAOH       = -1 ; INTEGER, parameter :: INM_IPEAOH     = 29
INTEGER :: INF_ME3BUOL      = -1 ; INTEGER, parameter :: INM_ME3BUOL    = 30
INTEGER :: INF_IPECOH       = -1 ; INTEGER, parameter :: INM_IPECOH     = 31
INTEGER :: INF_IPEBOH       = -1 ; INTEGER, parameter :: INM_IPEBOH     = 32
INTEGER :: INF_CYHEXOL      = -1 ; INTEGER, parameter :: INM_CYHEXOL    = 33
INTEGER :: INF_MIBKAOH      = -1 ; INTEGER, parameter :: INM_MIBKAOH    = 34
INTEGER :: INF_ETHGLY       = -1 ; INTEGER, parameter :: INM_ETHGLY     = 35
INTEGER :: INF_PROPGLY      = -1 ; INTEGER, parameter :: INM_PROPGLY    = 36
INTEGER :: INF_MBO          = -1 ; INTEGER, parameter :: INM_MBO        = 37
INTEGER :: INF_HCHO         = -1 ; INTEGER, parameter :: INM_HCHO       = 38
INTEGER :: INF_CH3CHO       = -1 ; INTEGER, parameter :: INM_CH3CHO     = 39
INTEGER :: INF_C2H5CHO      = -1 ; INTEGER, parameter :: INM_C2H5CHO    = 40
INTEGER :: INF_C3H7CHO      = -1 ; INTEGER, parameter :: INM_C3H7CHO    = 41
INTEGER :: INF_IPRCHO       = -1 ; INTEGER, parameter :: INM_IPRCHO     = 42
INTEGER :: INF_C4H9CHO      = -1 ; INTEGER, parameter :: INM_C4H9CHO    = 43
INTEGER :: INF_ACR          = -1 ; INTEGER, parameter :: INM_ACR        = 44
INTEGER :: INF_MACR         = -1 ; INTEGER, parameter :: INM_MACR       = 45
INTEGER :: INF_C4ALDB       = -1 ; INTEGER, parameter :: INM_C4ALDB     = 46
INTEGER :: INF_CH4          = -1 ; INTEGER, parameter :: INM_CH4        = 47
INTEGER :: INF_C2H6         = -1 ; INTEGER, parameter :: INM_C2H6       = 48
INTEGER :: INF_C3H8         = -1 ; INTEGER, parameter :: INM_C3H8       = 49
INTEGER :: INF_NC4H10       = -1 ; INTEGER, parameter :: INM_NC4H10     = 50
INTEGER :: INF_IC4H10       = -1 ; INTEGER, parameter :: INM_IC4H10     = 51
INTEGER :: INF_NC5H12       = -1 ; INTEGER, parameter :: INM_NC5H12     = 52
INTEGER :: INF_IC5H12       = -1 ; INTEGER, parameter :: INM_IC5H12     = 53
INTEGER :: INF_NEOP         = -1 ; INTEGER, parameter :: INM_NEOP       = 54
INTEGER :: INF_NC6H14       = -1 ; INTEGER, parameter :: INM_NC6H14     = 55
INTEGER :: INF_M2PE         = -1 ; INTEGER, parameter :: INM_M2PE       = 56
INTEGER :: INF_M3PE         = -1 ; INTEGER, parameter :: INM_M3PE       = 57
INTEGER :: INF_M22C4        = -1 ; INTEGER, parameter :: INM_M22C4      = 58
INTEGER :: INF_M23C4        = -1 ; INTEGER, parameter :: INM_M23C4      = 59
INTEGER :: INF_NC7H16       = -1 ; INTEGER, parameter :: INM_NC7H16     = 60
INTEGER :: INF_M2HEX        = -1 ; INTEGER, parameter :: INM_M2HEX      = 61
INTEGER :: INF_M3HEX        = -1 ; INTEGER, parameter :: INM_M3HEX      = 62
INTEGER :: INF_NC8H18       = -1 ; INTEGER, parameter :: INM_NC8H18     = 63
INTEGER :: INF_NC9H20       = -1 ; INTEGER, parameter :: INM_NC9H20     = 64
INTEGER :: INF_NC10H22      = -1 ; INTEGER, parameter :: INM_NC10H22    = 65
INTEGER :: INF_NC11H24      = -1 ; INTEGER, parameter :: INM_NC11H24    = 66
INTEGER :: INF_NC12H26      = -1 ; INTEGER, parameter :: INM_NC12H26    = 67
INTEGER :: INF_CHEX         = -1 ; INTEGER, parameter :: INM_CHEX       = 68
INTEGER :: INF_C2H4         = -1 ; INTEGER, parameter :: INM_C2H4       = 69
INTEGER :: INF_C3H6         = -1 ; INTEGER, parameter :: INM_C3H6       = 70
INTEGER :: INF_BUT1ENE      = -1 ; INTEGER, parameter :: INM_BUT1ENE    = 71
INTEGER :: INF_CBUT2ENE     = -1 ; INTEGER, parameter :: INM_CBUT2ENE   = 72
INTEGER :: INF_TBUT2ENE     = -1 ; INTEGER, parameter :: INM_TBUT2ENE   = 73
INTEGER :: INF_MEPROPENE    = -1 ; INTEGER, parameter :: INM_MEPROPENE  = 74
INTEGER :: INF_PENT1ENE     = -1 ; INTEGER, parameter :: INM_PENT1ENE   = 75
INTEGER :: INF_CPENT2ENE    = -1 ; INTEGER, parameter :: INM_CPENT2ENE  = 76
INTEGER :: INF_TPENT2ENE    = -1 ; INTEGER, parameter :: INM_TPENT2ENE  = 77
INTEGER :: INF_ME2BUT1ENE   = -1 ; INTEGER, parameter :: INM_ME2BUT1ENE = 78
INTEGER :: INF_ME3BUT1ENE   = -1 ; INTEGER, parameter :: INM_ME3BUT1ENE = 79
INTEGER :: INF_ME2BUT2ENE   = -1 ; INTEGER, parameter :: INM_ME2BUT2ENE = 80
INTEGER :: INF_HEX1ENE      = -1 ; INTEGER, parameter :: INM_HEX1ENE    = 81
INTEGER :: INF_CHEX2ENE     = -1 ; INTEGER, parameter :: INM_CHEX2ENE   = 82
INTEGER :: INF_THEX2ENE     = -1 ; INTEGER, parameter :: INM_THEX2ENE   = 83
INTEGER :: INF_DM23BU2ENE   = -1 ; INTEGER, parameter :: INM_DM23BU2ENE = 84
INTEGER :: INF_C2H2         = -1 ; INTEGER, parameter :: INM_C2H2       = 85
INTEGER :: INF_BENZENE      = -1 ; INTEGER, parameter :: INM_BENZENE    = 86
INTEGER :: INF_TOLUENE      = -1 ; INTEGER, parameter :: INM_TOLUENE    = 87
INTEGER :: INF_OXYL         = -1 ; INTEGER, parameter :: INM_OXYL       = 88
INTEGER :: INF_MXYL         = -1 ; INTEGER, parameter :: INM_MXYL       = 89
INTEGER :: INF_PXYL         = -1 ; INTEGER, parameter :: INM_PXYL       = 90
INTEGER :: INF_EBENZ        = -1 ; INTEGER, parameter :: INM_EBENZ      = 91
INTEGER :: INF_PBENZ        = -1 ; INTEGER, parameter :: INM_PBENZ      = 92
INTEGER :: INF_IPBENZ       = -1 ; INTEGER, parameter :: INM_IPBENZ     = 93
INTEGER :: INF_TM123B       = -1 ; INTEGER, parameter :: INM_TM123B     = 94
INTEGER :: INF_TM124B       = -1 ; INTEGER, parameter :: INM_TM124B     = 95
INTEGER :: INF_TM135B       = -1 ; INTEGER, parameter :: INM_TM135B     = 96
INTEGER :: INF_OETHTOL      = -1 ; INTEGER, parameter :: INM_OETHTOL    = 97
INTEGER :: INF_METHTOL      = -1 ; INTEGER, parameter :: INM_METHTOL    = 98
INTEGER :: INF_PETHTOL      = -1 ; INTEGER, parameter :: INM_PETHTOL    = 99
INTEGER :: INF_DIME35EB     = -1 ; INTEGER, parameter :: INM_DIME35EB   = 100
INTEGER :: INF_DIET35TOL    = -1 ; INTEGER, parameter :: INM_DIET35TOL  = 101
INTEGER :: INF_STYRENE      = -1 ; INTEGER, parameter :: INM_STYRENE    = 102
INTEGER :: INF_BENZAL       = -1 ; INTEGER, parameter :: INM_BENZAL     = 103
INTEGER :: INF_CH3CL        = -1 ; INTEGER, parameter :: INM_CH3CL      = 104
INTEGER :: INF_CH2CL2       = -1 ; INTEGER, parameter :: INM_CH2CL2     = 105
INTEGER :: INF_CHCL3        = -1 ; INTEGER, parameter :: INM_CHCL3      = 106
INTEGER :: INF_CH3CCL3      = -1 ; INTEGER, parameter :: INM_CH3CCL3    = 107
INTEGER :: INF_TCE          = -1 ; INTEGER, parameter :: INM_TCE        = 108
INTEGER :: INF_TRICLETH     = -1 ; INTEGER, parameter :: INM_TRICLETH   = 109
INTEGER :: INF_CDICLETH     = -1 ; INTEGER, parameter :: INM_CDICLETH   = 110
INTEGER :: INF_TDICLETH     = -1 ; INTEGER, parameter :: INM_TDICLETH   = 111
INTEGER :: INF_CH2CLCH2CL   = -1 ; INTEGER, parameter :: INM_CH2CLCH2CL = 112
INTEGER :: INF_CCL2CH2      = -1 ; INTEGER, parameter :: INM_CCL2CH2    = 113
INTEGER :: INF_CL12PROP     = -1 ; INTEGER, parameter :: INM_CL12PROP   = 114
INTEGER :: INF_CHCL2CH3     = -1 ; INTEGER, parameter :: INM_CHCL2CH3   = 115
INTEGER :: INF_CH3CH2CL     = -1 ; INTEGER, parameter :: INM_CH3CH2CL   = 116
INTEGER :: INF_CHCL2CHCL2   = -1 ; INTEGER, parameter :: INM_CHCL2CHCL2 = 117
INTEGER :: INF_CH2CLCHCL2   = -1 ; INTEGER, parameter :: INM_CH2CLCHCL2 = 118
INTEGER :: INF_VINCL        = -1 ; INTEGER, parameter :: INM_VINCL      = 119
INTEGER :: INF_C4H6         = -1 ; INTEGER, parameter :: INM_C4H6       = 120
INTEGER :: INF_C5H8         = -1 ; INTEGER, parameter :: INM_C5H8       = 121
INTEGER :: INF_CH3OCHO      = -1 ; INTEGER, parameter :: INM_CH3OCHO    = 122
INTEGER :: INF_METHACET     = -1 ; INTEGER, parameter :: INM_METHACET   = 123
INTEGER :: INF_ETHACET      = -1 ; INTEGER, parameter :: INM_ETHACET    = 124
INTEGER :: INF_NPROACET     = -1 ; INTEGER, parameter :: INM_NPROACET   = 125
INTEGER :: INF_IPROACET     = -1 ; INTEGER, parameter :: INM_IPROACET   = 126
INTEGER :: INF_NBUTACET     = -1 ; INTEGER, parameter :: INM_NBUTACET   = 127
INTEGER :: INF_SBUTACET     = -1 ; INTEGER, parameter :: INM_SBUTACET   = 128
INTEGER :: INF_TBUACET      = -1 ; INTEGER, parameter :: INM_TBUACET    = 129
INTEGER :: INF_CH3OCH3      = -1 ; INTEGER, parameter :: INM_CH3OCH3    = 130
INTEGER :: INF_DIETETHER    = -1 ; INTEGER, parameter :: INM_DIETETHER  = 131
INTEGER :: INF_MTBE         = -1 ; INTEGER, parameter :: INM_MTBE       = 132
INTEGER :: INF_DIIPRETHER   = -1 ; INTEGER, parameter :: INM_DIIPRETHER = 133
INTEGER :: INF_ETBE         = -1 ; INTEGER, parameter :: INM_ETBE       = 134
INTEGER :: INF_MO2EOL       = -1 ; INTEGER, parameter :: INM_MO2EOL     = 135
INTEGER :: INF_EOX2EOL      = -1 ; INTEGER, parameter :: INM_EOX2EOL    = 136
INTEGER :: INF_PR2OHMOX     = -1 ; INTEGER, parameter :: INM_PR2OHMOX   = 137
INTEGER :: INF_BUOX2ETOH    = -1 ; INTEGER, parameter :: INM_BUOX2ETOH  = 138
INTEGER :: INF_BOX2PROL     = -1 ; INTEGER, parameter :: INM_BOX2PROL   = 139
INTEGER :: INF_CH3BR        = -1 ; INTEGER, parameter :: INM_CH3BR      = 140
INTEGER :: INF_DIBRET       = -1 ; INTEGER, parameter :: INM_DIBRET     = 141
INTEGER :: INF_CH3COCH3     = -1 ; INTEGER, parameter :: INM_CH3COCH3   = 142
INTEGER :: INF_MEK          = -1 ; INTEGER, parameter :: INM_MEK        = 143
INTEGER :: INF_MPRK         = -1 ; INTEGER, parameter :: INM_MPRK       = 144
INTEGER :: INF_DIEK         = -1 ; INTEGER, parameter :: INM_DIEK       = 145
INTEGER :: INF_MIPK         = -1 ; INTEGER, parameter :: INM_MIPK       = 146
INTEGER :: INF_HEX2ONE      = -1 ; INTEGER, parameter :: INM_HEX2ONE    = 147
INTEGER :: INF_HEX3ONE      = -1 ; INTEGER, parameter :: INM_HEX3ONE    = 148
INTEGER :: INF_MIBK         = -1 ; INTEGER, parameter :: INM_MIBK       = 149
INTEGER :: INF_MTBK         = -1 ; INTEGER, parameter :: INM_MTBK       = 150
INTEGER :: INF_CYHEXONE     = -1 ; INTEGER, parameter :: INM_CYHEXONE   = 151
INTEGER :: INF_APINENE      = -1 ; INTEGER, parameter :: INM_APINENE    = 152
INTEGER :: INF_BPINENE      = -1 ; INTEGER, parameter :: INM_BPINENE    = 153
INTEGER :: INF_LIMONENE     = -1 ; INTEGER, parameter :: INM_LIMONENE   = 154
INTEGER :: INF_BCARY        = -1 ; INTEGER, parameter :: INM_BCARY      = 155
INTEGER :: INF_HCOOH        = -1 ; INTEGER, parameter :: INM_HCOOH      = 156
INTEGER :: INF_CH3CO2H      = -1 ; INTEGER, parameter :: INM_CH3CO2H    = 157
INTEGER :: INF_PROPACID     = -1 ; INTEGER, parameter :: INM_PROPACID   = 158
INTEGER :: INF_DMM          = -1 ; INTEGER, parameter :: INM_DMM        = 159
INTEGER :: INF_DMC          = -1 ; INTEGER, parameter :: INM_DMC        = 160
INTEGER :: INF_DMS          = -1 ; INTEGER, parameter :: INM_DMS        = 161
INTEGER :: INF_ETHOX        = -1 ; INTEGER, parameter :: INM_ETHOX      = 162

! Some reserves for EXTRAS

INTEGER, parameter :: inm_LAST  = 162 ! <-whatever the last index on previous line is

! INDICES
NAMELIST /NML_ICOLS/ inf_H2SO4,inf_NH3,inf_DMA,inf_CS,inf_swr,&
inf_RH,inf_pres,inf_temp,inf_SO2,inf_NO,inf_NO2,inf_CO,inf_H2,inf_O3,inf_IPR,&
inf_RES16,inf_RES17,inf_RES18,inf_RES19,&

INF_CH3OH,INF_C2H5OH,INF_NPROPOL,INF_IPROPOL,INF_NBUTOL,INF_BUT2OL,&
INF_IBUTOL,INF_TBUTOL,INF_PECOH,INF_IPEAOH,INF_ME3BUOL,INF_IPECOH,&
INF_IPEBOH,INF_CYHEXOL,INF_MIBKAOH,INF_ETHGLY,INF_PROPGLY,INF_MBO,&
INF_HCHO,INF_CH3CHO,INF_C2H5CHO,INF_C3H7CHO,INF_IPRCHO,INF_C4H9CHO,&
INF_ACR,INF_MACR,INF_C4ALDB,INF_CH4,INF_C2H6,INF_C3H8,INF_NC4H10,&
INF_IC4H10,INF_NC5H12,INF_IC5H12,INF_NEOP,INF_NC6H14,INF_M2PE,&
INF_M3PE,INF_M22C4,INF_M23C4,INF_NC7H16,INF_M2HEX,INF_M3HEX,&
INF_NC8H18,INF_NC9H20,INF_NC10H22,INF_NC11H24,INF_NC12H26,INF_CHEX,&
INF_C2H4,INF_C3H6,INF_BUT1ENE,INF_CBUT2ENE,INF_TBUT2ENE,INF_MEPROPENE,&
INF_PENT1ENE,INF_CPENT2ENE,INF_TPENT2ENE,INF_ME2BUT1ENE,INF_ME3BUT1ENE,&
INF_ME2BUT2ENE,INF_HEX1ENE,INF_CHEX2ENE,INF_THEX2ENE,INF_DM23BU2ENE,&
INF_C2H2,INF_BENZENE,INF_TOLUENE,INF_OXYL,INF_MXYL,INF_PXYL,INF_EBENZ,&
INF_PBENZ,INF_IPBENZ,INF_TM123B,INF_TM124B,INF_TM135B,INF_OETHTOL,&
INF_METHTOL,INF_PETHTOL,INF_DIME35EB,INF_DIET35TOL,INF_STYRENE,INF_BENZAL,&
INF_CH3CL,INF_CH2CL2,INF_CHCL3,INF_CH3CCL3,INF_TCE,INF_TRICLETH,INF_CDICLETH,&
INF_TDICLETH,INF_CH2CLCH2CL,INF_CCL2CH2,INF_CL12PROP,INF_CHCL2CH3,&
INF_CH3CH2CL,INF_CHCL2CHCL2,INF_CH2CLCHCL2,INF_VINCL,INF_C4H6,INF_C5H8,&
INF_CH3OCHO,INF_METHACET,INF_ETHACET,INF_NPROACET,INF_IPROACET,INF_NBUTACET,&
INF_SBUTACET,INF_TBUACET,INF_CH3OCH3,INF_DIETETHER,INF_MTBE,INF_DIIPRETHER,&
INF_ETBE,INF_MO2EOL,INF_EOX2EOL,INF_PR2OHMOX,INF_BUOX2ETOH,INF_BOX2PROL,&
INF_CH3BR,INF_DIBRET,INF_CH3COCH3,INF_MEK,INF_MPRK,INF_DIEK,INF_MIPK,&
INF_HEX2ONE,INF_HEX3ONE,INF_MIBK,INF_MTBK,INF_CYHEXONE,INF_APINENE,&
INF_BPINENE,INF_LIMONENE,INF_BCARY,INF_HCOOH,INF_CH3CO2H,INF_PROPACID,&
INF_DMM,INF_DMC,INF_DMS,INF_ETHOX

INTEGER :: INDRELAY(inm_LAST), INDRELAY_CH(inm_LAST) = 0

REAL(dp), allocatable, private :: INPUT_ENV(:,:)  ! will be of shape ( len(timevec) : inm_last+1 )
REAL(dp), allocatable, private :: INPUT_MCM(:,:)  ! will be of shape ( len(timevec) : inm_last+1 )
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

! MCM INPUT
character(len=256)  :: MCM_path = ''
character(len=256)  :: MCM_file = ''
NAMELIST /NML_MCM / MCM_path, MCM_file

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
  IF ((ENV_file /= '') .or. (MCM_file /= '')) THEN
    IF (ENV_file /= '') THEN
      OPEN(unit=51, File=TRIM(ENV_path) //'/'//TRIM(ENV_file), STATUS='OLD')
    ELSE
      OPEN(unit=51, File=TRIM(MCM_path) //'/'//TRIM(MCM_file), STATUS='OLD')
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

  if (MCM_file /= '') THEN
    ! READ MCM INPUT
    OPEN(unit=51, File=TRIM(MCM_path) //'/'//TRIM(MCM_file), STATUS='OLD')
    rowcol_count%rows = ROWCOUNT(51,'#')
    rowcol_count%cols = COLCOUNT(51)
    allocate(INPUT_MCM(rowcol_count%rows,rowcol_count%cols))
    INPUT_MCM = 0
    call fill_input_buff(51,rowcol_count,INPUT_MCM,MCM_file)
    timevec = INPUT_MCM(:,1)
    CLOSE(51)
  END IF
  CALL PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)


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
    IF (IOS(2) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Path, maybe some undefinded input?'
  READ(50,NML=NML_Flag,  IOSTAT= IOS(3)) ! flags
    IF (IOS(3) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Flag, maybe some undefinded input?'
  READ(50,NML=NML_TIME,  IOSTAT= IOS(4)) ! time related stuff
    IF (IOS(4) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_TIME, maybe some undefinded input?'
  READ(50,NML=NML_DMPS,  IOSTAT= IOS(5)) ! dmps_file information
    IF (IOS(5) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_DMPS, maybe some undefinded input?'
  READ(50,NML=NML_ENV,  IOSTAT= IOS(6)) ! environmental information
    IF (IOS(6) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_Temp, maybe some undefinded input?'
  READ(50,NML=NML_MCM,   IOSTAT= IOS(7)) ! MCM_file information
    IF (IOS(7) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MCM, maybe some undefinded input?'
  READ(50,NML=NML_MODS,  IOSTAT= IOS(8)) ! modification parameters
    IF (IOS(8) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MODS, maybe some undefinded input?'
  READ(50,NML=NML_MISC,  IOSTAT= IOS(9)) ! misc input
    IF (IOS(9) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_MISC, maybe some undefinded input?'
  READ(50,NML=NML_ICOLS,  IOSTAT= IOS(10)) ! indices input
    IF (IOS(10) /= 0) write(*,FMT_FAT0) 'Problem in INITFILE; NML_ICOLS, maybe some undefinded input?'
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

  i=i+1; MODS(inm_temp )%NAME = 'temperature'
  i=i+1; MODS(inm_pres )%NAME = 'pressure'
  i=i+1; MODS(inm_RH   )%NAME = 'relative_humid'
  i=i+1; MODS(inm_CS   )%NAME = 'condens_sink'
  i=i+1; MODS(inm_swr  )%NAME = 'shortwave_rad'
  i=i+1; MODS(inm_IPR  )%NAME = 'Ion_Prod_Rate'

  i=i+1; MODS(inm_H2SO4)%NAME = 'H2SO4'
  i=i+1; MODS(inm_NH3  )%NAME = 'NH3'
  i=i+1; MODS(inm_DMA  )%NAME = 'DMA'
  i=i+1; MODS(inm_SO2  )%NAME = 'SO2'
  i=i+1; MODS(inm_NO   )%NAME = 'NO'
  i=i+1; MODS(inm_NO2  )%NAME = 'NO2'
  i=i+1; MODS(inm_CO   )%NAME = 'CO'
  i=i+1; MODS(inm_H2   )%NAME = 'H2'
  i=i+1; MODS(inm_O3   )%NAME = 'O3'

  i=i+1; MODS(inm_RES16)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES17)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES18)%NAME = 'RESERVE'
  i=i+1; MODS(inm_RES19)%NAME = 'RESERVE'

  i=i+1; MODS(INM_CH3OH        )%NAME = "CH3OH"
  i=i+1; MODS(INM_C2H5OH       )%NAME = "C2H5OH"
  i=i+1; MODS(INM_NPROPOL      )%NAME = "NPROPOL"
  i=i+1; MODS(INM_IPROPOL      )%NAME = "IPROPOL"
  i=i+1; MODS(INM_NBUTOL       )%NAME = "NBUTOL"
  i=i+1; MODS(INM_BUT2OL       )%NAME = "BUT2OL"
  i=i+1; MODS(INM_IBUTOL       )%NAME = "IBUTOL"
  i=i+1; MODS(INM_TBUTOL       )%NAME = "TBUTOL"
  i=i+1; MODS(INM_PECOH        )%NAME = "PECOH"
  i=i+1; MODS(INM_IPEAOH       )%NAME = "IPEAOH"
  i=i+1; MODS(INM_ME3BUOL      )%NAME = "ME3BUOL"
  i=i+1; MODS(INM_IPECOH       )%NAME = "IPECOH"
  i=i+1; MODS(INM_IPEBOH       )%NAME = "IPEBOH"
  i=i+1; MODS(INM_CYHEXOL      )%NAME = "CYHEXOL"
  i=i+1; MODS(INM_MIBKAOH      )%NAME = "MIBKAOH"
  i=i+1; MODS(INM_ETHGLY       )%NAME = "ETHGLY"
  i=i+1; MODS(INM_PROPGLY      )%NAME = "PROPGLY"
  i=i+1; MODS(INM_MBO          )%NAME = "MBO"
  i=i+1; MODS(INM_HCHO         )%NAME = "HCHO"
  i=i+1; MODS(INM_CH3CHO       )%NAME = "CH3CHO"
  i=i+1; MODS(INM_C2H5CHO      )%NAME = "C2H5CHO"
  i=i+1; MODS(INM_C3H7CHO      )%NAME = "C3H7CHO"
  i=i+1; MODS(INM_IPRCHO       )%NAME = "IPRCHO"
  i=i+1; MODS(INM_C4H9CHO      )%NAME = "C4H9CHO"
  i=i+1; MODS(INM_ACR          )%NAME = "ACR"
  i=i+1; MODS(INM_MACR         )%NAME = "MACR"
  i=i+1; MODS(INM_C4ALDB       )%NAME = "C4ALDB"
  i=i+1; MODS(INM_CH4          )%NAME = "CH4"
  i=i+1; MODS(INM_C2H6         )%NAME = "C2H6"
  i=i+1; MODS(INM_C3H8         )%NAME = "C3H8"
  i=i+1; MODS(INM_NC4H10       )%NAME = "NC4H10"
  i=i+1; MODS(INM_IC4H10       )%NAME = "IC4H10"
  i=i+1; MODS(INM_NC5H12       )%NAME = "NC5H12"
  i=i+1; MODS(INM_IC5H12       )%NAME = "IC5H12"
  i=i+1; MODS(INM_NEOP         )%NAME = "NEOP"
  i=i+1; MODS(INM_NC6H14       )%NAME = "NC6H14"
  i=i+1; MODS(INM_M2PE         )%NAME = "M2PE"
  i=i+1; MODS(INM_M3PE         )%NAME = "M3PE"
  i=i+1; MODS(INM_M22C4        )%NAME = "M22C4"
  i=i+1; MODS(INM_M23C4        )%NAME = "M23C4"
  i=i+1; MODS(INM_NC7H16       )%NAME = "NC7H16"
  i=i+1; MODS(INM_M2HEX        )%NAME = "M2HEX"
  i=i+1; MODS(INM_M3HEX        )%NAME = "M3HEX"
  i=i+1; MODS(INM_NC8H18       )%NAME = "NC8H18"
  i=i+1; MODS(INM_NC9H20       )%NAME = "NC9H20"
  i=i+1; MODS(INM_NC10H22      )%NAME = "NC10H22"
  i=i+1; MODS(INM_NC11H24      )%NAME = "NC11H24"
  i=i+1; MODS(INM_NC12H26      )%NAME = "NC12H26"
  i=i+1; MODS(INM_CHEX         )%NAME = "CHEX"
  i=i+1; MODS(INM_C2H4         )%NAME = "C2H4"
  i=i+1; MODS(INM_C3H6         )%NAME = "C3H6"
  i=i+1; MODS(INM_BUT1ENE      )%NAME = "BUT1ENE"
  i=i+1; MODS(INM_CBUT2ENE     )%NAME = "CBUT2ENE"
  i=i+1; MODS(INM_TBUT2ENE     )%NAME = "TBUT2ENE"
  i=i+1; MODS(INM_MEPROPENE    )%NAME = "MEPROPENE"
  i=i+1; MODS(INM_PENT1ENE     )%NAME = "PENT1ENE"
  i=i+1; MODS(INM_CPENT2ENE    )%NAME = "CPENT2ENE"
  i=i+1; MODS(INM_TPENT2ENE    )%NAME = "TPENT2ENE"
  i=i+1; MODS(INM_ME2BUT1ENE   )%NAME = "ME2BUT1ENE"
  i=i+1; MODS(INM_ME3BUT1ENE   )%NAME = "ME3BUT1ENE"
  i=i+1; MODS(INM_ME2BUT2ENE   )%NAME = "ME2BUT2ENE"
  i=i+1; MODS(INM_HEX1ENE      )%NAME = "HEX1ENE"
  i=i+1; MODS(INM_CHEX2ENE     )%NAME = "CHEX2ENE"
  i=i+1; MODS(INM_THEX2ENE     )%NAME = "THEX2ENE"
  i=i+1; MODS(INM_DM23BU2ENE   )%NAME = "DM23BU2ENE"
  i=i+1; MODS(INM_C2H2         )%NAME = "C2H2"
  i=i+1; MODS(INM_BENZENE      )%NAME = "BENZENE"
  i=i+1; MODS(INM_TOLUENE      )%NAME = "TOLUENE"
  i=i+1; MODS(INM_OXYL         )%NAME = "OXYL"
  i=i+1; MODS(INM_MXYL         )%NAME = "MXYL"
  i=i+1; MODS(INM_PXYL         )%NAME = "PXYL"
  i=i+1; MODS(INM_EBENZ        )%NAME = "EBENZ"
  i=i+1; MODS(INM_PBENZ        )%NAME = "PBENZ"
  i=i+1; MODS(INM_IPBENZ       )%NAME = "IPBENZ"
  i=i+1; MODS(INM_TM123B       )%NAME = "TM123B"
  i=i+1; MODS(INM_TM124B       )%NAME = "TM124B"
  i=i+1; MODS(INM_TM135B       )%NAME = "TM135B"
  i=i+1; MODS(INM_OETHTOL      )%NAME = "OETHTOL"
  i=i+1; MODS(INM_METHTOL      )%NAME = "METHTOL"
  i=i+1; MODS(INM_PETHTOL      )%NAME = "PETHTOL"
  i=i+1; MODS(INM_DIME35EB     )%NAME = "DIME35EB"
  i=i+1; MODS(INM_DIET35TOL    )%NAME = "DIET35TOL"
  i=i+1; MODS(INM_STYRENE      )%NAME = "STYRENE"
  i=i+1; MODS(INM_BENZAL       )%NAME = "BENZAL"
  i=i+1; MODS(INM_CH3CL        )%NAME = "CH3CL"
  i=i+1; MODS(INM_CH2CL2       )%NAME = "CH2CL2"
  i=i+1; MODS(INM_CHCL3        )%NAME = "CHCL3"
  i=i+1; MODS(INM_CH3CCL3      )%NAME = "CH3CCL3"
  i=i+1; MODS(INM_TCE          )%NAME = "TCE"
  i=i+1; MODS(INM_TRICLETH     )%NAME = "TRICLETH"
  i=i+1; MODS(INM_CDICLETH     )%NAME = "CDICLETH"
  i=i+1; MODS(INM_TDICLETH     )%NAME = "TDICLETH"
  i=i+1; MODS(INM_CH2CLCH2CL   )%NAME = "CH2CLCH2CL"
  i=i+1; MODS(INM_CCL2CH2      )%NAME = "CCL2CH2"
  i=i+1; MODS(INM_CL12PROP     )%NAME = "CL12PROP"
  i=i+1; MODS(INM_CHCL2CH3     )%NAME = "CHCL2CH3"
  i=i+1; MODS(INM_CH3CH2CL     )%NAME = "CH3CH2CL"
  i=i+1; MODS(INM_CHCL2CHCL2   )%NAME = "CHCL2CHCL2"
  i=i+1; MODS(INM_CH2CLCHCL2   )%NAME = "CH2CLCHCL2"
  i=i+1; MODS(INM_VINCL        )%NAME = "VINCL"
  i=i+1; MODS(INM_C4H6         )%NAME = "C4H6"
  i=i+1; MODS(INM_C5H8         )%NAME = "C5H8"
  i=i+1; MODS(INM_CH3OCHO      )%NAME = "CH3OCHO"
  i=i+1; MODS(INM_METHACET     )%NAME = "METHACET"
  i=i+1; MODS(INM_ETHACET      )%NAME = "ETHACET"
  i=i+1; MODS(INM_NPROACET     )%NAME = "NPROACET"
  i=i+1; MODS(INM_IPROACET     )%NAME = "IPROACET"
  i=i+1; MODS(INM_NBUTACET     )%NAME = "NBUTACET"
  i=i+1; MODS(INM_SBUTACET     )%NAME = "SBUTACET"
  i=i+1; MODS(INM_TBUACET      )%NAME = "TBUACET"
  i=i+1; MODS(INM_CH3OCH3      )%NAME = "CH3OCH3"
  i=i+1; MODS(INM_DIETETHER    )%NAME = "DIETETHER"
  i=i+1; MODS(INM_MTBE         )%NAME = "MTBE"
  i=i+1; MODS(INM_DIIPRETHER   )%NAME = "DIIPRETHER"
  i=i+1; MODS(INM_ETBE         )%NAME = "ETBE"
  i=i+1; MODS(INM_MO2EOL       )%NAME = "MO2EOL"
  i=i+1; MODS(INM_EOX2EOL      )%NAME = "EOX2EOL"
  i=i+1; MODS(INM_PR2OHMOX     )%NAME = "PR2OHMOX"
  i=i+1; MODS(INM_BUOX2ETOH    )%NAME = "BUOX2ETOH"
  i=i+1; MODS(INM_BOX2PROL     )%NAME = "BOX2PROL"
  i=i+1; MODS(INM_CH3BR        )%NAME = "CH3BR"
  i=i+1; MODS(INM_DIBRET       )%NAME = "DIBRET"
  i=i+1; MODS(INM_CH3COCH3     )%NAME = "CH3COCH3"
  i=i+1; MODS(INM_MEK          )%NAME = "MEK"
  i=i+1; MODS(INM_MPRK         )%NAME = "MPRK"
  i=i+1; MODS(INM_DIEK         )%NAME = "DIEK"
  i=i+1; MODS(INM_MIPK         )%NAME = "MIPK"
  i=i+1; MODS(INM_HEX2ONE      )%NAME = "HEX2ONE"
  i=i+1; MODS(INM_HEX3ONE      )%NAME = "HEX3ONE"
  i=i+1; MODS(INM_MIBK         )%NAME = "MIBK"
  i=i+1; MODS(INM_MTBK         )%NAME = "MTBK"
  i=i+1; MODS(INM_CYHEXONE     )%NAME = "CYHEXONE"
  i=i+1; MODS(INM_APINENE      )%NAME = "APINENE"
  i=i+1; MODS(INM_BPINENE      )%NAME = "BPINENE"
  i=i+1; MODS(INM_LIMONENE     )%NAME = "LIMONENE"
  i=i+1; MODS(INM_BCARY        )%NAME = "BCARY"
  i=i+1; MODS(INM_HCOOH        )%NAME = "HCOOH"
  i=i+1; MODS(INM_CH3CO2H      )%NAME = "CH3CO2H"
  i=i+1; MODS(INM_PROPACID     )%NAME = "PROPACID"
  i=i+1; MODS(INM_DMM          )%NAME = "DMM"
  i=i+1; MODS(INM_DMC          )%NAME = "DMC"
  i=i+1; MODS(INM_DMS          )%NAME = "DMS"
  i=i+1; MODS(INM_ETHOX        )%NAME = "ETHOX"
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
  INDRELAY(inm_temp ) = inf_temp  ;i=i+1
  INDRELAY(inm_pres ) = inf_pres  ;i=i+1
  INDRELAY(inm_RH   ) = inf_RH    ;i=i+1
  INDRELAY(inm_CS   ) = inf_CS    ;i=i+1
  INDRELAY(inm_swr  ) = inf_swr   ;i=i+1
  INDRELAY(inm_IPR  ) = inf_IPR   ;i=i+1

  INDRELAY(inm_H2SO4) = inf_H2SO4 ;i=i+1
  INDRELAY(inm_NH3  ) = inf_NH3   ;i=i+1
  INDRELAY(inm_DMA  ) = inf_DMA   ;i=i+1
  INDRELAY(inm_SO2  ) = inf_SO2   ;i=i+1
  INDRELAY(inm_NO   ) = inf_NO    ;i=i+1
  INDRELAY(inm_NO2  ) = inf_NO2   ;i=i+1
  INDRELAY(inm_CO   ) = inf_CO    ;i=i+1
  INDRELAY(inm_H2   ) = inf_H2    ;i=i+1
  INDRELAY(inm_O3   ) = inf_O3    ;i=i+1

  INDRELAY(inm_RES16) = inf_RES16 ;i=i+1
  INDRELAY(inm_RES17) = inf_RES17 ;i=i+1
  INDRELAY(inm_RES18) = inf_RES18 ;i=i+1
  INDRELAY(inm_RES19) = inf_RES19 ;i=i+1

  INDRELAY(INM_CH3OH       ) = INF_CH3OH      ; i=i+1
  INDRELAY(INM_C2H5OH      ) = INF_C2H5OH     ; i=i+1
  INDRELAY(INM_NPROPOL     ) = INF_NPROPOL    ; i=i+1
  INDRELAY(INM_IPROPOL     ) = INF_IPROPOL    ; i=i+1
  INDRELAY(INM_NBUTOL      ) = INF_NBUTOL     ; i=i+1
  INDRELAY(INM_BUT2OL      ) = INF_BUT2OL     ; i=i+1
  INDRELAY(INM_IBUTOL      ) = INF_IBUTOL     ; i=i+1
  INDRELAY(INM_TBUTOL      ) = INF_TBUTOL     ; i=i+1
  INDRELAY(INM_PECOH       ) = INF_PECOH      ; i=i+1
  INDRELAY(INM_IPEAOH      ) = INF_IPEAOH     ; i=i+1
  INDRELAY(INM_ME3BUOL     ) = INF_ME3BUOL    ; i=i+1
  INDRELAY(INM_IPECOH      ) = INF_IPECOH     ; i=i+1
  INDRELAY(INM_IPEBOH      ) = INF_IPEBOH     ; i=i+1
  INDRELAY(INM_CYHEXOL     ) = INF_CYHEXOL    ; i=i+1
  INDRELAY(INM_MIBKAOH     ) = INF_MIBKAOH    ; i=i+1
  INDRELAY(INM_ETHGLY      ) = INF_ETHGLY     ; i=i+1
  INDRELAY(INM_PROPGLY     ) = INF_PROPGLY    ; i=i+1
  INDRELAY(INM_MBO         ) = INF_MBO        ; i=i+1
  INDRELAY(INM_HCHO        ) = INF_HCHO       ; i=i+1
  INDRELAY(INM_CH3CHO      ) = INF_CH3CHO     ; i=i+1
  INDRELAY(INM_C2H5CHO     ) = INF_C2H5CHO    ; i=i+1
  INDRELAY(INM_C3H7CHO     ) = INF_C3H7CHO    ; i=i+1
  INDRELAY(INM_IPRCHO      ) = INF_IPRCHO     ; i=i+1
  INDRELAY(INM_C4H9CHO     ) = INF_C4H9CHO    ; i=i+1
  INDRELAY(INM_ACR         ) = INF_ACR        ; i=i+1
  INDRELAY(INM_MACR        ) = INF_MACR       ; i=i+1
  INDRELAY(INM_C4ALDB      ) = INF_C4ALDB     ; i=i+1
  INDRELAY(INM_CH4         ) = INF_CH4        ; i=i+1
  INDRELAY(INM_C2H6        ) = INF_C2H6       ; i=i+1
  INDRELAY(INM_C3H8        ) = INF_C3H8       ; i=i+1
  INDRELAY(INM_NC4H10      ) = INF_NC4H10     ; i=i+1
  INDRELAY(INM_IC4H10      ) = INF_IC4H10     ; i=i+1
  INDRELAY(INM_NC5H12      ) = INF_NC5H12     ; i=i+1
  INDRELAY(INM_IC5H12      ) = INF_IC5H12     ; i=i+1
  INDRELAY(INM_NEOP        ) = INF_NEOP       ; i=i+1
  INDRELAY(INM_NC6H14      ) = INF_NC6H14     ; i=i+1
  INDRELAY(INM_M2PE        ) = INF_M2PE       ; i=i+1
  INDRELAY(INM_M3PE        ) = INF_M3PE       ; i=i+1
  INDRELAY(INM_M22C4       ) = INF_M22C4      ; i=i+1
  INDRELAY(INM_M23C4       ) = INF_M23C4      ; i=i+1
  INDRELAY(INM_NC7H16      ) = INF_NC7H16     ; i=i+1
  INDRELAY(INM_M2HEX       ) = INF_M2HEX      ; i=i+1
  INDRELAY(INM_M3HEX       ) = INF_M3HEX      ; i=i+1
  INDRELAY(INM_NC8H18      ) = INF_NC8H18     ; i=i+1
  INDRELAY(INM_NC9H20      ) = INF_NC9H20     ; i=i+1
  INDRELAY(INM_NC10H22     ) = INF_NC10H22    ; i=i+1
  INDRELAY(INM_NC11H24     ) = INF_NC11H24    ; i=i+1
  INDRELAY(INM_NC12H26     ) = INF_NC12H26    ; i=i+1
  INDRELAY(INM_CHEX        ) = INF_CHEX       ; i=i+1
  INDRELAY(INM_C2H4        ) = INF_C2H4       ; i=i+1
  INDRELAY(INM_C3H6        ) = INF_C3H6       ; i=i+1
  INDRELAY(INM_BUT1ENE     ) = INF_BUT1ENE    ; i=i+1
  INDRELAY(INM_CBUT2ENE    ) = INF_CBUT2ENE   ; i=i+1
  INDRELAY(INM_TBUT2ENE    ) = INF_TBUT2ENE   ; i=i+1
  INDRELAY(INM_MEPROPENE   ) = INF_MEPROPENE  ; i=i+1
  INDRELAY(INM_PENT1ENE    ) = INF_PENT1ENE   ; i=i+1
  INDRELAY(INM_CPENT2ENE   ) = INF_CPENT2ENE  ; i=i+1
  INDRELAY(INM_TPENT2ENE   ) = INF_TPENT2ENE  ; i=i+1
  INDRELAY(INM_ME2BUT1ENE  ) = INF_ME2BUT1ENE ; i=i+1
  INDRELAY(INM_ME3BUT1ENE  ) = INF_ME3BUT1ENE ; i=i+1
  INDRELAY(INM_ME2BUT2ENE  ) = INF_ME2BUT2ENE ; i=i+1
  INDRELAY(INM_HEX1ENE     ) = INF_HEX1ENE    ; i=i+1
  INDRELAY(INM_CHEX2ENE    ) = INF_CHEX2ENE   ; i=i+1
  INDRELAY(INM_THEX2ENE    ) = INF_THEX2ENE   ; i=i+1
  INDRELAY(INM_DM23BU2ENE  ) = INF_DM23BU2ENE ; i=i+1
  INDRELAY(INM_C2H2        ) = INF_C2H2       ; i=i+1
  INDRELAY(INM_BENZENE     ) = INF_BENZENE    ; i=i+1
  INDRELAY(INM_TOLUENE     ) = INF_TOLUENE    ; i=i+1
  INDRELAY(INM_OXYL        ) = INF_OXYL       ; i=i+1
  INDRELAY(INM_MXYL        ) = INF_MXYL       ; i=i+1
  INDRELAY(INM_PXYL        ) = INF_PXYL       ; i=i+1
  INDRELAY(INM_EBENZ       ) = INF_EBENZ      ; i=i+1
  INDRELAY(INM_PBENZ       ) = INF_PBENZ      ; i=i+1
  INDRELAY(INM_IPBENZ      ) = INF_IPBENZ     ; i=i+1
  INDRELAY(INM_TM123B      ) = INF_TM123B     ; i=i+1
  INDRELAY(INM_TM124B      ) = INF_TM124B     ; i=i+1
  INDRELAY(INM_TM135B      ) = INF_TM135B     ; i=i+1
  INDRELAY(INM_OETHTOL     ) = INF_OETHTOL    ; i=i+1
  INDRELAY(INM_METHTOL     ) = INF_METHTOL    ; i=i+1
  INDRELAY(INM_PETHTOL     ) = INF_PETHTOL    ; i=i+1
  INDRELAY(INM_DIME35EB    ) = INF_DIME35EB   ; i=i+1
  INDRELAY(INM_DIET35TOL   ) = INF_DIET35TOL  ; i=i+1
  INDRELAY(INM_STYRENE     ) = INF_STYRENE    ; i=i+1
  INDRELAY(INM_BENZAL      ) = INF_BENZAL     ; i=i+1
  INDRELAY(INM_CH3CL       ) = INF_CH3CL      ; i=i+1
  INDRELAY(INM_CH2CL2      ) = INF_CH2CL2     ; i=i+1
  INDRELAY(INM_CHCL3       ) = INF_CHCL3      ; i=i+1
  INDRELAY(INM_CH3CCL3     ) = INF_CH3CCL3    ; i=i+1
  INDRELAY(INM_TCE         ) = INF_TCE        ; i=i+1
  INDRELAY(INM_TRICLETH    ) = INF_TRICLETH   ; i=i+1
  INDRELAY(INM_CDICLETH    ) = INF_CDICLETH   ; i=i+1
  INDRELAY(INM_TDICLETH    ) = INF_TDICLETH   ; i=i+1
  INDRELAY(INM_CH2CLCH2CL  ) = INF_CH2CLCH2CL ; i=i+1
  INDRELAY(INM_CCL2CH2     ) = INF_CCL2CH2    ; i=i+1
  INDRELAY(INM_CL12PROP    ) = INF_CL12PROP   ; i=i+1
  INDRELAY(INM_CHCL2CH3    ) = INF_CHCL2CH3   ; i=i+1
  INDRELAY(INM_CH3CH2CL    ) = INF_CH3CH2CL   ; i=i+1
  INDRELAY(INM_CHCL2CHCL2  ) = INF_CHCL2CHCL2 ; i=i+1
  INDRELAY(INM_CH2CLCHCL2  ) = INF_CH2CLCHCL2 ; i=i+1
  INDRELAY(INM_VINCL       ) = INF_VINCL      ; i=i+1
  INDRELAY(INM_C4H6        ) = INF_C4H6       ; i=i+1
  INDRELAY(INM_C5H8        ) = INF_C5H8       ; i=i+1
  INDRELAY(INM_CH3OCHO     ) = INF_CH3OCHO    ; i=i+1
  INDRELAY(INM_METHACET    ) = INF_METHACET   ; i=i+1
  INDRELAY(INM_ETHACET     ) = INF_ETHACET    ; i=i+1
  INDRELAY(INM_NPROACET    ) = INF_NPROACET   ; i=i+1
  INDRELAY(INM_IPROACET    ) = INF_IPROACET   ; i=i+1
  INDRELAY(INM_NBUTACET    ) = INF_NBUTACET   ; i=i+1
  INDRELAY(INM_SBUTACET    ) = INF_SBUTACET   ; i=i+1
  INDRELAY(INM_TBUACET     ) = INF_TBUACET    ; i=i+1
  INDRELAY(INM_CH3OCH3     ) = INF_CH3OCH3    ; i=i+1
  INDRELAY(INM_DIETETHER   ) = INF_DIETETHER  ; i=i+1
  INDRELAY(INM_MTBE        ) = INF_MTBE       ; i=i+1
  INDRELAY(INM_DIIPRETHER  ) = INF_DIIPRETHER ; i=i+1
  INDRELAY(INM_ETBE        ) = INF_ETBE       ; i=i+1
  INDRELAY(INM_MO2EOL      ) = INF_MO2EOL     ; i=i+1
  INDRELAY(INM_EOX2EOL     ) = INF_EOX2EOL    ; i=i+1
  INDRELAY(INM_PR2OHMOX    ) = INF_PR2OHMOX   ; i=i+1
  INDRELAY(INM_BUOX2ETOH   ) = INF_BUOX2ETOH  ; i=i+1
  INDRELAY(INM_BOX2PROL    ) = INF_BOX2PROL   ; i=i+1
  INDRELAY(INM_CH3BR       ) = INF_CH3BR      ; i=i+1
  INDRELAY(INM_DIBRET      ) = INF_DIBRET     ; i=i+1
  INDRELAY(INM_CH3COCH3    ) = INF_CH3COCH3   ; i=i+1
  INDRELAY(INM_MEK         ) = INF_MEK        ; i=i+1
  INDRELAY(INM_MPRK        ) = INF_MPRK       ; i=i+1
  INDRELAY(INM_DIEK        ) = INF_DIEK       ; i=i+1
  INDRELAY(INM_MIPK        ) = INF_MIPK       ; i=i+1
  INDRELAY(INM_HEX2ONE     ) = INF_HEX2ONE    ; i=i+1
  INDRELAY(INM_HEX3ONE     ) = INF_HEX3ONE    ; i=i+1
  INDRELAY(INM_MIBK        ) = INF_MIBK       ; i=i+1
  INDRELAY(INM_MTBK        ) = INF_MTBK       ; i=i+1
  INDRELAY(INM_CYHEXONE    ) = INF_CYHEXONE   ; i=i+1
  INDRELAY(INM_APINENE     ) = INF_APINENE    ; i=i+1
  INDRELAY(INM_BPINENE     ) = INF_BPINENE    ; i=i+1
  INDRELAY(INM_LIMONENE    ) = INF_LIMONENE   ; i=i+1
  INDRELAY(INM_BCARY       ) = INF_BCARY      ; i=i+1
  INDRELAY(INM_HCOOH       ) = INF_HCOOH      ; i=i+1
  INDRELAY(INM_CH3CO2H     ) = INF_CH3CO2H    ; i=i+1
  INDRELAY(INM_PROPACID    ) = INF_PROPACID   ; i=i+1
  INDRELAY(INM_DMM         ) = INF_DMM        ; i=i+1
  INDRELAY(INM_DMC         ) = INF_DMC        ; i=i+1
  INDRELAY(INM_DMS         ) = INF_DMS        ; i=i+1
  INDRELAY(INM_ETHOX       ) = INF_ETHOX      ; i=i+1

  ! As safety feature, check that here nothing is left out
  IF (i /= INM_LAST) THEN
    print FMT_FAT0, 'Trouble in FILL_INDRELAY_WITH_INDICES (input.f90). Number of input variables differs'
    print FMT_MSG, 'from what is assumed in INM_LAST. Fix and recompile.'
    print FMT_LEND
    print*,
    STOP
  END IF
end subroutine FILL_INDRELAY_WITH_INDICES

SUBROUTINE PUT_INPUT_IN_THEIR_PLACES(INPUT_ENV,INPUT_MCM,CONC_MAT)
  implicit none
  REAL(DP), intent(inout) :: CONC_MAT(:,:)
  REAL(DP), intent(in)    :: INPUT_ENV(:,:)
  REAL(DP), intent(in)    :: INPUT_MCM(:,:)
  integer::i
  i = 0

  IF (inf_temp  > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_temp ) = input_ENV(:,inf_temp ) ;i=i+1
  IF (inf_pres  > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_pres ) = input_ENV(:,inf_pres ) ;i=i+1
  IF (inf_RH    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_RH   ) = input_ENV(:,inf_RH   ) ;i=i+1
  IF (inf_CS    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_CS   ) = input_ENV(:,inf_CS   ) ;i=i+1
  IF (inf_swr   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_swr  ) = input_ENV(:,inf_swr  ) ;i=i+1
  IF (inf_IPR   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_IPR  ) = input_ENV(:,inf_IPR  ) ;i=i+1

  IF (inf_H2SO4 > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_H2SO4) = input_ENV(:,inf_H2SO4) ;i=i+1
  IF (inf_NH3   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_NH3  ) = input_ENV(:,inf_NH3  ) ;i=i+1
  IF (inf_DMA   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_DMA  ) = input_ENV(:,inf_DMA  ) ;i=i+1
  IF (inf_SO2   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_SO2  ) = input_ENV(:,inf_SO2  ) ;i=i+1
  IF (inf_NO    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_NO   ) = input_ENV(:,inf_NO   ) ;i=i+1
  IF (inf_NO2   > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_NO2  ) = input_ENV(:,inf_NO2  ) ;i=i+1
  IF (inf_CO    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_CO   ) = input_ENV(:,inf_CO   ) ;i=i+1
  IF (inf_H2    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_H2   ) = input_ENV(:,inf_H2   ) ;i=i+1
  IF (inf_O3    > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_O3   ) = input_ENV(:,inf_O3   ) ;i=i+1
  IF (inf_RES16 > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_RES16) = INPUT_MCM(:,inf_RES16) ;i=i+1
  IF (inf_RES17 > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_RES17) = INPUT_MCM(:,inf_RES17) ;i=i+1
  IF (inf_RES18 > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_RES18) = INPUT_MCM(:,inf_RES18) ;i=i+1
  IF (inf_RES19 > 0 .and. (ENV_file /= '')) CONC_MAT(:,inm_RES19) = INPUT_MCM(:,inf_RES19) ;i=i+1

  IF (INF_CH3OH        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3OH     ) = input_MCM(:,INF_CH3OH     ); i=i+1
  IF (INF_C2H5OH       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C2H5OH    ) = input_MCM(:,INF_C2H5OH    ); i=i+1
  IF (INF_NPROPOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NPROPOL   ) = input_MCM(:,INF_NPROPOL   ); i=i+1
  IF (INF_IPROPOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPROPOL   ) = input_MCM(:,INF_IPROPOL   ); i=i+1
  IF (INF_NBUTOL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NBUTOL    ) = input_MCM(:,INF_NBUTOL    ); i=i+1
  IF (INF_BUT2OL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BUT2OL    ) = input_MCM(:,INF_BUT2OL    ); i=i+1
  IF (INF_IBUTOL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IBUTOL    ) = input_MCM(:,INF_IBUTOL    ); i=i+1
  IF (INF_TBUTOL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TBUTOL    ) = input_MCM(:,INF_TBUTOL    ); i=i+1
  IF (INF_PECOH        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PECOH     ) = input_MCM(:,INF_PECOH     ); i=i+1
  IF (INF_IPEAOH       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPEAOH    ) = input_MCM(:,INF_IPEAOH    ); i=i+1
  IF (INF_ME3BUOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ME3BUOL   ) = input_MCM(:,INF_ME3BUOL   ); i=i+1
  IF (INF_IPECOH       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPECOH    ) = input_MCM(:,INF_IPECOH    ); i=i+1
  IF (INF_IPEBOH       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPEBOH    ) = input_MCM(:,INF_IPEBOH    ); i=i+1
  IF (INF_CYHEXOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CYHEXOL   ) = input_MCM(:,INF_CYHEXOL   ); i=i+1
  IF (INF_MIBKAOH      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MIBKAOH   ) = input_MCM(:,INF_MIBKAOH   ); i=i+1
  IF (INF_ETHGLY       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ETHGLY    ) = input_MCM(:,INF_ETHGLY    ); i=i+1
  IF (INF_PROPGLY      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PROPGLY   ) = input_MCM(:,INF_PROPGLY   ); i=i+1
  IF (INF_MBO          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MBO       ) = input_MCM(:,INF_MBO       ); i=i+1
  IF (INF_HCHO         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_HCHO      ) = input_MCM(:,INF_HCHO      ); i=i+1
  IF (INF_CH3CHO       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3CHO    ) = input_MCM(:,INF_CH3CHO    ); i=i+1
  IF (INF_C2H5CHO      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C2H5CHO   ) = input_MCM(:,INF_C2H5CHO   ); i=i+1
  IF (INF_C3H7CHO      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C3H7CHO   ) = input_MCM(:,INF_C3H7CHO   ); i=i+1
  IF (INF_IPRCHO       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPRCHO    ) = input_MCM(:,INF_IPRCHO    ); i=i+1
  IF (INF_C4H9CHO      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C4H9CHO   ) = input_MCM(:,INF_C4H9CHO   ); i=i+1
  IF (INF_ACR          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ACR       ) = input_MCM(:,INF_ACR       ); i=i+1
  IF (INF_MACR         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MACR      ) = input_MCM(:,INF_MACR      ); i=i+1
  IF (INF_C4ALDB       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C4ALDB    ) = input_MCM(:,INF_C4ALDB    ); i=i+1
  IF (INF_CH4          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH4       ) = input_MCM(:,INF_CH4       ); i=i+1
  IF (INF_C2H6         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C2H6      ) = input_MCM(:,INF_C2H6      ); i=i+1
  IF (INF_C3H8         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C3H8      ) = input_MCM(:,INF_C3H8      ); i=i+1
  IF (INF_NC4H10       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC4H10    ) = input_MCM(:,INF_NC4H10    ); i=i+1
  IF (INF_IC4H10       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IC4H10    ) = input_MCM(:,INF_IC4H10    ); i=i+1
  IF (INF_NC5H12       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC5H12    ) = input_MCM(:,INF_NC5H12    ); i=i+1
  IF (INF_IC5H12       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IC5H12    ) = input_MCM(:,INF_IC5H12    ); i=i+1
  IF (INF_NEOP         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NEOP      ) = input_MCM(:,INF_NEOP      ); i=i+1
  IF (INF_NC6H14       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC6H14    ) = input_MCM(:,INF_NC6H14    ); i=i+1
  IF (INF_M2PE         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M2PE      ) = input_MCM(:,INF_M2PE      ); i=i+1
  IF (INF_M3PE         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M3PE      ) = input_MCM(:,INF_M3PE      ); i=i+1
  IF (INF_M22C4        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M22C4     ) = input_MCM(:,INF_M22C4     ); i=i+1
  IF (INF_M23C4        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M23C4     ) = input_MCM(:,INF_M23C4     ); i=i+1
  IF (INF_NC7H16       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC7H16    ) = input_MCM(:,INF_NC7H16    ); i=i+1
  IF (INF_M2HEX        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M2HEX     ) = input_MCM(:,INF_M2HEX     ); i=i+1
  IF (INF_M3HEX        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_M3HEX     ) = input_MCM(:,INF_M3HEX     ); i=i+1
  IF (INF_NC8H18       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC8H18    ) = input_MCM(:,INF_NC8H18    ); i=i+1
  IF (INF_NC9H20       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC9H20    ) = input_MCM(:,INF_NC9H20    ); i=i+1
  IF (INF_NC10H22      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC10H22   ) = input_MCM(:,INF_NC10H22   ); i=i+1
  IF (INF_NC11H24      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC11H24   ) = input_MCM(:,INF_NC11H24   ); i=i+1
  IF (INF_NC12H26      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NC12H26   ) = input_MCM(:,INF_NC12H26   ); i=i+1
  IF (INF_CHEX         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CHEX      ) = input_MCM(:,INF_CHEX      ); i=i+1
  IF (INF_C2H4         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C2H4      ) = input_MCM(:,INF_C2H4      ); i=i+1
  IF (INF_C3H6         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C3H6      ) = input_MCM(:,INF_C3H6      ); i=i+1
  IF (INF_BUT1ENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BUT1ENE   ) = input_MCM(:,INF_BUT1ENE   ); i=i+1
  IF (INF_CBUT2ENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CBUT2ENE  ) = input_MCM(:,INF_CBUT2ENE  ); i=i+1
  IF (INF_TBUT2ENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TBUT2ENE  ) = input_MCM(:,INF_TBUT2ENE  ); i=i+1
  IF (INF_MEPROPENE    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MEPROPENE ) = input_MCM(:,INF_MEPROPENE ); i=i+1
  IF (INF_PENT1ENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PENT1ENE  ) = input_MCM(:,INF_PENT1ENE  ); i=i+1
  IF (INF_CPENT2ENE    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CPENT2ENE ) = input_MCM(:,INF_CPENT2ENE ); i=i+1
  IF (INF_TPENT2ENE    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TPENT2ENE ) = input_MCM(:,INF_TPENT2ENE ); i=i+1
  IF (INF_ME2BUT1ENE   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ME2BUT1ENE) = input_MCM(:,INF_ME2BUT1ENE); i=i+1
  IF (INF_ME3BUT1ENE   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ME3BUT1ENE) = input_MCM(:,INF_ME3BUT1ENE); i=i+1
  IF (INF_ME2BUT2ENE   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ME2BUT2ENE) = input_MCM(:,INF_ME2BUT2ENE); i=i+1
  IF (INF_HEX1ENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_HEX1ENE   ) = input_MCM(:,INF_HEX1ENE   ); i=i+1
  IF (INF_CHEX2ENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CHEX2ENE  ) = input_MCM(:,INF_CHEX2ENE  ); i=i+1
  IF (INF_THEX2ENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_THEX2ENE  ) = input_MCM(:,INF_THEX2ENE  ); i=i+1
  IF (INF_DM23BU2ENE   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DM23BU2ENE) = input_MCM(:,INF_DM23BU2ENE); i=i+1
  IF (INF_C2H2         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C2H2      ) = input_MCM(:,INF_C2H2      ); i=i+1
  IF (INF_BENZENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BENZENE   ) = input_MCM(:,INF_BENZENE   ); i=i+1
  IF (INF_TOLUENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TOLUENE   ) = input_MCM(:,INF_TOLUENE   ); i=i+1
  IF (INF_OXYL         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_OXYL      ) = input_MCM(:,INF_OXYL      ); i=i+1
  IF (INF_MXYL         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MXYL      ) = input_MCM(:,INF_MXYL      ); i=i+1
  IF (INF_PXYL         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PXYL      ) = input_MCM(:,INF_PXYL      ); i=i+1
  IF (INF_EBENZ        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_EBENZ     ) = input_MCM(:,INF_EBENZ     ); i=i+1
  IF (INF_PBENZ        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PBENZ     ) = input_MCM(:,INF_PBENZ     ); i=i+1
  IF (INF_IPBENZ       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPBENZ    ) = input_MCM(:,INF_IPBENZ    ); i=i+1
  IF (INF_TM123B       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TM123B    ) = input_MCM(:,INF_TM123B    ); i=i+1
  IF (INF_TM124B       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TM124B    ) = input_MCM(:,INF_TM124B    ); i=i+1
  IF (INF_TM135B       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TM135B    ) = input_MCM(:,INF_TM135B    ); i=i+1
  IF (INF_OETHTOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_OETHTOL   ) = input_MCM(:,INF_OETHTOL   ); i=i+1
  IF (INF_METHTOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_METHTOL   ) = input_MCM(:,INF_METHTOL   ); i=i+1
  IF (INF_PETHTOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PETHTOL   ) = input_MCM(:,INF_PETHTOL   ); i=i+1
  IF (INF_DIME35EB     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIME35EB  ) = input_MCM(:,INF_DIME35EB  ); i=i+1
  IF (INF_DIET35TOL    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIET35TOL ) = input_MCM(:,INF_DIET35TOL ); i=i+1
  IF (INF_STYRENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_STYRENE   ) = input_MCM(:,INF_STYRENE   ); i=i+1
  IF (INF_BENZAL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BENZAL    ) = input_MCM(:,INF_BENZAL    ); i=i+1
  IF (INF_CH3CL        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3CL     ) = input_MCM(:,INF_CH3CL     ); i=i+1
  IF (INF_CH2CL2       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH2CL2    ) = input_MCM(:,INF_CH2CL2    ); i=i+1
  IF (INF_CHCL3        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CHCL3     ) = input_MCM(:,INF_CHCL3     ); i=i+1
  IF (INF_CH3CCL3      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3CCL3   ) = input_MCM(:,INF_CH3CCL3   ); i=i+1
  IF (INF_TCE          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TCE       ) = input_MCM(:,INF_TCE       ); i=i+1
  IF (INF_TRICLETH     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TRICLETH  ) = input_MCM(:,INF_TRICLETH  ); i=i+1
  IF (INF_CDICLETH     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CDICLETH  ) = input_MCM(:,INF_CDICLETH  ); i=i+1
  IF (INF_TDICLETH     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TDICLETH  ) = input_MCM(:,INF_TDICLETH  ); i=i+1
  IF (INF_CH2CLCH2CL   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH2CLCH2CL) = input_MCM(:,INF_CH2CLCH2CL); i=i+1
  IF (INF_CCL2CH2      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CCL2CH2   ) = input_MCM(:,INF_CCL2CH2   ); i=i+1
  IF (INF_CL12PROP     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CL12PROP  ) = input_MCM(:,INF_CL12PROP  ); i=i+1
  IF (INF_CHCL2CH3     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CHCL2CH3  ) = input_MCM(:,INF_CHCL2CH3  ); i=i+1
  IF (INF_CH3CH2CL     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3CH2CL  ) = input_MCM(:,INF_CH3CH2CL  ); i=i+1
  IF (INF_CHCL2CHCL2   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CHCL2CHCL2) = input_MCM(:,INF_CHCL2CHCL2); i=i+1
  IF (INF_CH2CLCHCL2   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH2CLCHCL2) = input_MCM(:,INF_CH2CLCHCL2); i=i+1
  IF (INF_VINCL        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_VINCL     ) = input_MCM(:,INF_VINCL     ); i=i+1
  IF (INF_C4H6         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C4H6      ) = input_MCM(:,INF_C4H6      ); i=i+1
  IF (INF_C5H8         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_C5H8      ) = input_MCM(:,INF_C5H8      ); i=i+1
  IF (INF_CH3OCHO      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3OCHO   ) = input_MCM(:,INF_CH3OCHO   ); i=i+1
  IF (INF_METHACET     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_METHACET  ) = input_MCM(:,INF_METHACET  ); i=i+1
  IF (INF_ETHACET      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ETHACET   ) = input_MCM(:,INF_ETHACET   ); i=i+1
  IF (INF_NPROACET     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NPROACET  ) = input_MCM(:,INF_NPROACET  ); i=i+1
  IF (INF_IPROACET     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_IPROACET  ) = input_MCM(:,INF_IPROACET  ); i=i+1
  IF (INF_NBUTACET     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_NBUTACET  ) = input_MCM(:,INF_NBUTACET  ); i=i+1
  IF (INF_SBUTACET     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_SBUTACET  ) = input_MCM(:,INF_SBUTACET  ); i=i+1
  IF (INF_TBUACET      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_TBUACET   ) = input_MCM(:,INF_TBUACET   ); i=i+1
  IF (INF_CH3OCH3      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3OCH3   ) = input_MCM(:,INF_CH3OCH3   ); i=i+1
  IF (INF_DIETETHER    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIETETHER ) = input_MCM(:,INF_DIETETHER ); i=i+1
  IF (INF_MTBE         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MTBE      ) = input_MCM(:,INF_MTBE      ); i=i+1
  IF (INF_DIIPRETHER   > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIIPRETHER) = input_MCM(:,INF_DIIPRETHER); i=i+1
  IF (INF_ETBE         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ETBE      ) = input_MCM(:,INF_ETBE      ); i=i+1
  IF (INF_MO2EOL       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MO2EOL    ) = input_MCM(:,INF_MO2EOL    ); i=i+1
  IF (INF_EOX2EOL      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_EOX2EOL   ) = input_MCM(:,INF_EOX2EOL   ); i=i+1
  IF (INF_PR2OHMOX     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PR2OHMOX  ) = input_MCM(:,INF_PR2OHMOX  ); i=i+1
  IF (INF_BUOX2ETOH    > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BUOX2ETOH ) = input_MCM(:,INF_BUOX2ETOH ); i=i+1
  IF (INF_BOX2PROL     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BOX2PROL  ) = input_MCM(:,INF_BOX2PROL  ); i=i+1
  IF (INF_CH3BR        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3BR     ) = input_MCM(:,INF_CH3BR     ); i=i+1
  IF (INF_DIBRET       > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIBRET    ) = input_MCM(:,INF_DIBRET    ); i=i+1
  IF (INF_CH3COCH3     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3COCH3  ) = input_MCM(:,INF_CH3COCH3  ); i=i+1
  IF (INF_MEK          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MEK       ) = input_MCM(:,INF_MEK       ); i=i+1
  IF (INF_MPRK         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MPRK      ) = input_MCM(:,INF_MPRK      ); i=i+1
  IF (INF_DIEK         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DIEK      ) = input_MCM(:,INF_DIEK      ); i=i+1
  IF (INF_MIPK         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MIPK      ) = input_MCM(:,INF_MIPK      ); i=i+1
  IF (INF_HEX2ONE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_HEX2ONE   ) = input_MCM(:,INF_HEX2ONE   ); i=i+1
  IF (INF_HEX3ONE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_HEX3ONE   ) = input_MCM(:,INF_HEX3ONE   ); i=i+1
  IF (INF_MIBK         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MIBK      ) = input_MCM(:,INF_MIBK      ); i=i+1
  IF (INF_MTBK         > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_MTBK      ) = input_MCM(:,INF_MTBK      ); i=i+1
  IF (INF_CYHEXONE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CYHEXONE  ) = input_MCM(:,INF_CYHEXONE  ); i=i+1
  IF (INF_APINENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_APINENE   ) = input_MCM(:,INF_APINENE   ); i=i+1
  IF (INF_BPINENE      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BPINENE   ) = input_MCM(:,INF_BPINENE   ); i=i+1
  IF (INF_LIMONENE     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_LIMONENE  ) = input_MCM(:,INF_LIMONENE  ); i=i+1
  IF (INF_BCARY        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_BCARY     ) = input_MCM(:,INF_BCARY     ); i=i+1
  IF (INF_HCOOH        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_HCOOH     ) = input_MCM(:,INF_HCOOH     ); i=i+1
  IF (INF_CH3CO2H      > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_CH3CO2H   ) = input_MCM(:,INF_CH3CO2H   ); i=i+1
  IF (INF_PROPACID     > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_PROPACID  ) = input_MCM(:,INF_PROPACID  ); i=i+1
  IF (INF_DMM          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DMM       ) = input_MCM(:,INF_DMM       ); i=i+1
  IF (INF_DMC          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DMC       ) = input_MCM(:,INF_DMC       ); i=i+1
  IF (INF_DMS          > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_DMS       ) = input_MCM(:,INF_DMS       ); i=i+1
  IF (INF_ETHOX        > 0 .and. (MCM_file /= '')) CONC_MAT(:,INM_ETHOX     ) = input_MCM(:,INF_ETHOX     ); i=i+1

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
