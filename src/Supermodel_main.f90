PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)

USE second_Precision,  ONLY : dp ! KPP Numerical type
USE constants
USE AUXILLARIES
#ifdef ISACDC
  USE ACDC_NH3
  USE ACDC_DMA
#endif
use OUTPUT

IMPLICIT NONE

type(input_mod) :: conc_MODS(5)
integer   :: ind_Temp  = 1
integer   :: ind_SA   = 2
integer   :: ind_NH3 = 3
integer   :: ind_DMA  = 4
integer   :: ind_CS   = 5

!Variable declaration time
TYPE(type_ambient) :: Gases
REAL(dp), allocatable :: conc(:,:)
REAL(dp) :: c_acid = 1e7*1d6
REAL(dp) :: c_base = 1d8*1d6
REAL(dp) :: c_dma
REAL(dp) :: c_org = 0*1d12*1d6
REAL(dp) :: cs_H2SO4 = 0.0003d0
REAL(dp) :: TempK = K0
REAL(dp) :: IPR = 3d6
REAL(dp) :: J_ACDC_NH3
REAL(dp) :: J_ACDC_DMA
REAL(dp) :: J_NH3_BY_IONS(3)
REAL(dp) :: acdc_cluster_diam = 2.17d-9
LOGICAL :: ACDC_solve_ss = .true. , NUCLEATION = .true., ACDC = .true.
INTEGER :: rows, cols, i,j
CHARACTER(222) :: CASE, PATH,JDch, Description
type(input_mod) :: M_NH3 = input_mod(min=1.000000d4,max=1.000000d10,sig=2.348050d0,mju=12.000000d0,fv=0.966667d0,ph=-2.280000d0,am=1.000000d0, MODE=2)
print FMT_LEND,

!Create/Open outup file netcdf
!Variable initialization
!Species properties in gas/particle phase
!Gas phase
!Aerosol phase(distribution & composition)
!Boundary conditions(dilution, losses, light,...)

CALL SET_MODIFIERS()
MODELTIME%JD = 1
Case = 'NJA'
Description = 'Nanjing day 1 Nominal values'
write(JDch,'(i0.3)') MODELTIME%JD
write(path, '(a,i0,a)') 'input/NANJING/case_Nanjing',MODELTIME%JD,'.txt'
print FMT_SUB, 'Input file: '//Trim(path)

open(522, file = TRIM(path))
CALL OPEN_GASFILE(('output/'//TRIM(JDch)//TRIM(CASE)//'.nc'), CONC_MODS, Description)!, model_options)

rows = ROWCOUNT(522)
cols = COLCOUNT(522)
ALLOCATE(conc(rows, cols))
do i=1,rows
  read(522, *) (conc(i,j), j=1,cols)
end do
conc(:,3) = conc(:,3) + K0
conc(:,5) = conc(:,5)*1d6
conc(:,6) = conc(:,6)*1d6
conc(:,7) = conc(:,7)*1d6
CONC_MODS(ind_temp)%shift = 0d0
CONC_MODS(ind_NH3)%multi = 2d0
CONC_MODS(ind_NH3)%MODE = 1
CONC_MODS(ind_NH3)%sig = 4.5
!CONC_MODS(ind_DMA)%am = 100d0
!CONC_MODS(ind_CS)%am = 5d-1

! Before the main loop starts, tell the user if any modifiers differ from default values
CALL CHECK_MODIFIERS()

!Main loop time: Eulerian forward integration
DO WHILE (MODELTIME%sec < MODELTIME%SIM_TIME_S)

  if (MODELTIME%printnow) print FMT_TIME, MODELTIME%hms
  conc_MODS(ind_NH3) = input_mod(min=1.000000d4,max=1.000000d10,sig=2.348050d0,mju=12.000000d0,fv=0.966667d0,ph=-2.280000d0,am=1.000000d0, MODE=2)

  ! TempK    = PERIODICAL(conc_MODS(ind_temp))
  ! c_base   = NORMALD(conc_MODS(ind_NH3))
  TempK    = interp(conc(:,1), conc(:,3))  .mod. conc_MODS(ind_Temp)
  CS_H2SO4 = interp(conc(:,1), conc(:,4))  .mod. conc_MODS(ind_CS)
  c_acid   = interp(conc(:,1), conc(:,5))  .mod. conc_MODS(ind_SA)
  c_base   = interp(conc(:,1), conc(:,6))  .mod. conc_MODS(ind_NH3)
  c_dma    = interp(conc(:,1), conc(:,7))  .mod. conc_MODS(ind_dma)

  if (NUCLEATION) THEN

#ifdef ISACDC
  if (ACDC) CALL ACDC_J() ! SUBROUTINE in CONTAINS of this file
#endif
#ifndef ISACDC
  if (MODELTIME%printnow) print*, NORMALD(MODELTIME%sec,Jout_par)
#endif

  END if ! NUCLEATION

  ! Condensation
  ! Coagulation
  ! Deposition


  ! Write printouts to screen and outputs to netcdf-file
  if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION()
  if (MODELTIME%savenow) CALL SAVE_GASES(TempK, C_acid, C_base, C_DMA, J_ACDC_NH3, J_ACDC_DMA, CS_H2SO4, Gases, CONC_MODS)

  MODELTIME = ADD(MODELTIME)
  if (MODELTIME%printnow) print *
END DO	! Main loop time: Eulerian forward integration

print FMT_TIME, MODELTIME%hms
CALL PRINT_KEY_INFORMATION()
call SAVE_GASES(TempK, C_acid, C_base, C_DMA, J_ACDC_NH3, J_ACDC_DMA, CS_H2SO4, Gases, CONC_MODS)
!Close output file netcdf
CALL CLOSE_FILES()

CONTAINS

SUBROUTINE ACDC_J()
! =================================================================================================
! ACDC Nucleation. See that all values here are in SI-UNITS: CUBIC METERS, KELVINS AND PASCALS.
! Written by Tinja Olenius
! .................................................................................................
! Input for get_acdc_J:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_base:            base (ammonia) concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric amincid [1/s]
! TempK:             Temperature in Kelvins
! IPR:               Ion production rate in ion pairs per second [1/m3/s]. 3d6 is a good guestimate
! dt:                Main time step [s]
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_NH3:        Particle formation rate due to ammonia [1/s/m3]. Sum of J_NH3_BY_IONS
! acdc_cluster_diam: Outgrowing cluster diameter [m]
! J_NH3_BY_IONS:     Particle formation rate by ions (neutral, negative and positive) [1/s/m3].
!                    The outgrowing cluster typically has 5 to 6 H2SO4 in it
! .................................................................................................
! Input for get_acdc_D:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_dma:             DMA concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
! TempK:             Temperature in Kelvins
! dt:                Main time step [s]
! time:              model time; used only for outputting information in cluster
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_DMA:        Particle formation rate due to DMA [1/s/m3]
! acdc_cluster_diam: Outgrowing cluster diameter [m]. The cluster typically has 5 to 6 H2SO4 in it
! =================================================================================================

! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!
  ! Speed up program by ignoring nucleation when there is none
  if (c_base > 1d12 .or. J_ACDC_NH3 > 1d-6) THEN

    ! from cc to m3:
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,IPR,MODELTIME,&
          ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (MODELTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
  END IF

  ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  ! if (c_dma > 1d6 .or. J_ACDC_DMA > 1d-6) THEN
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,MODELTIME,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
  ! ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
  !   if (EVENMIN(time, 180)) print FMT_SUB, 'DMA IGNORED'
  ! END IF

  ! ................................................................................................
  ! END ACDC Nucleation
  ! ================================================================================================
END SUBROUTINE ACDC_J

SUBROUTINE SET_MODIFIERS()

  conc_MODS(ind_Temp) %Name  = "Temperature     "
  conc_MODS(ind_SA)   %Name  = "H2SO4           "
  conc_MODS(ind_NH3)  %Name  = "Base_NH3        "
  conc_MODS(ind_DMA)  %Name  = "DMA             "
  conc_MODS(ind_CS)   %Name  = "C_sink          "

END SUBROUTINE SET_MODIFIERS

SUBROUTINE SOME_OTHER_NUCLEATION_TYPE(J1, J2)
  IMPLICIT NONE
  REAL(dp) :: J1, J2
  J1 = 1d0
  J2 = 1d0
END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE

SUBROUTINE PRINT_KEY_INFORMATION()
  IMPLICIT NONE
  print FMT10_2CVU,'ACID C: ', c_acid*1d-6, ' [1/cm3]', 'Temp:', TempK, 'Kelvin'
  print FMT10_2CVU, 'NH3 C:', c_base*1d-6, ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
  print FMT10_2CVU, 'DMA C:', c_dma*1d-6 , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
  print FMT10_3CVU, 'Jion1:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion1:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion1:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
  print FMT10_CVU, 'C-sink:', CS_H2SO4 , ' [1/s]'
  print FMT_LEND,
END SUBROUTINE PRINT_KEY_INFORMATION

SUBROUTINE CHECK_MODIFIERS()
  IMPLICIT NONE
  type(input_mod) :: test
  integer                 :: i
  do i=1,size(CONC_MODS)
  IF (conc_MODS(i)%MODE > 0) THEN
    print FMT_WARN0, 'Replacing input for '//TRIM(conc_MODS(i)%name)//' with parametrized function.'
  ELSE
    if (ABS(conc_MODS(i)%shift - test%shift) > 1d-9) THEN
      print FMT_WARN1, 'Adding a constant to '//TRIM(conc_MODS(i)%name)//', value is: ',conc_MODS(i)%shift
    ELSEIF (ABS(conc_MODS(i)%multi - test%multi) > 1d-9) THEN
      print FMT_WARN1, 'Multiplying '//TRIM(conc_MODS(i)%name)//' with: ',conc_MODS(i)%multi
    END IF
  END IF
  END DO
END SUBROUTINE CHECK_MODIFIERS

END PROGRAM SUPERMODEL
