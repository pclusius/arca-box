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

CHARACTER(16)  :: modifier_names(5)
real(dp)  :: cons_multipliers(5)
real(dp)  :: cons_shifters(5)

integer   :: ind_Temp_shift  = 1
integer   :: ind_SA_factor   = 2
integer   :: ind_base_factor = 3
integer   :: ind_DMA_factor  = 4
integer   :: ind_CS_factor   = 5

!Variable declaration
TYPE(timetype) :: time
TYPE(type_ambient) :: Gases
REAL(dp), allocatable :: conc(:,:)
REAL(dp) :: c_acid =1e7*1d6
REAL(dp) :: c_base =1d8*1d6
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
type(parametered_input) :: MOD_NH3   = parametered_input(sigma=2d0,min_c=10d0,max_c=200d0,peaktime=8.30,omega=0d0,amplitude=-1.1d0, LOGSCALE=.true.)
type(parametered_input) :: T_par  = parametered_input(sigma=2d0,min_c=25d0+K0,max_c=35d0+K0,peaktime=9,omega=0d0,amplitude=-1.1d0, LOGSCALE=.false.)
! type(parametered_input) :: MOD_H2SO4 = parametered_input(sigma=0.5d0,min_c=4.8d05,max_c=5.48d07,peaktime=8d0,omega=-1.1d0,amplitude=1.1d0, LOGSCALE=.true.)
! type(parametered_input) :: Jout_par  = parametered_input(sigma=2d0,min_c=1.d-4,max_c=5d0,peaktime=9,omega=0d0,amplitude=-1.1d0, LOGSCALE=.false.)
INTEGER :: rows, cols, i,j

!Create/Open outup file netcdf
!Variable initialization
!Species properties in gas/particle phase
!Gas phase
!Aerosol phase(distribution & composition)
!Boundary conditions(dilution, losses, light,...)
! open(9333,FILE='Nanjing2_CS05_DMA10_K-70.dat')
CALL SET_MODIFIERS()

open(522, file = 'input/NANJING/case_Nanjing2.txt')

open(9333,FILE='output/xxx.dat')
write(9333,'(8(a20))') 'time(s)', 'c_acid(1/m3)', 'c_base(1/m3)', 'J_ACDC_NH3(1/m3/s)', 'c_dma(1/m3)', 'J_ACDC_DMA(1/m3/s)', 'TempK(K)', 'CS(1/s)'

CALL OPEN_GASFILE('output/OutputGas.nc', modifier_names, cons_multipliers, cons_shifters)!, model_options)

rows = ROWCOUNT(522)
cols = COLCOUNT(522)
ALLOCATE(conc(rows, cols))
do i=1,rows
  read(522, *) (conc(i,j), j=1,cols)
end do

print FMT_Tm, time%hms

!Main loop time: Eulerian forward integration
DO WHILE (time%sec < time%SIM_TIME_S)
  if (time%printnow) THEN
    print*,
    print FMT_Tm, time%hms
  end if
  !Photolysis
  !Chemistry

  TempK    = interp(time%sec, conc(:,1), conc(:,3))  * cons_multipliers(ind_Temp_shift)  + cons_shifters(ind_Temp_shift)
  CS_H2SO4 = interp(time%sec, conc(:,1), conc(:,4))  * cons_multipliers(ind_CS_factor)   + cons_shifters(ind_CS_factor)
  c_acid   = interp(time%sec, conc(:,1), conc(:,5))  * cons_multipliers(ind_SA_factor)   + cons_shifters(ind_SA_factor)
  c_base   = interp(time%sec, conc(:,1), conc(:,6))  * cons_multipliers(ind_base_factor) + cons_shifters(ind_base_factor)
  c_dma    = interp(time%sec, conc(:,1), conc(:,7))  * cons_multipliers(ind_DMA_factor)  + cons_shifters(ind_DMA_factor)

  if (NUCLEATION) THEN

#ifdef ISACDC
  if (ACDC) CALL ACDC_J() ! SUBROUTINE in CONTAINS of this file
#endif

#ifndef ISACDC
  if (time%printnow) print*, NORMALD(time%sec,Jout_par)
#endif

  if (time%printnow) THEN
    print FMT10_2CVU,'ACID C: ', c_acid*1d-6, ' [1/cm3]', 'Temp:', TempK, 'Kelvin'
    print FMT10_2CVU, 'NH3 C:', c_base*1d-6, ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
    print FMT10_2CVU, 'DMA C:', c_dma*1d-6 , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
    print FMT10_3CVU, 'Jion1:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion1:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion1:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    print FMT10_CVU, 'C-sink:', CS_H2SO4 , ' [1/s]'
    print FMT_LEND,
  END IF

  END if ! NUCLEATION

  !Condensation
  !Coagulation
  !Deposition

  !Write output to file
  if (time%savenow) THEN
    write(9333,'(8(es20.8))') time%sec, c_acid, c_base, J_ACDC_NH3, c_dma, J_ACDC_DMA, TempK, CS_H2SO4
    call SAVE_GASES(time, TempK, C_acid, C_base, C_DMA, J_ACDC_NH3, J_ACDC_DMA, CS_H2SO4, Gases, cons_multipliers, cons_shifters)
    ! call output_ambient(time, ambient)
    !call output_particles(time, particles)


  END IF

  time = time + time%dt

END DO	! Main loop time: Eulerian forward integration

print FMT_Tm, time%hms

write(9333,'(8(es20.8))') time%sec, c_acid, c_base, J_ACDC_NH3, c_dma, J_ACDC_DMA, TempK, CS_H2SO4
call SAVE_GASES(time, TempK, C_acid, C_base, C_DMA, J_ACDC_NH3, J_ACDC_DMA, CS_H2SO4, Gases, cons_multipliers, cons_shifters)
! call output_time(time)
! call output_ambient(time, ambient)


  !Close output file netcdf
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
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
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
  ! if (c_base > 1d12 .or. J_ACDC_NH3 > 1d-6) THEN
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,IPR,time,&
          ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
  ! ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
  !   if (EVENMIN(time, 15)) print FMT_SUB, 'NH3 IGNORED'
  ! END IF

  ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  ! if (c_dma > 1d6 .or. J_ACDC_DMA > 1d-6) THEN
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,time,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
  ! ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
  !   if (EVENMIN(time, 180)) print FMT_SUB, 'DMA IGNORED'
  ! END IF

  ! ................................................................................................
  ! END ACDC Nucleation
  ! ================================================================================================
END SUBROUTINE ACDC_J

SUBROUTINE SET_MODIFIERS()
  modifier_names(ind_Temp_shift)    = "Temperature     "
  modifier_names(ind_SA_factor)     = "H2SO4           "
  modifier_names(ind_base_factor)   = "Base_NH3        "
  modifier_names(ind_DMA_factor)    = "DMA             "
  modifier_names(ind_CS_factor)     = "C_sink          "

  cons_multipliers(ind_Temp_shift)  = 1d0
  cons_multipliers(ind_SA_factor)   = 1d0
  cons_multipliers(ind_base_factor) = 1d0
  cons_multipliers(ind_DMA_factor)  = 1d0
  cons_multipliers(ind_CS_factor)   = 1d0

  cons_shifters(ind_Temp_shift)     = 0d0 + K0
  cons_shifters(ind_SA_factor)      = 0d0
  cons_shifters(ind_base_factor)    = 0d0
  cons_shifters(ind_DMA_factor)     = 0d0
  cons_shifters(ind_CS_factor)      = 0d0
END SUBROUTINE SET_MODIFIERS

SUBROUTINE SOME_OTHER_NUCLEATION_TYPE(J1, J2)
  REAL(dp) :: J1, J2
  J1 = 1d0
  J2 = 1d0
END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE

END PROGRAM SUPERMODEL
