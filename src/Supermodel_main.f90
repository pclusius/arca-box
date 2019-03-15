PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)

!USE module statements
USE second_Precision,  ONLY : dp ! KPP Numerical type
USE constants
USE AUXILLARIES
#ifdef ISACDC
  USE ACDC_NH3
  USE ACDC_DMA
#endif
USE ACDC_DMA

IMPLICIT NONE


  !Variable declaration
  REAL(dp) ::			&
    dt,				&	!integration time step [s]
    time,     & !simulation time [s]
    sim_time

    REAL(dp) :: c_acid =1e7*1d6
    REAL(dp) :: c_base =1d8*1d6
    REAL(dp) :: c_dma
    REAL(dp) :: c_org = 1d12*1d6
    REAL(dp) :: cs_H2SO4 = 0.000d6
    REAL(dp) :: TempK = K0 + 15d0
    REAL(dp) :: IPR = 3d6
    REAL(dp) :: J_ACDC_NH3
    REAL(dp) :: J_ACDC_DMA
    REAL(dp) :: J_NH3_BY_IONS(3)
    REAL(dp) :: acdc_cluster_diam = 1.17d-9
    LOGICAL :: ACDC_solve_ss = .false., NUCLEATION=.true., ACDC=.true.
    type(parametered_input) :: MOD_H2SO4 = parametered_input(sigma=3d0,min_c=5d7,max_c=5d7,peaktime=12,omega=1d0,amplitude=-1.1d0, LOGSCALE=.true.)
    type(parametered_input) :: MOD_NH3   = parametered_input(sigma=5d0,min_c=4d2,max_c=5d10,peaktime=8,omega=0d0,amplitude=-1.1d0, LOGSCALE=.true.)
    type(parametered_input) :: Jout_par  = parametered_input(sigma=2d0,min_c=1.d-4,max_c=5d0,peaktime=9,omega=0d0,amplitude=-1.1d0, LOGSCALE=.true.)
    ! real(dp) :: sigma, min_c, max_c, peaktime, omega, amplitude
    ! LOGICAL  :: LOGSCALE


    dt = 10
    time = 0
    sim_time=24.
  !Create/Open outup file netcdf
  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
open(111,FILE='output.dat')
write(111,'(5(a20))') 'c_acid', 'c_base', 'J_ACDC_NH3', 'c_dma', 'J_ACDC_DMA'

print FMT_Tm, time_hms(int(time))
  !Main loop time: Eulerian forward integration
  DO WHILE (time < sim_time*hour_s)

    !Photolysis
    !Chemistry

c_acid = NORMALD(time,MOD_H2SO4)*1d6
c_base = NORMALD(time,MOD_NH3)*1d6
  if (NUCLEATION) THEN

#ifdef ISACDC
  if (ACDC) CALL ACDC_J() ! SUBROUTINE in CONTAINS of this file
!  ELSE CALL SOME_OTHER_NUCLEATION_TYPE() ! SUBROUTINE in CONTAINS of this file
#endif

#ifndef ISACDC
  if (EVENMIN(time,15)) print*, NORMALD(time,Jout_par)
#endif


  END if ! NUCLEATION


    !Condensation
    !Coagulation
    !Deposition



    time = time + dt

    !Write output to file

  END DO	!Main loop time: Eulerian forward integration
  print FMT_Tm, time_hms(int(time))

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
! c_org:             Nucleating organic concentration [1/m3]. not in AFAIK
! cs_H2SO4:          Condensationnh3_to_ACDC sink of sulfuric acid [1/s/m3]
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
! c_org:             Nucleating organic concentration [1/m3]. not in AFAIK
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s/m3]
! TempK:             Temperature in Kelvins
! dt:                Main time step [s]
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_DMA:        Particle formation rate due to DMA [1/s/m3]
! acdc_cluster_diam: Outgrowing cluster diameter [m]. The cluster typically has 5 to 6 H2SO4 in it
! =================================================================================================

! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  if (c_base > 1d12 .or. J_ACDC_NH3 > 1d-6) THEN
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,IPR,dt,&
          ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (EVENMIN(time, 15)) print FMT_SUB, 'NH3 IGNORED'
  END IF

  ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  if (c_dma > 1d6 .or. J_ACDC_DMA > 1d-6) THEN
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,dt,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (EVENMIN(time, 15)) print FMT_SUB, 'DMA IGNORED'
  END IF

    if (EVENMIN(time, 15)) THEN
      print*,
      print FMT_Tm, time_hms(int(time))
      print FMT10_CVU,'ACID C: ', c_acid*1d-6, ' [1/cm3]'
      print FMT10_2CVU, 'NH3 C:', c_base*1d-6, ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
      print FMT10_2CVU, 'DMA C:', c_dma*1d-6 , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
      print FMT_LEND,
    END IF
    if (EVENMIN(time, 5)) write(111,'(5(es20.3))') c_acid, c_base, J_ACDC_NH3, c_dma, J_ACDC_DMA
  ! ................................................................................................
  ! END ACDC Nucleation
  ! ================================================================================================
END SUBROUTINE ACDC_J

SUBROUTINE SOME_OTHER_NUCLEATION_TYPE(J1, J2)
  REAL(dp) :: J1, J2
  J1 = 1d0
  J2 = 1d0
END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE

END PROGRAM SUPERMODEL
