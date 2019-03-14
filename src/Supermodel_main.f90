PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)

!USE module statements
USE second_Precision,  ONLY : dp ! KPP Numerical type
USE constants
USE AUXILLARIES
#ifdef ACDC
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

    REAL(dp), dimension(:),allocatable :: modelled_SA
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
    LOGICAL :: ACDC_solve_ss = .false., NUCLEATION=.true.

    dt = 10
    time = 0
    sim_time=24.
  !Create/Open outup file netcdf
    ALLOCATE(modelled_SA(ceiling(sim_time)*6+1))

  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
open(111,FILE='output.dat')
write(111,'(5(a20))') 'c_acid', 'c_base', 'J_ACDC_NH3', 'c_dma', 'J_ACDC_DMA'

  !Main loop time: Eulerian forward integration
  DO WHILE (time < sim_time*hour_s)

    !Photolysis
    !Chemistry

  if (NUCLEATION) THEN
#ifdef ACDC
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
  c_dma  = 0!NORMALD_SINGLE(time,3d0, 4d4, 2.5d7, 11.d0)*1d6*0
  !c_base = NORMALD(time,6d0, 5d4, 2.5d9, 24*SIN(time/3600d0), .true.)*1d6
  c_base = NORMALD(time,3d0, 2d5, 9.5d5, 11.d0, .false.)*1d6
  c_acid = NORMALD(time,9d0, 2d5, 9.5d5, 11.d0, .true.)*1d6

  ! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  if (c_base > 1d-1 .or. J_ACDC_NH3 > 1d-6) THEN
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,IPR,dt,&
          ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (EVENMIN(time, 15)) print*, 'NH3 IGNORED'
  END IF

  ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!

  ! Speed up program by ignoring nucleation when there is none
  if (c_dma > 1d1 .or. J_ACDC_DMA > 1d-6) THEN
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,dt,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (EVENMIN(time, 15)) print*, 'DMA IGNORED'
  END IF

    if (EVENMIN(time, 15)) THEN
            print*,
            print'(a,a)', '* Time: ', time_hms(int(time))
            print'(a,es10.3,a)', '* ACID C:', c_acid*1d-6, '/cm3'
            print'(a,es10.3,a,es10.3,a)', '* NH3 C:', c_base*1d-6, '/cm3      DMA C:', c_dma*1d-6, '/cm3'
            print'(a,es10.3,a,es10.3,a)', '* J_NH3 :', J_ACDC_NH3*1d-6, '/cm3      J_DMA:', J_ACDC_DMA*1d-6, '/cm3'
    END IF
    if (EVENMIN(time, 5)) write(111,'(5(es20.3))') c_acid, c_base, J_ACDC_NH3, c_dma, J_ACDC_DMA
  ! ................................................................................................
  ! END ACDC Nucleation
  ! ================================================================================================
#endif
  END if ! NUCLEATION


    !Condensation
    !Coagulation
    !Deposition



    time = time + dt

    !Write output to file

  END DO	!Main loop time: Eulerian forward integration

  !Close output file netcdf

END PROGRAM SUPERMODEL
