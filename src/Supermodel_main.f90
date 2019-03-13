PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)


  !USE module statements
  USE second_Precision,  ONLY : dp, ik    ! KPP Numerical type
  USE constants
  USE AUXILLARIES

  IMPLICIT NONE


  !Variable declaration
  REAL(dp) ::			&
    dt,				&	!integration time step [s]
    time,     & !simulation time [s]
    sim_time
<<<<<<< HEAD
    REAL(dp) :: c_acid =1e7*1d6
    REAL(dp) :: c_base =1d8*1d6
    REAL(dp) :: c_dma
    REAL(dp) :: c_org = 1d12*1d6
    REAL(dp) :: cs_H2SO4 = 0.000d6
=======
    REAL(dp) :: c_acid =1e6*1d6
    REAL(dp) :: c_base =1d8*1d6
    REAL(dp) :: c_dma = 10 ! in [ppt]
    REAL(dp) :: c_org = 1d6*1d6
    REAL(dp) :: cs_H2SO4 = 0.003d0
>>>>>>> bdd8694825940d04cd97e2446bdf98e234def27d
    REAL(dp) :: TempK = K0 + 15d0
    REAL(dp) :: ION_RATE = 3d6
    REAL(dp) :: J_ACDC_NH3
    REAL(dp) :: J_ACDC_DMA
    REAL(dp) :: J_NH3_BY_IONS(3)
    REAL(dp) :: acdc_cluster_diam = 1.17d-9
    LOGICAL(dp) :: ACDC_solve_ss = .false.

    dt = 10
    time = 0
<<<<<<< HEAD
    sim_time=30
=======
    sim_time=1
>>>>>>> bdd8694825940d04cd97e2446bdf98e234def27d
  !Create/Open outup file netcdf

  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
<<<<<<< HEAD
=======
    print*, c_acid,c_base,c_org,cs_H2SO4,TempK,ION_RATE,dt,ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS

>>>>>>> bdd8694825940d04cd97e2446bdf98e234def27d

  !Main loop time: Eulerian forward integration
  DO WHILE (time < sim_time*hour_s)


<<<<<<< HEAD
  ! ================================================================================================
  ! Call for Nucleation. See that all values entering here are in CUBIC METERS, KELVINS AND PASCALS.
  ! ................................................................................................
  ! Input for get_acdc_J:
  ! c_acid: Sulfuric acid concentration [1/m3]
  ! c_base: base (ammonia) concentration [1/m3]
  ! c_org: nucleating organic concentration [1/m3]
  ! cs_H2SO4: Condensation sink of sulfuric acid [1/s/m3]
  ! TempK: temperature in Kelvins
  ! ION_RATE: Ion production rate in ion pairs per second [1/m3/s]. 3d6 is a good guestimate
  ! dt: Main time step [s]
  ! ACDC_solve_ss: Solve steady state or only to timestep duration (generally makes no difference)
  ! J_ACDC_NH3: particle formation rate due to ammonia [1/s/m3]. Sum of J_NH3_BY_IONS
  ! acdc_cluster_diam: outgrowing cluster diameter [m]
  ! J_NH3_BY_IONS: particle formation rate by ions (neutral, negative and positive) [1/s/m3]
  ! ................................................................................................
  ! Input for get_acdc_D:
  ! c_acid: Sulfuric acid concentration [1/m3]
  ! c_dma: DMA concentration [1/m3]
  ! c_org: nucleating organic concentration [1/m3]
  ! cs_H2SO4: Condensation sink of sulfuric acid [1/s/m3]
  ! TempK: temperature in Kelvins
  ! dt: Main time step [s]
  ! ACDC_solve_ss: Solve steady state or only to timestep duration (generally makes no difference)
  ! J_ACDC_DMA: particle formation rate due to DMA [1/s/m3]
  ! acdc_cluster_diam: outgrowing cluster diameter [m]
  ! ================================================================================================
    c_dma = 10*1d-12*C_AIR_m3(101325d0,TempK)
    c_base = 200*1d-12*C_AIR_m3(101325d0,TempK)


=======
    c_dma = 0.2*1d-12*C_AIR_m3(101325d0,TempK)
    c_base = 200*1d-12*C_AIR_m3(101325d0,TempK)
    print*, C_AIR_m3(101325d0,TempK)
>>>>>>> bdd8694825940d04cd97e2446bdf98e234def27d
    ! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,ION_RATE,dt,ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
    ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,dt,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)

<<<<<<< HEAD
    print*, 'ACID C:', c_acid*1d-6, '/cm3      DMA C:', c_dma*1d-6, '/cm3'
    print*, 'J_NH3:', J_ACDC_NH3*1d-6, '/cm3      J_DMA:', J_ACDC_DMA*1d-6, '/cm3'
    print*, time_hms(int(time))
    print*, '#############################'
  ! ................................................................................................
  ! END Nucleation
  ! ================================================================================================
=======
    print*, '#############################'
    print*, J_ACDC_NH3, J_NH3_BY_IONS, J_ACDC_DMA
    print*, time_hms(int(time))
>>>>>>> bdd8694825940d04cd97e2446bdf98e234def27d

    !Photolysis
    !Chemistry
    !Nucleation
    !Condensation
    !Coagulation
    !Deposition



    time = time + dt

    !Write output to file

  END DO	!Main loop time: Eulerian forward integration

  !Close output file netcdf

END PROGRAM SUPERMODEL
