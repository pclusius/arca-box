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
    REAL(dp) :: c_acid =1e6*1d6
    REAL(dp) :: c_base =1d8*1d6
    REAL(dp) :: c_dma = 10 ! in [ppt]
    REAL(dp) :: c_org = 1d6*1d6
    REAL(dp) :: cs_H2SO4 = 0.003d0
    REAL(dp) :: TempK = K0 + 15d0
    REAL(dp) :: ION_RATE = 3d6
    REAL(dp) :: J_ACDC_NH3
    REAL(dp) :: J_ACDC_DMA
    REAL(dp) :: J_NH3_BY_IONS(3)
    REAL(dp) :: acdc_cluster_diam = 1.17d-9
    LOGICAL(dp) :: ACDC_solve_ss = .false.

    dt = 10
    time = 0
    sim_time=1
  !Create/Open outup file netcdf

  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
    print*, c_acid,c_base,c_org,cs_H2SO4,TempK,ION_RATE,dt,ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS


  !Main loop time: Eulerian forward integration
  DO WHILE (time < sim_time*hour_s)


    c_dma = 0.2*1d-12*C_AIR_m3(101325d0,TempK)
    c_base = 200*1d-12*C_AIR_m3(101325d0,TempK)
    print*, C_AIR_m3(101325d0,TempK)
    ! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!
    CALL get_acdc_J(c_acid,c_base,c_org,cs_H2SO4,TempK,ION_RATE,dt,ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
    ! NUCLEATION BY S-ACID AND DMA - NOTE: ingoing concentrations are assumed to be in 1/m3!!
    CALL get_acdc_D(c_acid,c_dma,c_org,cs_H2SO4,TempK,dt,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)

    print*, '#############################'
    print*, J_ACDC_NH3, J_NH3_BY_IONS, J_ACDC_DMA
    print*, time_hms(int(time))

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
