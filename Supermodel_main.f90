PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)


  !USE module statements
!  USE second_Precision,  ONLY : dp    ! KPP Numerical type
  USE Read_init
  USE Inputs

  IMPLICIT NONE


  !Variable declaration
  REAL ::       &
    dt,         &!integration time step [s]
    time         !simulation time [s]

  !Create/Open outup file netcdf

  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
 
  Call read_init_file 

 
 write(*,*), 'The flags are', Aerosol_flag
 write(*,*), 'The Working directory is', Work_dir
 write(*,*), 'The CASE directory is',    case_dir
 write(*,*), 'The Case name is',         case_name
 write(*,*), 'The Input file is',        Input_file

 
CALL read_input_data
  
  !Main loop time: Eulerian forward integration 
  DO WHILE (time < sim_time)

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

END PROGRAM 
