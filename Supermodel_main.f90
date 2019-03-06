PROGRAM Supermodel
!We keep things easy and clear (no goto commands,...)


  !USE module statements
  USE second_Precision,  ONLY : dp    ! KPP Numerical type
  USE Read_init
  USE Inputs
  USE Constants

  IMPLICIT NONE


  !Variable declaration

integer   :: counter
real (dp) :: time_start
real,dimension(:,:), allocatable:: apin_value,dmps_value
!!! none at the moment
 !Create/Open outup file netcdf

  !Variable initialization
    !Species properties in gas/particle phase
    !Gas phase
    !Aerosol phase(distribution & composition)
    !Boundary conditions(dilution, losses, light,...)
 
  Call read_init_file 
  CALL read_input_data(apin_value,dmps_value)

counter = 0 
! write(*,*), 'The flags in order are', Aerosol_flag, Chemistry_flag, particle_flag, extra_data
! write(*,*), 'The Working directory is','', Work_dir
! write(*,*), 'The CASE directory is','',    case_dir
! write(*,*), 'The Case name is','',         case_name
! write(*,*), 'The Input file is', '',       Input_file
 

write(*,*), dmps_value
  !Main loop time: Eulerian forward integration 
DO WHILE (ctime < sim_time)

    !Photolysis
    !Chemistry
 if (Chemistry_flag) then 
 !    WRITE(*,*), 'Execute chemistry part here'
 end if

   ! Aerosol part

 if (Aerosol_flag) then 
 !    WRITE(*,*), 'Execute Aerosol part here'
 
 
    !Nucleation
    !Condensation
    !Coagulation
    !Deposition
 end if





   ctime = ctime + dt
   counter= counter+1  
    !Write output to file
 ! write(*,*), counter, 'time steps'
  END DO	!Main loop time: Eulerian forward integration 

  !Close output file netcdf

END PROGRAM 
