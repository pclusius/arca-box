!!!! module to read init file

Module Read_init

USE second_precision

Implicit none
!Private

!!!! path names
character(len=100):: WORK_DIR
character(len=100):: CASE_DIR
character(len=100):: CASE_NAME
character(len=100):: INPUT_FILE
NAMELIST /NML_Path/ Work_dir, Case_Dir, Case_name, Input_file

!!!!! Flags needed
Logical :: Aerosol_flag         ,&
           Chemistry_flag       ,&
           Particle_flag        ,& 
           Extra_data 
NAMELIST /NML_Flag/ Aerosol_flag, chemistry_flag, particle_flag, Extra_data

!!!!! integration timesteps
REAL(dp) :: dt           ,& ! integration time step (sec)
            sim_time     ,& ! total simulation time (sec)
            ctime           ! current time (sec)
NAMELIST /NML_Time/ dt, sim_time, ctime

!!!!! file shapes
integer :: rows, cols
NAMELIST /NML_Shape/ rows, cols

!!!!! init file name
character(len=100):: Fname_init ! init file names

!!!!! reading DMPS file 
character(len=100) :: DMPS_file
logical :: dmps_file_check
integer :: dmps_rows, dmps_cols
NAMELIST /NML_DMPS/ DMPS_file, dmps_file_check, dmps_rows, dmps_cols


contains

subroutine read_init_file

!!! reads the .init file !!!!1


integer:: io_status, ioinit, iofloat

CALL GETARG(1,Fname_init)

Write(*,*) 'Reading initializtion variables from file',TRIM(ADJUSTL(Fname_init))

OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD')

READ(50,NML=NML_Path,  IOSTAT= io_status) !directories and test cases
READ(50,NML=NML_Flag,  IOSTAT= ioinit) !flags
READ(50,NML=NML_Time,  IOSTAT= iofloat) !time
READ(50,NML=NML_SHAPE, IOSTAT= io_status) !shape
READ(50,NML=NML_DMPS,  IOSTAT= io_status) ! dmps_file information

CLOSE(50)

if (ioinit /= 0) then
  backspace(ioinit)
  write(*,*) 'Exit error:: Check the flags'
  STOP
elseif (iofloat /=0) then
  backspace(iofloat)
  write(*,*) 'Exit error:: Check the floating point number'
  Stop
end if 
 

end subroutine read_init_file



end module Read_init
