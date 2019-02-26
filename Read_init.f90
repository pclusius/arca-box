!!!! module to read init file

Module Read_init

Implicit none
!Private

!!!! path names
character(len=100):: WORK_DIR
character(len=100):: CASE_DIR
character(len=100):: CASE_NAME
character(len=100):: INPUT_FILE
NAMELIST /NML_Path/ Work_dir, Case_Dir, Case_name, Input_file

!!!!! Flags needed
Integer :: Aerosol_flag
Integer :: Chemistry_flag
Integer :: Particle_flag
Integer :: Extra_data 
NAMELIST /NML_Flag/ Aerosol_flag, chemistry_flag, particle_flag, Extra_data

!!!!! integration timesteps
REAL :: dt
NAMELIST /NML_Time/ dt

!!!!! file shapes
integer :: rows, cols
NAMELIST /NML_Shape/ rows, cols

!!!!! init file name
character(len=100):: Fname_init ! init file names


!!!! check
integer :: ioint, iofloat


contains

subroutine read_init_file

!!! reads the .init file !!!!1


integer:: io_status, ioinit, iofloat

CALL GETARG(1,Fname_init)

Write(*,*) 'Reading initializtion variables from file',TRIM(ADJUSTL(Fname_init))

OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD')

READ(50,NML=NML_Path,  IOSTAT= io_status) !directories and test cases
READ(50,NML=NML_Flag,  IOSTAT= ioint) !flags
READ(50,NML=NML_Time,  IOSTAT= iofloat) !flags
READ(50,NML=NML_SHAPE, IOSTAT= io_status) !flags

CLOSE(50)

if (ioint /= 0 .or. iofloat/=0) then
 write(*,*) 'Exit error:: Check the flags'
 STOP
end if 



end subroutine read_init_file



end module Read_init
