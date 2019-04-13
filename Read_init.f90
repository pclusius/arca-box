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
           Extra_data           ,&
           Current_case
NAMELIST /NML_Flag/ Aerosol_flag, chemistry_flag, particle_flag, Extra_data, Current_case


!!!!! init file name
character(len=100):: Fname_init ! init file names

!!!!! reading DMPS file, temp
character(len=100) :: DMPS_file, TEMP_file, VOC_file, TEMP_path
logical :: dmps_file_check, temp_file_check,voc_file_check

NAMELIST /NML_DMPS/ DMPS_file, dmps_file_check
NAMELIST /NML_TEMP/ TEMP_path, TEMP_file, temp_file_check
NAMELIST /NML_VOC / VOC_file, voc_file_check


contains

subroutine read_init_file

!!! reads the .init file !!!!1

integer:: io_status, ioinit, iofloat

CALL GETARG(1,Fname_init)

Write(*,*) 'Reading initializtion variables from file','',TRIM(ADJUSTL(Fname_init))

OPEN(UNIT=50, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD')

READ(50,NML=NML_Path,  IOSTAT= io_status) !directories and test cases
READ(50,NML=NML_Flag,  IOSTAT= ioinit) !flags
READ(50,NML=NML_DMPS,  IOSTAT= io_status) ! dmps_file information
READ(50,NML=NML_TEMP,  IOSTAT= io_status) ! temp_file information
READ(50,NML=NML_VOC,  IOSTAT= io_status) ! temp_file information


CLOSE(50)

if (ioinit /= 0) then
  backspace(ioinit)
  write(*,*) 'Exit error:: Check the flags'
  STOP
end if


end subroutine read_init_file



end module Read_init
