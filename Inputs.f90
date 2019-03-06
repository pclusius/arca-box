!!!!!! Module contains all the inputs used in "Superbox-model"

module Inputs

USE Read_init


implicit none

 character(len=100):: data_dir
 real, dimension(:,:) ,allocatable:: apin_val,dmps_data


contains

subroutine read_input_data(apin_val,dmps_data)

implicit none

character(len=100)::  data_dir
integer :: io_status, i, j, line
real, dimension(:,:) ,allocatable,intent(out):: apin_val,dmps_data


call read_init_file

data_dir =Trim(ADJUSTL(Case_dir)) //  Case_name  

OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' // Input_file,STATUS='OLD')
OPEN(unit=52, File=TRIM(ADJUSTL(data_dir)) // '/' // DMPS_file ,STATUS='OLD')

allocate(apin_val(rows,cols))


!!!! reading data into the variable
DO i=1,rows
 READ(51,*) (apin_val(i,j),j=1,cols)
 !WRITE(*,*) (apin_val(i,j),j=1,cols)
End do 

!! check if dmps data is used or not. If no then do nothing
if (dmps_file_check) then
 write(*,*), 'DMPS FILE PRESENT'
 allocate(dmps_data(dmps_rows,dmps_cols))
 
 DO i=1,dmps_rows
  READ(52,*) (dmps_data(i,j),j=1, dmps_cols)
 END DO


end if

end subroutine read_input_data

end module
