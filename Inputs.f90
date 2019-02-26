!!!!!! Module contains all the inputs used in "Superbox-model"

module Inputs

USE Read_init


implicit none

 character(len=100):: data_dir

contains


subroutine read_input_data

implicit none

character(len=100)::  data_dir
integer :: io_status, i, j, line
!character(len=120):: line
integer, dimension(:,:) ,allocatable:: apin_val!, val2(:), val3(:), val4(:), val5(:), val6(:), val7(:), val8(:)

call read_init_file

data_dir =Trim(ADJUSTL(Case_dir)) //  Case_name  

OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' // Input_file,STATUS='OLD')

allocate(apin_val(rows,cols))

!!!! reading data into the variable
DO i=1,rows
 READ(51,*) (apin_val(i,j),j=1,cols)
 WRITE(*,*) (apin_val(i,j),j=1,cols)
End do 

end subroutine read_input_data

end module
