!!!!!! Module contains all the inputs used in "Superbox-model"

module Inputs
USE second_precision
USE Read_init
Use constants
Use auxillaries

implicit none

 !character(len=100):: data_dir
 real, dimension(:,:) ,allocatable:: apin_val,dmps_data,temp_data,voc_data
 character(len=100)::  data_dir, temp_dir
 integer :: io_status, i, j, line
 type(nrowcol) :: rowcol_count

contains

!!  assigning the allocatables and reading data

subroutine read_input_data()!(apin_val,dmps_data,temp_data,voc_data)

call read_init_file

data_dir =Trim(ADJUSTL(Case_dir)) //  Case_name
temp_dir =Trim(Adjustl(Temp_path))

OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' // Input_file,STATUS='OLD')
rowcol_count%rows = ROWCOUNT(51)
rowcol_count%cols = COLCOUNT(51)

allocate(apin_val(rowcol_count%rows,rowcol_count%cols))

!!!! reading data into the variable
DO i=1,rowcol_count%rows
 READ(51,*) (apin_val(i,j),j=1,rowcol_count%cols)
End do

call clear_values() !!! clears row and colum variables

!! check if dmps data is used or not. If no then do nothing
if (dmps_file_check) then
 write(*,*),'DMPS FILE PRESENT'
 OPEN(unit=52, File=TRIM(ADJUSTL(data_dir)) // '/' // DMPS_file ,STATUS='OLD')
 rowcol_count%rows = ROWCOUNT(52)
 rowcol_count%cols = COLCOUNT(52)
 allocate(dmps_data(rowcol_count%rows,rowcol_count%cols))
 call clear_values()
end if


end subroutine read_input_data

!!! clears the values allcoated
subroutine clear_values()
 rowcol_count%rows=0
 rowcol_count%cols=0
end subroutine clear_values

end module
