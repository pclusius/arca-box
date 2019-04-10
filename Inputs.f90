!!!!!! Module contains all the inputs used in "Superbox-model"

module Inputs
USE second_precision
USE Read_init
Use constants
Use auxillaries

implicit none

 character(len=100):: data_dir
 real, dimension(:,:) ,allocatable:: apin_val,dmps_data,temp_data,voc_data


contains

subroutine read_input_data(apin_val,dmps_data,temp_data,voc_data)
implicit none

character(len=100)::  data_dir, temp_dir
integer :: io_status, i, j, line
real, dimension(:,:) ,allocatable,intent(out):: apin_val,dmps_data,temp_data,voc_data
integer :: nrows, ncols!,dmps_rows,dmps_cols, temp_rows,temp_cols

call read_init_file

data_dir =Trim(ADJUSTL(Case_dir)) //  Case_name
temp_dir =Trim(Adjustl(Temp_path))

OPEN(unit=51, File=TRIM(ADJUSTL(data_dir)) // '/' // Input_file,STATUS='OLD')
nrows = ROWCOUNT(51)
ncols = COLCOUNT(51)

allocate(apin_val(nrows,ncols))


!!!! reading data into the variable
DO i=1,nrows
 READ(51,*) (apin_val(i,j),j=1,ncols)
End do

!! check if dmps data is used or not. If no then do nothing
print*, dmps_file_check, '', dmps_rows,'',nrows,'',ncols
if (dmps_file_check) then
 write(*,*),'DMPS FILE PRESENT'
 OPEN(unit=52, File=TRIM(ADJUSTL(data_dir)) // '/' // DMPS_file ,STATUS='OLD')
 dmps_rows = ROWCOUNT(52)
 dmps_cols = COLCOUNT(52)
 allocate(dmps_data(dmps_rows,dmps_cols))
end if

! if (temp_file_check) then
!  write(*,*), 'Temperature FILE PRESENT'
!  OPEN(unit=53, File=TRIM(ADJUSTL(temp_dir)) // TEMP_file ,STATUS='OLD')
!  allocate(temp_data(temp_rows,temp_cols))
!  print*, 'temp read in'
! end if
!
! if (voc_file_check) then
!  write(*,*), 'VOC  FILE PRESENT'
!  OPEN(unit=54, File=TRIM(ADJUSTL(temp_dir)) // VOC_file ,STATUS='OLD')
!  allocate(voc_data(voc_rows,voc_cols))
! end if
!
!
!  DO i=1,dmps_rows
!   READ(52,*) (dmps_data(i,j),j=1, dmps_cols)
!  END DO
!
!  DO i=1,temp_rows
!   READ(53,*) (temp_data(i,j),j=1, temp_cols)
!  END DO
!
!  DO i=1,voc_rows
!   READ(54,*) (voc_data(i,j),j=1, voc_cols)
!  END DO


end subroutine read_input_data

end module
