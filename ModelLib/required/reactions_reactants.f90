program RetR
use second_precision
use second_Monitor
use second_Parameters

implicit none
integer :: i

OPEN(600,file="SPC_NAMES.txt",status='replace',action='write')
do i=1,size(SPC_NAMES,1)
  WRITE(600,*) TRIM(SPC_NAMES(i))
end do
CLOSE(600)

OPEN(600,file="EQN_NAMES.txt",status='replace',action='write')
do i=1,size(EQN_NAMES,1)
  WRITE(600,*) TRIM(EQN_NAMES(i))
end do
CLOSE(600)
end program RetR
