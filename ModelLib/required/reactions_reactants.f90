program RetR
use second_precision
use second_Monitor
use second_Parameters

implicit none
integer :: i
character(len=512) :: path
path = ''

CALL GETARG(1,path)
if (TRIM(path) == '') STOP 'Give the path to chemistry directory as cmdline option'


OPEN(600,file=TRIM(path)//"/SPC_NAMES.txt",status='replace',action='write')
do i=1,size(SPC_NAMES,1)
  WRITE(600,*) TRIM(SPC_NAMES(i))
end do
CLOSE(600)

OPEN(600,file=TRIM(path)//"/EQN_NAMES.txt",status='replace',action='write')
do i=1,size(EQN_NAMES,1)
  WRITE(600,*) TRIM(EQN_NAMES(i))
end do
CLOSE(600)
end program RetR
