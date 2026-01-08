program RetR
use second_precision
use second_Monitor
use second_Parameters

implicit none
integer :: i
character(len=512) :: path
character(len=40)  :: prefix(9)
path = ''

prefix = [&
'TEMPK                                 ', &
'PRESSURE                              ', &
'REL_HUMIDITY                          ', &
'CONDENS_SINK                          ', &
'CON_SIN_NITR                          ', &
'SW_RADIATION                          ', &
'ION_PROD_RATE                         ', &
'NUC_RATE_IN                           ', &
'# DONT_CHANGE Anything ABOVE THIS LINE']

CALL GETARG(1,path)
if (TRIM(path) == '') STOP 'Give the path to chemistry directory as cmdline option'


OPEN(600,file=TRIM(path)//"/SPC_NAMES.txt",status='replace',action='write')
do i=1,size(SPC_NAMES,1)
  WRITE(600,'(a)') TRIM(SPC_NAMES(i))
end do
CLOSE(600)

OPEN(600,file=TRIM(path)//"/NAMES.dat",status='replace',action='write')
do i=1,9
WRITE(600,'(a)') TRIM(prefix(i))
end do

do i=1,size(SPC_NAMES,1)
  IF ((TRIM(SPC_NAMES(i)) /= 'DUMMY').and.(TRIM(SPC_NAMES(i)) /= 'dummy')) &
    WRITE(600,'(a)') TRIM(SPC_NAMES(i))
end do
CLOSE(600)

OPEN(600,file=TRIM(path)//"/EQN_NAMES.txt",status='replace',action='write')
do i=1,size(EQN_NAMES,1)
  WRITE(600,'(a)') TRIM(EQN_NAMES(i))
end do
CLOSE(600)
end program RetR
