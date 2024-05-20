program save_reactions
  use second_Monitor, only: EQN_NAMES
  integer :: i
  real :: RCONST(size(EQN_NAMES))
  CHARACTER(len=100) :: TEXT(size(EQN_NAMES)) = ''
  RCONST = -1.0
! Do not change anything in this file, it is a template used by create_chemistry
! Do not change anything in this file, it is a template used by create_chemistry
! Do not change anything in this file, it is a template used by create_chemistry
! Do not change anything in this file, it is a template used by create_chemistry
! Do not change anything in this file, it is a template used by create_chemistry
  print*, '  order   index   constant_R                       reaction'
  do i=1,size(EQN_NAMES)
  if (RCONST(i)>-1.0) print'(2i8,es12.3,a)',0,i,RCONST(i), EQN_NAMES(i)
  if (RCONST(i)< 0.0) print'(2i8,a12,a,a)',0,i,'    -   ', EQN_NAMES(i),TEXT(i)
  end do
end program save_reactions
