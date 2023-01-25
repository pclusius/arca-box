program save_reactions
  use second_Monitor, only: EQN_NAMES
  integer :: i
  print*, '       order       index                         reaction'
  do i=1,size(EQN_NAMES)
	print*,1,i, EQN_NAMES(i)
  end do
end program save_reactions
