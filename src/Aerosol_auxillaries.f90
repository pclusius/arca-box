Module Aerosol_auxillaries

Use constants
Use SECOND_PRECISION, only:dp
Use AUXILLARIES

implicit NONE
PUBLIC

real(dp),dimension(3)::N_speed !! nucleation, condesnation,coagulation
real(dp) :: precision
logical :: first_time

contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed

subroutine
End module Aerosol_auxillaries
