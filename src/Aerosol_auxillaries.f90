Module Aerosol_auxillaries

Use constants
Use SECOND_PRECISION, only:dp
Use AUXILLARIES
USE INPUT

implicit NONE
PUBLIC


real(dp),dimension(3)::N_speed !! nucleation, condensation,coagulation
real(dp) :: precision
logical :: first_time



!print*, vap_details%vapour_number, '', vap_details%vbs_bins


contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed

!subroutine
End module Aerosol_auxillaries
