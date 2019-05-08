Module aerosol_auxillaries

Use constants
Use SECOND_PRECISION, only:dp
Use AUXILLARIES


implicit NONE
PUBLIC


real(dp),dimension(3)::N_speed !! nucleation, condensation,coagulation
real(dp) :: precision
logical :: first_time

!!! This datatype contains all parameters for input vapours
type :: vapour_ambient
  real(dp), allocatable :: molar_mass(:), parameter_A(:), parameter_B(:)
  character(len=256), allocatable :: vapour_names(:)
  real(dp) :: alpha        = 1.0
  real(dp) :: density      = 1400.0
  real(dp) :: surf_tension = 0.05
  integer :: vapour_number
  integer :: vbs_bins
end type vapour_ambient


contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed


End module aerosol_auxillaries
