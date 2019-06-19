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
  integer  :: vapour_number
  integer  :: vbs_bins
  real(dp),allocatable :: molec_mass(:), molec_volume(:) !!! molecule mass and molecule volume
end type vapour_ambient

type aerosol_setup
!!!! Contains setup for size bins, number of condensable compounds, distribution method (moving avg fully stationary, moving sections )
 integer             ::  sections,   &    !! size sections used
                         n_noncond,  &    !! number of non-condensable species
                         n_cond,     &    !! number of condensable species
                         tot_spec          !!! total number of species, cond+noncond

 character(len=4)    ::  dist_type
end type aerosol_setup
!
! type particle_properties
! !!! Contains particle properties such as particle number_concentration
!  real(dp) ::

! end type particle_properties

contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed

!!!! Calculate molecular mass in KG
!!! input molar_mass
pure elemental function calculate_molecular_mass(molecular_mass) result(mass)
  real(dp), intent(in) :: molecular_mass
  real(dp) :: mass
  mass = molecular_mass / Na *1D-3 !! convert to kg/#
end function calculate_molecular_mass

!!! calculate molecular volume
pure elemental function calculate_molecular_volume(density, molecule_mass) result(volume)
  real(dp), intent(in) :: molecule_mass, density
  real(dp) :: volume
  volume = molecule_mass / density !!
end function calculate_molecular_volume


!!!! calculate saturation vapour PRESSURE
!!! input parameter_A, parameter_B and temperature
pure elemental function calculate_saturation_vp(A,B, Temperature) result(Vapour_concentration)
  real(dp), intent(in) :: A, B, temperature
  real(dp) :: Vapour_concentration, vapour_pressure

!Using antoine equation log_10(p) = A- (B/T)
  vapour_pressure      = 10 ** (A - (B/temperature)) ! in atm
  Vapour_concentration = (vapour_pressure*101325)/(kb * temperature) ! #/m3

end function calculate_saturation_vp

End module aerosol_auxillaries
