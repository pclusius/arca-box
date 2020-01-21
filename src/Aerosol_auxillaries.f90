Module aerosol_auxillaries

Use constants
Use SECOND_PRECISION, only:dp
Use AUXILLARIES

implicit NONE
PUBLIC

real(dp),dimension(3)::N_speed !! nucleation, condensation,coagulation
real(dp) :: precision
logical :: first_time
integer :: ii

!!! This datatype contains all parameters for input vapours
type :: vapour_ambient
  real(dp), allocatable :: molar_mass(:), parameter_A(:), parameter_B(:)
  character(len=256), allocatable :: vapour_names(:)
  real(dp),allocatable :: alpha(:) !       = 1.0
  real(dp),allocatable  :: density(:)
  real(dp),allocatable :: surf_tension(:)
  integer  :: vapour_number
  integer  :: vbs_bins
  integer,allocatable  :: cond_type(:)
  real(dp),allocatable :: molec_dia(:)
  real(dp),allocatable :: molec_mass(:), molec_volume(:) !!! molecule mass and molecule volume
  real(dp),allocatable :: c_sat(:), vap_conc(:)!, vapour_type(:), condensing_type(:)
  real(dp),allocatable :: mfractions(:)        !!! dimension(tot_spec) mole fractions
end type vapour_ambient


type ambient_properties
  !!! contains variables related to vapors and ambient conditions
   real(dp) ::  temperature = 293.15
   real(dp) ::  pressure = 1D5
   real(dp) ::  rh   = 60.0
end type ambient_properties

type atoms  ! for reading in molar mass of each atom. WIll be used to calculate diffusion
  real(dp), allocatable :: N_Carbon(:)
  real(dp), allocatable :: N_Oxygen(:)
  real(dp), allocatable :: N_Hydrogen(:)
  real(dp), allocatable :: N_Nitrogen(:)
  REAL(dp), allocatable :: comp_prop(:,:)
end type atoms





contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed

!!!! Calculate molecular mass in KG
!!! input molar_mass
pure elemental function calculate_molecular_mass(molar_mass) result(mass)
  real(dp), intent(in) :: molar_mass
  real(dp) :: mass
  mass = molar_mass / Na !*1D-3 !! convert to kg/#
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
