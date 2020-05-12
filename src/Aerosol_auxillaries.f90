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


contains

subroutine set_speed()

if (First_time) THEN
  N_speed(:) = 1.0
  First_time = .False.
end if

end subroutine set_speed


! =====================================================================================================
! Calculate molecular mass in KG
! input molar_mass
! .....................................................................................................
! pure elemental function calculate_molecular_mass(molar_mass) result(mass)
!   real(dp), intent(in) :: molar_mass
!   real(dp) :: mass
!   mass = molar_mass / Na ! kg/molecule
! end function calculate_molecular_mass

! calculate molecular volume
pure elemental function calculate_molecular_volume(density, molecule_mass) result(volume)
  real(dp), intent(in) :: molecule_mass, density
  real(dp) :: volume
  volume = molecule_mass / density !!
end function calculate_molecular_volume



SUBROUTINE Calculate_SaturationVapourConcentration(VAPOUR_PROP,TEMPK, MLT)
  IMPLICIT NONE
  type(vapour_ambient), INTENT(INOUT) :: VAPOUR_PROP
  real(dp), INTENT(IN) :: TEMPK
  real(dp), INTENT(IN) :: MLT
  ! Update saturation concentrations to current temperature, omit sulfuric acid
  VAPOUR_PROP%c_sat(1:VAPOUR_PROP%vapour_number) =  calculate_saturation_vp( &
  VAPOUR_PROP%parameter_A(1:VAPOUR_PROP%vapour_number),&
  VAPOUR_PROP%parameter_B(1:VAPOUR_PROP%vapour_number),&
  TEMPK) * MLT

END SUBROUTINE Calculate_SaturationVapourConcentration

! calculate saturation vapour PRESSURE
! input parameter_A, parameter_B and temperature
pure elemental function calculate_saturation_vp(A,B, Temperature) result(Vapour_concentration)
  real(dp), intent(in) :: A, B, temperature
  real(dp) :: Vapour_concentration, vapour_pressure

  ! Using antoine equation log_10(p) = A- (B/T)
  vapour_pressure      = 10 ** (A - (B/temperature)) ! in atm
  Vapour_concentration = (vapour_pressure*101325)/(kb * temperature) ! #/m3

end function calculate_saturation_vp


End module aerosol_auxillaries
