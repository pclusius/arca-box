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
  real(dp) :: alpha        = 1.0
  real(dp),allocatable  :: density(:)
  real(dp), allocatable :: surf_tension(:)
  integer  :: vapour_number
  integer  :: vbs_bins
  real(dp),allocatable :: molec_mass(:), molec_volume(:) !!! molecule mass and molecule volume
  real(dp),allocatable :: c_sat(:), vap_conc(:)!, vapout_type(:), condensing_type(:)
end type vapour_ambient

type aerosol_setup
!!!! Contains setup for size bins, number of condensable compounds, distribution method (moving avg fully stationary, moving sections )
 integer             ::  sections,   &    !! size sections used
                         n_noncond,  &    !! number of non-condensable species
                         n_cond,     &    !! number of condensable species
                         tot_spec,   &    !!! total number of species, cond+noncond
                         init_modes       !!! initial particle modes

 character(len=4)    ::  dist_type
 logical             ::  Use_Raoult, quasistationary
end type aerosol_setup

type particle_properties
!!! Contains particle properties such as particle number_concentration
!!! the following n_conc ----> orifinal radius is of dimension sections (dimension(sections))
!!! volume concentration is of dimension (sections, n_noncond+n_cond) ---> (dimension(sections, n_noncond+n_cond))
!!! order is of size sections (dimension(sections))
 real(dp),dimension(:), allocatable :: n_conc,       &       !!! number conc
                                       radius,       &       !!! radius
                                       rdry,         &       !!! dry radius
                                       rdry_orig,    &       !!! dry radius original
                                       mass,         &       !!! mass in sections
                                       core,         &       !!! core volume in sections
                                       gr,           &       !!! growth rate
                                       orig_radius           !!! original radius

 real(dp),dimension(:,:), allocatable :: vol_conc ! um^3/m^3

 integer, dimension(:), allocatable   :: order !!! used for rearranging in bins


end type particle_properties


type ambient_properties
  !!! contains variables related to vapors and ambient conditions
   real(dp) ::         temp,      &
                       pres,      &
                       rh
end type ambient_properties

type initial_distribution
!!!contains parameters for initial particle size distribution
!!! n_mode, sigma_mode, r_mode is of dimension(init_mode)
!!! mfractions is of dimension tot_spec
    real(dp) :: r_min,r_max                                       !!! dry minimum and maximum radius
    real(dp),dimension(:),allocatable :: n_mode,sigma_mode,r_mode   !!! dimension(init mode) mode conc, mode deivation and radius of mode
    real(dp),dimension(:),allocatable :: mfractions                 !!! dimension(tot_spec) mole fractions
end type initial_distribution


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


subroutine Aerosol_intialization(AER_setup_mod, AER_par_prop_mod, vapours_mod,AER_init_dist_mod, ambient_mod)
  !!!This subroutine allocates dimensions and setup options such as distribution, number and radius modes, section size
  !!!  Setup number of sections and condesable and non-condesable compounds
  IMPLICIT none
  type(aerosol_setup)       :: AER_setup_mod
  type(particle_properties) :: AER_par_prop_mod
  type(vapour_ambient)      :: vapours_mod
  type(initial_distribution):: AER_init_dist_mod
  type(ambient_properties)  :: ambient_mod

 !!!!!! Setup aerosol
  AER_setup_mod%sections        = 40
  AER_setup_mod%n_noncond       = 0
  AER_setup_mod%n_cond          = vapours_mod%vbs_bins
  AER_setup_mod%tot_spec        = AER_setup_mod%n_noncond + AER_setup_mod%n_cond
  AER_setup_mod%dist_type       = 'MVFS'
  AER_setup_mod%Use_Raoult      = .True.
  AER_setup_mod%quasistationary = .True.
  AER_setup_mod%init_modes      = 2

  !!! Now allocate particle properties
  allocate(AER_par_prop_mod%n_conc(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%radius(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%rdry(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%rdry_orig(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%core(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%mass(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%gr(AER_setup_mod%sections))
  allocate(AER_par_prop_mod%orig_radius(AER_setup_mod%sections))

  allocate(AER_par_prop_mod%order(AER_setup_mod%sections))

  allocate(AER_par_prop_mod%vol_conc(AER_setup_mod%sections, AER_setup_mod%tot_spec))
  print * ,'Number of sections available for  particle concentration is=', size(AER_par_prop_mod%n_conc)
  print * ,'Number of sections x compounds for volume concentration is=', shape(AER_par_prop_mod%vol_conc)

  !!! Now allocating the initial distribution
  allocate(AER_init_dist_mod%n_mode(AER_setup_mod%init_modes))
  allocate(AER_init_dist_mod%sigma_mode(AER_setup_mod%init_modes))
  allocate(AER_init_dist_mod%r_mode(AER_setup_mod%init_modes))
  allocate(AER_init_dist_mod%mfractions(AER_setup_mod%tot_spec))

  print * ,'Size of n_modes is =',     size(AER_init_dist_mod%n_mode)
  print * ,'Size of mfractions is =',  shape(AER_init_dist_mod%mfractions)

  AER_init_dist_mod%r_min     = 0.5d-9
  AER_init_dist_mod%r_max     = 1.0d-5
  AER_init_dist_mod%r_mode    = (/2.0d-8,1.0d-7/)/2.0
  AER_init_dist_mod%n_mode    = (/1.0d7, 1.0d7/)
  AER_init_dist_mod%sigma_mode= (/1.5, 1.3/)

  !!! ambient Properties
  ambient_mod%temp = 293.15 !! K
  ambient_mod%rh   = 60.0   !! %

  !!! saturation vapour pressure
  do ii = 1, AER_setup_mod%tot_spec
   vapours_mod%c_sat(ii) = calculate_saturation_vp(vapours_mod%parameter_A(ii),vapours_mod%parameter_B(ii), ambient_mod%temp)
  end do

end SUBROUTINE Aerosol_intialization

End module aerosol_auxillaries
