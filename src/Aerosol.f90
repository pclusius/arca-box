Module Aerosol

Use CONSTANTS
USE INPUT
Use SECOND_PRECISION, only:dp
Use AUXILLARIES
USE Aerosol_auxillaries

implicit NONE



Contains

Subroutine condensation_routine(values,timestep, end_time, setup, particle,ambient,vapours)
!!! This contains the condensation scheme
implicit NONE

 real(dp),intent(in) :: timestep, end_time

 type(aerosol_setup),intent(in)       :: setup
 type(particle_properties),intent(in) :: particle
 type(ambient_properties),intent(in)  :: ambient
 type(vapour_ambient)  :: vapours
 logical :: first_call=.True.

 !!! Vector containing important output VARIABLES
 real(dp),dimension(:), allocatable,intent(inout), target :: values!(setup%n_cond+setup%sections+setup%sections*setup%n_cond) !INTENT(in)
 real(dp), pointer, dimension(:) :: n_conc, vap_conc
 real(dp), pointer, dimension(:,:) :: vol_conc

 !!!!!!! Local VARIABLES
 !!! mole_fraction(sections, number of condensable compounds), volume concetration
 !!! gas phase vapour conc , vap_conc_gp(number of condesable compounds), vapour gas phase total
 real(dp),dimension(setup%sections,setup%n_cond) :: mole_fractions
 real(dp),dimension(setup%n_cond) :: vap_conc_gp,vap_total, vapor_volume
 real(dp) :: organics, time

 !!!!!! Output of Subroutines call
 !!! CR = Collision rate, KT = Kelvin-Taylor formulae for Kelvin effect
  real(dp),dimension(setup%sections) :: CR, KT, bin_volume !! CR = Collision rate, KT = Kelvin-Taylor formulae for Kelvin effect

!!! intiailize
   time           = 0.


   if (first_call) THEN
     first_call=.False.
     allocate(values(setup%n_cond+setup%sections+setup%sections*setup%n_cond))
     print*,'in first_call'
     print*,sum(vol_conc,1)
   end if

  !Associate array segments (pointers) with the correct variables (target)
  vap_conc(1:setup%n_cond) => values(1:setup%n_cond)
  n_conc(1:setup%sections) => values(setup%n_cond+1:setup%n_cond+setup%sections)
  vol_conc(1:setup%sections, 1:setup%n_cond) => values(setup%n_cond+setup%sections:)

 if(time <=end_time) THEN
     !!! total vapour concentrations
     vap_total = vap_conc + sum(vol_conc,1)*um3_to_m3/vapours%molec_volume ! # m^-3
     mole_fractions = 0.
     ! print*, 'Here in loop',sum(vol_conc,1)!vapours%molec_volume
     bin_volume = sum(vol_conc,2)
     vapor_volume = sum(vol_conc, 1)

 end if




end subroutine condensation_routine

End module Aerosol
