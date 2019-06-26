Module Aerosol

Use CONSTANTS
USE INPUT
Use SECOND_PRECISION, only:dp
Use AUXILLARIES
USE Aerosol_auxillaries

implicit NONE


Contains

Subroutine condensation_routine(values,timestep, end_time, setup, particle,ambient)
!!! This contains the condensation scheme
implicit NONE

 real(dp), dimension(:), allocatable:: Values
 real(dp),intent(in) :: timestep, end_time

 allocate(values(vapours%vapour_number))

 type(aerosol_setup),intent(in)       :: setup
 type(particle_properties),intent(in) :: particle
 type(ambient_properties),intent(in)  :: ambient

 print*,'In condensation maybe', vapours%vapour_number
 print*,'In condensation', size(vapours%molar_mass)

end subroutine condensation_routine

End module Aerosol
