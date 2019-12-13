module monomer_settings

use acdc_system, only : n_monomer_types, monomer_indices
implicit none

contains

subroutine sources_and_constants(nclust,source,isconst,fitted)
	implicit none
	
	real(kind(1.d0)) :: source(nclust)
	logical :: isconst(nclust)
	integer :: nclust, fitted(nclust,0:nclust), n_monomers(n_monomer_types), n, i
	
	source = 0.d0		! initialize all source terms to 0
	isconst = .false.	! initialize all concentrations to not be constants ...
	fitted = 0			! ... or fitted
	
	! get the cluster numbers corresponding to the monomers
	!call monomer_indices(n_monomers)

	! use constant concentrations for all monomers
	!isconst(n_monomers) = .true.


!   (Everything related to fitting or setting something constant commented out
!   ---> all concentrations are allowed to vary according to the equations)

end subroutine sources_and_constants


end module monomer_settings
