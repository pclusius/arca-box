module monomer_settingsD

use acdc_systemD, only : neq, nclust, n_mol_types, monomer_indices
use acdc_systemD, only : n1A, n_1A_clusters, clusters_with_1_A

implicit none

contains

subroutine sources_and_constants(source,isconst,fitted,ipar)
	implicit none

	real(kind(1.d0)) :: source(neq)
	logical :: isconst(neq)
	integer :: fitted(nclust,0:nclust), ipar(4), n_monomers(n_mol_types), n
	integer :: cluster_numbers(n_1A_clusters)

	source = 0.d0											! Initialize all source terms to 0
	isconst = .false.										! Initialize all concentrations to
	fitted = 0												! change freely according to the eqs.

  !---Setting all monomer concentrations/sources----------------------------------------------------
	! Get the cluster numbers corresponding to the monomers
	call monomer_indices(n_monomers)

!   Use constant concentrations for all monomers
	isconst(n_monomers) = .true.							! All monomer concentrations are constant

!   ...or fit concentrations (here for 1A)...
	isconst(n1A) = .false.
	fitted(1,0) = 1											! 1 concentration is fitted
	call clusters_with_1_A(cluster_numbers)
	n = n_1A_clusters
	fitted(2,0:n) = (/n1A, n-1, cluster_numbers(2:n)/)		! Concentration of 1A is fitted using
															!  n-1 other concentrations: [1A]+... = const.
	! Other monomer concentrations will remain constant

!   (Everything related to fitting or setting something constant commented out
!   ---> all concentrations are allowed to vary according to the equations)

!---Setting sources of generic ions----------------------------------------------------------------
!	source((/nneg,npos/)) = 3.d6


end subroutine sources_and_constants


end module monomer_settingsD
