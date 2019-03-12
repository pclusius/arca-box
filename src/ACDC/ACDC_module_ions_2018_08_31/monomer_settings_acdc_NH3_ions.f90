module monomer_settings

use acdc_system, only : neq, nclust
use acdc_system, only : monomer_arrays, n_neutral_monomers, n_negative_monomers, n_positive_monomers
use acdc_system, only : n1A, n_1A_clusters, clusters_with_1_A

implicit none

contains

subroutine sources_and_constants(source,isconst,fitted,ipar)
	implicit none

	real(kind(1.d0)) :: source(neq)
	logical :: isconst(neq)
	integer :: fitted(nclust,0:nclust), ipar(4), n
	integer :: cluster_numbers(n_1A_clusters)
	integer :: neutral_monomers(n_neutral_monomers), negative_monomers(n_negative_monomers), positive_monomers(n_positive_monomers)
	
	source = 0.d0											! Initialize all source terms to 0
	isconst = .false.										! Initialize all concentrations to 
	fitted = 0												! change freely according to the eqs.

!---Setting all monomer concentrations/sources----------------------------------------------------

	! Get the cluster numbers corresponding to the monomers
	call monomer_arrays(neutral_monomers, negative_monomers, positive_monomers)

!   Use constant concentrations for monomers	
	isconst(neutral_monomers) = .true.						! Neutral monomer concentrations are constant
															! (Ionic monomers can vary according to the ion production rate)
	
!   ...or fit concentrations (here for 1A)...
!	isconst(n1A) = .false.
!	fitted(1,0) = 1											! 1 concentration is fitted
!	call clusters_with_1_A(cluster_numbers)
!	n = n_1A_clusters
!	fitted(2,0:n) = (/n1A, n-1, cluster_numbers(2:n)/)		! Concentration of 1A is fitted using
															!  n-1 other concentrations: [1A]+... = const.
	! Other monomer concentrations will remain as set above

!   (Everything related to fitting or setting something constant commented out
!   ---> all concentrations are allowed to vary according to the equations)

!---Setting sources of generic ions----------------------------------------------------------------
!	source((/nneg,npos/)) = 3.d6

	
end subroutine sources_and_constants


end module monomer_settings
