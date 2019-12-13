module acdc_system

implicit none

integer, parameter :: nclust_max = 21041						! max. number of clusters, molecules and ions
integer, parameter :: neq_max = 21042							! max. number of equations
integer, parameter :: nclust = -1, neq = -1						! not pre-determined in the restrict mode
logical, parameter :: restrict_system = .true.							! exclude given cluster compositions
logical, parameter :: small_set_mode = .false.							! very large set of clusters (loop mode)

real(kind(1.d0)), parameter :: temp = 298.15d0				! temperature in K
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+00			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-03

integer, parameter :: n_monomer_types = 3
integer, parameter :: n_charges = 1			! number of charging states
integer, parameter :: nout_all(1) = -1			! not pre-determined

integer, parameter :: n_mol_types = 3
integer, parameter :: nmolA = 1, nmolN = 2, nmolO = 3			! molecule indices for the used species

integer, parameter :: ij_ind_max(3) = (/6,5,500/)
real(kind(1.d0)), parameter :: diameter_max = 6.99

contains

subroutine monomer_names(clust_mon)
	implicit none
	character(len=11), dimension(3) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1O'

end subroutine monomer_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(3) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'N'
	labels(3)(:) = 'O'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(3)
	integer :: clust_from_indices(0:12,0:10,0:1000), i, ij_ind(3)

	call get_cluster_numbers(clust_from_indices)
	n_monomers = 0

	do i = 1,3
		ij_ind = 0
		ij_ind(i) = 1
		n_monomers(i) = clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3))
	end do

end subroutine monomer_indices


end module acdc_system

