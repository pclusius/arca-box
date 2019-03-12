module acdc_system

implicit none

integer, parameter :: nclust = 15						! number of clusters, molecules and ions
integer, parameter :: neq = 22							! number of equations

real(kind(1.d0)), parameter :: temp = 298.15d0				! temperature in K
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+000			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-003

integer, parameter :: n1A = 1, n1D = 3			! cluster indices for monomers and ions
integer, parameter :: nout_neu = 21

integer, parameter :: n_mol_types = 2
integer, parameter :: nmolA = 1, nmolD = 2			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 3				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 15			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! number of negative molecules and clusters
integer, parameter :: n_positives = 0			! number of positive molecules and clusters

integer, parameter :: n_neutral_clusters = 13			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! number of negative clusters
integer, parameter :: n_positive_clusters = 0			! number of positive clusters

real(kind(1.d0)), parameter :: mass_max = 572.64
real(kind(1.d0)), parameter :: diameter_max = 1.15
real(kind(1.d0)), parameter :: mob_diameter_max = 1.45

integer, parameter :: nmols_out_neutral(2, 2) = reshape((/5, 4, 4, 5/),(/2, 2/))			! criteria for outgrowing neutrals


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(15)

	n_A = (/1, 2, 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 3, 4/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(3)

	cluster_numbers = (/1, 4, 7/)

end subroutine clusters_with_1_A

subroutine cluster_arrays(neutral_clusters)
	implicit none
	integer :: neutral_clusters(13)

	neutral_clusters = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/)

end subroutine cluster_arrays

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(15)

	mass = (/98.08, 196.16, 45.08, 143.16, 241.24, 339.32, 188.24, 286.32, 384.40, 482.48, &
		&331.40, 429.48, 527.56, 474.56, 572.64/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(15)

	 diameter = (/0.55, 0.70, 0.59, 0.72, 0.82, 0.90, 0.84, 0.91, 0.98, 1.03, &
		&0.99, 1.04, 1.09, 1.11, 1.15/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(15)

	 mob_diameter = (/0.85, 1.00, 0.89, 1.02, 1.12, 1.20, 1.14, 1.21, 1.28, 1.33, &
		&1.29, 1.34, 1.39, 1.41, 1.45/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(22) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '1D'
	clust(4)(:) = '1A1D'
	clust(5)(:) = '2A1D'
	clust(6)(:) = '3A1D'
	clust(7)(:) = '1A2D'
	clust(8)(:) = '2A2D'
	clust(9)(:) = '3A2D'
	clust(10)(:) = '4A2D'
	clust(11)(:) = '2A3D'
	clust(12)(:) = '3A3D'
	clust(13)(:) = '4A3D'
	clust(14)(:) = '3A4D'
	clust(15)(:) = '4A4D'
	clust(16)(:) = 'source'
	clust(17)(:) = 'coag'
	clust(18)(:) = 'wall'
	clust(19)(:) = 'dil'
	clust(20)(:) = 'rec'
	clust(21)(:) = 'out_neu'
	clust(22)(:) = 'bound'

end subroutine cluster_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(2) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'D'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(2)

	n_monomers = (/1, 3/)

end subroutine monomer_indices


end module acdc_system

