module acdc_systemD

implicit none

integer, parameter :: nclust = 24						! number of clusters, molecules and ions
integer, parameter :: neq = 31							! number of equations

real(kind(1.d0)), parameter :: temp = 298.15d0				! temperature in K
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+00			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-03

integer, parameter :: n1A = 1, n1D = 5			! cluster indices for monomers and ions
integer, parameter :: nout_neu = 30

integer, parameter :: n_mol_types = 2
integer, parameter :: nmolA = 1, nmolD = 2			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 5				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! number of negative molecules and clusters
integer, parameter :: n_positives = 0			! number of positive molecules and clusters

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! number of negative clusters
integer, parameter :: n_positive_clusters = 0			! number of positive clusters

real(kind(1.d0)), parameter :: mass_max = 572.64
real(kind(1.d0)), parameter :: diameter_max = 1.15
real(kind(1.d0)), parameter :: mob_diameter_max = 1.45

integer, parameter :: nmols_out_neutral(1, 2) = reshape((/5, 4/),(/1, 2/))			! criteria for outgrowing neutrals


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(24)

	n_A = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(5)

	cluster_numbers = (/1, 6, 11, 16, 21/)

end subroutine clusters_with_1_A

subroutine cluster_arrays(neutral_clusters)
	implicit none
	integer :: neutral_clusters(22)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)

end subroutine cluster_arrays

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(24)

	mass = (/98.08, 196.16, 294.24, 392.32, 45.08, 143.16, 241.24, 339.32, 437.40, 90.16, &
		&188.24, 286.32, 384.40, 482.48, 135.24, 233.32, 331.40, 429.48, 527.56, 180.32, &
		&278.40, 376.48, 474.56, 572.64/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(24)

	 diameter = (/0.55, 0.70, 0.80, 0.88, 0.59, 0.72, 0.82, 0.90, 0.96, 0.75, &
		&0.84, 0.91, 0.98, 1.03, 0.86, 0.93, 0.99, 1.04, 1.09, 0.94, &
		&1.00, 1.06, 1.11, 1.15/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(24)

	 mob_diameter = (/0.85, 1.00, 1.10, 1.18, 0.89, 1.02, 1.12, 1.20, 1.26, 1.05, &
		&1.14, 1.21, 1.28, 1.33, 1.16, 1.23, 1.29, 1.34, 1.39, 1.24, &
		&1.30, 1.36, 1.41, 1.45/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(31) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '3A'
	clust(4)(:) = '4A'
	clust(5)(:) = '1D'
	clust(6)(:) = '1A1D'
	clust(7)(:) = '2A1D'
	clust(8)(:) = '3A1D'
	clust(9)(:) = '4A1D'
	clust(10)(:) = '2D'
	clust(11)(:) = '1A2D'
	clust(12)(:) = '2A2D'
	clust(13)(:) = '3A2D'
	clust(14)(:) = '4A2D'
	clust(15)(:) = '3D'
	clust(16)(:) = '1A3D'
	clust(17)(:) = '2A3D'
	clust(18)(:) = '3A3D'
	clust(19)(:) = '4A3D'
	clust(20)(:) = '4D'
	clust(21)(:) = '1A4D'
	clust(22)(:) = '2A4D'
	clust(23)(:) = '3A4D'
	clust(24)(:) = '4A4D'
	clust(25)(:) = 'source'
	clust(26)(:) = 'coag'
	clust(27)(:) = 'wall'
	clust(28)(:) = 'dil'
	clust(29)(:) = 'rec'
	clust(30)(:) = 'out_neu'
	clust(31)(:) = 'bound'

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

	n_monomers = (/1, 5/)

end subroutine monomer_indices


end module acdc_systemD

