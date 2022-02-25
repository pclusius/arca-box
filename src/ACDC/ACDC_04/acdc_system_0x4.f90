! Input files used to create this file:
! System name: DLPNO NH3 system
! System file: Perl_input/Besel_H2SO4_NH3_TO/INPUT.inp
! Energy file: Perl_input/Besel_H2SO4_NH3_TO/HS_YesSymYesQuasiHCorr.txt
! Dipole file: Perl_input/Besel_H2SO4_NH3_TO/dipoles_YesQuasiHCorr.txt
module acdc_system_0x4

implicit none

integer, parameter :: nclust = 65						! number of clusters, molecules and ions
integer, parameter :: neq = 74							! number of equations
integer, parameter :: nclust_max = 65, neq_max = 74		! same as nclust and neq in the small set mode
logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: n_variables = 4
logical, parameter :: variable_temp = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+00			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-03

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1N = 4, n1B = 25, n1P1N = 47, nneg = 64, npos = 65			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/71, 72, 73/), nout_neu = 71, nout_neg = 72, nout_pos = 73			! indices for outgoing fluxes

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolP = 3, nmolN = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 3				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 23			! negative
integer, parameter :: n_positives = 18			! positive

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 4/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/25/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/47/)			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 21			! negative
integer, parameter :: n_positive_clusters = 16			! positive

real(kind(1.d0)), parameter :: mass_max = 691.72
real(kind(1.d0)), parameter :: diameter_max = 1.14
real(kind(1.d0)), parameter :: mob_diameter_max = 1.44
integer, parameter :: ij_ind_max(4) = (/6, 1, 1, 6/)		! maximum molecular content
integer, parameter :: n_bound = 11		! number of clusters at the system boundary

integer, parameter :: nmols_out_neutral(1, 4) = reshape((/7, 0, 0, 6/),(/1, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/6, 1, 0, 3/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/6, 0, 1, 7/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(65)

	n_A = (/1, 2, 3, 0, 1, 2, 3, 1, 2, 3, &
		&4, 2, 3, 4, 5, 3, 4, 5, 6, 4, &
		&5, 6, 5, 6, 1, 2, 3, 4, 1, 2, &
		&3, 4, 2, 3, 4, 5, 3, 4, 5, 6, &
		&4, 5, 6, 5, 6, 6, 0, 1, 0, 1, &
		&2, 1, 2, 3, 2, 3, 4, 3, 4, 5, &
		&4, 5, 6, 0, 0/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(3)

	cluster_numbers = (/1, 5, 8/)

end subroutine clusters_with_1_A

subroutine arrays(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(24), negatives(23), positives(18)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negatives = (/25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46, 64/)
	positives = (/47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63, 65/)

end subroutine arrays

subroutine cluster_arrays(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(22), negative_clusters(21), positive_clusters(16)

	neutral_clusters = (/2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negative_clusters = (/26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46/)
	positive_clusters = (/48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63/)

end subroutine cluster_arrays

subroutine get_charging_state(charging_state)
	implicit none
	integer :: charging_state(65)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 2, 3/)

end subroutine get_charging_state

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(65)

	mass = (/98.08, 196.16, 294.24, 17.04, 115.12, 213.20, 311.28, 132.16, 230.24, 328.32, &
		&426.40, 247.28, 345.36, 443.44, 541.52, 362.40, 460.48, 558.56, 656.64, 477.52, &
		&575.60, 673.68, 592.64, 690.72, 97.08, 195.16, 293.24, 391.32, 114.12, 212.20, &
		&310.28, 408.36, 229.24, 327.32, 425.40, 523.48, 344.36, 442.44, 540.52, 638.60, &
		&459.48, 557.56, 655.64, 574.60, 672.68, 689.72, 18.04, 116.12, 35.08, 133.16, &
		&231.24, 150.20, 248.28, 346.36, 265.32, 363.40, 461.48, 380.44, 478.52, 576.60, &
		&495.56, 593.64, 691.72, 32.00, 19.02/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(65)

	 diameter = (/0.55, 0.70, 0.80, 0.43, 0.63, 0.75, 0.84, 0.69, 0.79, 0.87, &
		&0.94, 0.83, 0.91, 0.97, 1.03, 0.94, 1.00, 1.05, 1.10, 1.02, &
		&1.07, 1.12, 1.10, 1.14, 0.55, 0.70, 0.80, 0.88, 0.63, 0.75, &
		&0.84, 0.91, 0.79, 0.87, 0.94, 1.00, 0.90, 0.97, 1.03, 1.08, &
		&1.00, 1.05, 1.10, 1.07, 1.12, 1.14, 0.43, 0.63, 0.54, 0.69, &
		&0.79, 0.74, 0.83, 0.91, 0.87, 0.94, 1.00, 0.96, 1.02, 1.07, &
		&1.05, 1.10, 1.14, 0.45, 0.39/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(65)

	 mob_diameter = (/0.85, 1.00, 1.10, 0.73, 0.93, 1.05, 1.14, 0.99, 1.09, 1.17, &
		&1.24, 1.13, 1.21, 1.27, 1.33, 1.24, 1.30, 1.35, 1.40, 1.32, &
		&1.37, 1.42, 1.40, 1.44, 0.85, 1.00, 1.10, 1.18, 0.93, 1.05, &
		&1.14, 1.21, 1.09, 1.17, 1.24, 1.30, 1.20, 1.27, 1.33, 1.38, &
		&1.30, 1.35, 1.40, 1.37, 1.42, 1.44, 0.73, 0.93, 0.84, 0.99, &
		&1.09, 1.04, 1.13, 1.21, 1.17, 1.24, 1.30, 1.26, 1.32, 1.37, &
		&1.35, 1.40, 1.44, 0.75, 0.69/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(74) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '3A'
	clust(4)(:) = '1N'
	clust(5)(:) = '1A1N'
	clust(6)(:) = '2A1N'
	clust(7)(:) = '3A1N'
	clust(8)(:) = '1A2N'
	clust(9)(:) = '2A2N'
	clust(10)(:) = '3A2N'
	clust(11)(:) = '4A2N'
	clust(12)(:) = '2A3N'
	clust(13)(:) = '3A3N'
	clust(14)(:) = '4A3N'
	clust(15)(:) = '5A3N'
	clust(16)(:) = '3A4N'
	clust(17)(:) = '4A4N'
	clust(18)(:) = '5A4N'
	clust(19)(:) = '6A4N'
	clust(20)(:) = '4A5N'
	clust(21)(:) = '5A5N'
	clust(22)(:) = '6A5N'
	clust(23)(:) = '5A6N'
	clust(24)(:) = '6A6N'
	clust(25)(:) = '1B'
	clust(26)(:) = '1A1B'
	clust(27)(:) = '2A1B'
	clust(28)(:) = '3A1B'
	clust(29)(:) = '1B1N'
	clust(30)(:) = '1A1B1N'
	clust(31)(:) = '2A1B1N'
	clust(32)(:) = '3A1B1N'
	clust(33)(:) = '1A1B2N'
	clust(34)(:) = '2A1B2N'
	clust(35)(:) = '3A1B2N'
	clust(36)(:) = '4A1B2N'
	clust(37)(:) = '2A1B3N'
	clust(38)(:) = '3A1B3N'
	clust(39)(:) = '4A1B3N'
	clust(40)(:) = '5A1B3N'
	clust(41)(:) = '3A1B4N'
	clust(42)(:) = '4A1B4N'
	clust(43)(:) = '5A1B4N'
	clust(44)(:) = '4A1B5N'
	clust(45)(:) = '5A1B5N'
	clust(46)(:) = '5A1B6N'
	clust(47)(:) = '1P1N'
	clust(48)(:) = '1A1P1N'
	clust(49)(:) = '1P2N'
	clust(50)(:) = '1A1P2N'
	clust(51)(:) = '2A1P2N'
	clust(52)(:) = '1A1P3N'
	clust(53)(:) = '2A1P3N'
	clust(54)(:) = '3A1P3N'
	clust(55)(:) = '2A1P4N'
	clust(56)(:) = '3A1P4N'
	clust(57)(:) = '4A1P4N'
	clust(58)(:) = '3A1P5N'
	clust(59)(:) = '4A1P5N'
	clust(60)(:) = '5A1P5N'
	clust(61)(:) = '4A1P6N'
	clust(62)(:) = '5A1P6N'
	clust(63)(:) = '6A1P6N'
	clust(64)(:) = 'neg'
	clust(65)(:) = 'pos'
	clust(66)(:) = 'source'
	clust(67)(:) = 'coag'
	clust(68)(:) = 'wall'
	clust(69)(:) = 'dil'
	clust(70)(:) = 'rec'
	clust(71)(:) = 'out_neu'
	clust(72)(:) = 'out_neg'
	clust(73)(:) = 'out_pos'
	clust(74)(:) = 'bound'

end subroutine cluster_names

subroutine monomer_names(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1P1N'

end subroutine monomer_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'P'
	labels(4)(:) = 'N'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 4, 25, 47/)

end subroutine monomer_indices

subroutine get_bound(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(11), nmols_bound(11,4)

	nmols_bound(1,:) = (/6, 0, 0, 4/)
	nmols_bound(2,:) = (/6, 0, 0, 5/)
	nmols_bound(3,:) = (/5, 0, 0, 6/)
	nmols_bound(4,:) = (/6, 0, 0, 6/)
	nmols_bound(5,:) = (/5, 1, 0, 3/)
	nmols_bound(6,:) = (/5, 1, 0, 4/)
	nmols_bound(7,:) = (/5, 1, 0, 5/)
	nmols_bound(8,:) = (/5, 1, 0, 6/)
	nmols_bound(9,:) = (/4, 0, 1, 6/)
	nmols_bound(10,:) = (/5, 0, 1, 6/)
	nmols_bound(11,:) = (/6, 0, 1, 6/)

	bound_clusters = (/19, 22, 23, 24, 40, 43, 45, 46, 61, 62, 63/)

end subroutine get_bound


end module acdc_system_0x4

