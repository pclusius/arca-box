module acdc_size_param
use acdc_system
! use DATA_FORMAT, only: n_bins_acdc
implicit none

integer, parameter :: nbins = 17					! number of size bins for size-classifying the clusters

contains

subroutine get_system_size(neq_syst,nclust_syst,nout_syst,diameter_max_syst)
use restriction_criteria, only : diameter_max_loop

	implicit none
	integer :: nclust_syst, neq_syst			! numbers of clusters and equations
	integer :: nout_syst(n_charges)				! indices of outgoing fluxes
	real(kind(1.d0)) :: diameter_max_syst		! max. diameter in system
	integer, allocatable :: clust_comp(:,:)		! cluster composition in the loop mode
   
	if (small_set_mode) then

		neq_syst = neq
		nclust_syst = nclust

		nout_syst = nout_all

		diameter_max_syst = diameter_max*1.d-9

	else

		allocate(clust_comp(nclust_max,n_mol_types))
		if (n_mol_types .gt. 1) then
			call get_molecule_numbers(nclust_syst,clust_comp)
		else
			nclust_syst = nclust_max
		end if
		neq_syst = nclust_syst+(neq_max-nclust_max)

		nout_syst = nclust_syst+1				! loop mode: last, additional index

		diameter_max_syst = diameter_max_loop

	end if

end subroutine get_system_size


subroutine get_vapor_indices(n1acid,n1base,n1org)
	implicit none
	integer :: n1acid, n1base, n1org
	integer :: n1vapor(4), n_monomers(n_monomer_types)
	character(len=3), dimension(4) :: names_vapor_def
	character(len=11), dimension(n_monomer_types) :: names_monomers
	integer :: i, j

	! Find the cluster numbers for acid, base and possible organic compound
	n1acid = 0
	n1base = 0
	n1org = 0

	n1vapor = 0

	names_vapor_def(1)(:) = '1A'
	names_vapor_def(2)(:) = '1N'
	names_vapor_def(3)(:) = '1D'
	names_vapor_def(4)(:) = '1O'

	call monomer_names(names_monomers)
	call monomer_indices(n_monomers)

	do i = 1,size(names_vapor_def)
		do j = 1,size(names_monomers)
			if (trim(names_vapor_def(i)(:)) .eq. trim(names_monomers(j)(:))) then
				n1vapor(i) = n_monomers(j)
				exit
			end if
		end do
	end do

	n1acid = n1vapor(1)
	n1base = max(n1vapor(2),n1vapor(3))
	n1org = n1vapor(4)

	if (min(n1vapor(2),n1vapor(3)) .ne. 0) then
		write (*,*) 'Found more than one base molecule?'
		stop
	end if

end subroutine get_vapor_indices


subroutine get_bin_limits(bin_limits)
	implicit none
	real(kind(1.d0)) :: bin_limits(nbins+1)   ! vector of the bin limits (m)

	! Bin limits according to mobility diameter
	bin_limits = (/0.951d0, 1.05d0, 1.16d0, 1.28d0, 1.42d0, 1.57d0, 1.73d0, 1.92d0, &
                  2.12d0, 2.34d0, 2.59d0, 2.86d0, 3.16d0, 3.49d0, 3.86d0, 4.27d0, &
                  4.71d0, 5.21d0 /)*1.d-9 ! m
	! bin_limits = (/1.05d0, 1.28d0, 1.73d0, 2.59d0, 4.27d0, 6.36d0/)*1.d-9 ! m

end subroutine get_bin_limits


subroutine group_size_bins(grouped_indices,nclust_per_bin)
	implicit none
	integer :: grouped_indices(0:nbins,nclust_max)			! indices of the clusters belonging to each bin
	integer :: nclust_per_bin(0:nbins)						! total number of clusters belonging to each bin
	real(kind(1.d0)) :: mi, ri, ivalue
	real(kind(1.d0)) :: bin_limits(nbins+1)					! vector of the bin limits (m)
	integer :: i, j, nclust_syst, indices(nclust_max,n_mol_types)
	logical :: lfound

	grouped_indices = 0
	nclust_per_bin =  0
	lfound = .false.

	if(n_mol_types .gt. 1) then
		call get_molecule_numbers(nclust_syst,indices)
	else
		nclust_syst = nclust_max
		indices = reshape((/((i,i=1,nclust_syst),j=1,n_mol_types)/), shape(indices))
	end if

	call get_bin_limits(bin_limits)
	do i=1,nclust_syst
		if (sum(indices(i,:)) .eq. 1) then
			cycle ! Don't add monomers to the size bins
		end if
		call get_masses_and_radii(mi,ri,indices(i,:))
		ivalue = 2.d0*ri + 0.3d0*1.d-9	! assuming d_mob = d_mass + 0.3 nm
		if (ivalue .ge. bin_limits(1)) then
			lfound = .false.
			do j=1,nbins
				if ((ivalue .ge. bin_limits(j)) .and. (ivalue .lt. bin_limits(j+1))) then
					nclust_per_bin(j) = nclust_per_bin(j) + 1
					grouped_indices(j,nclust_per_bin(j)) = i
					lfound = .true.
					exit
				end if
			end do
			if (.not. lfound) then
				write(*,*) 'Could not group cluster ', i, ', molecules: ', indices(i,:)
				stop
			end if
		else
			nclust_per_bin(0) = nclust_per_bin(0) + 1
			grouped_indices(0,nclust_per_bin(0)) = i
			!write(*,*) 'Cluster below the lowest bin limit: ' , i, ', molecules: ', indices(i,:)
		end if
	end do

end subroutine group_size_bins


end module acdc_size_param
