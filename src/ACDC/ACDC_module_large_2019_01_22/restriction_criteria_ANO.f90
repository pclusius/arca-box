module restriction_criteria

implicit none

real(kind(1.d0)), parameter :: diameter_max_loop = 3d0*1.d-9	! size of the extended system (loop mode) (m)

contains

subroutine check_validity(indices,valid)
use acdc_system, only : small_set_mode, ij_ind_max, n_mol_types, molecule_names
use acdc_hsdata, only : get_hsdata

	implicit none
	integer :: indices(n_mol_types), boundary_indices(n_mol_types), imol, ind(2)
	integer :: valid
	integer, save :: nmolacid=0, nmolbase=0, nmolorg=0
	logical, save :: first_call = .true.
	real(kind(1.d0)), save :: r_max
	real(kind(1.d0)) :: xmol, xmols(n_mol_types), r_min, r, m, deltah, deltas
	character(len=11), dimension(n_mol_types) :: mol_names
	logical :: lfound, first_mol, lhs

	if (.not. small_set_mode) then ! the conditional needs to be included to compile also in the small set mode

		if (first_call) then

			first_call = .false.

			call molecule_names(mol_names)

			do imol = 1,n_mol_types
				if (trim(mol_names(imol)(:)) .eq. 'A') then
					nmolacid = imol
				else if ((trim(mol_names(imol)(:)) .eq. 'N') .or. (trim(mol_names(imol)(:)) .eq. 'D')) then
					nmolbase = imol
				else if (trim(mol_names(imol)(:)) .eq. 'O') then
					nmolorg = imol
				end if
			end do
			! The system is restricted to sizes <= r_max
			call get_masses_and_radii(m,r,ij_ind_max) ! Assuming AN-seed and growth by O
			r_max = r
			if (diameter_max_loop .gt. 2.d0*r_max) then
				write(*,*) 'Using the max. possible diameter of ',2.d0*r_max*1.d9,' nm instead of ',diameter_max_loop*1.d9,' nm set by the user'
			else
				r_max = 0.5d0*diameter_max_loop
				write(*,*) 'Using a max. diameter of ',2.d0*r_max*1.d9,' nm'
			end if

		end if

		! See if the input composition fulfills the criteria
		valid = 1
		if ((sum(indices) .eq. 0) .or. (any(indices .gt. ij_ind_max))) then
			valid = -1
			return
		end if

		! First, see if the composition is stable
		if (sum(indices) .eq. 1) then
			valid = 1
			return
		end if

! ------ The following part sets the allowed chemical compositions ------
		! Regime for which there is QC data
		!if (all(indices .le. 5)) then ! Allow scavenging of these by larger ACDC-particles (often negligible)
		if ((indices(nmolacid) .le. 5) .and. (indices(nmolbase) .le. 5)) then
			if (sum(indices) .gt. indices(nmolacid)+indices(nmolbase)) then
				! The smallest clusters are assumed to contain only acid and base
				valid = -1
				return
			end if
			ind(1) = indices(nmolacid)
			ind(2) = indices(nmolbase)
			call get_hsdata(ind,lhs,deltah,deltas)
			if (.not. lhs) then
				valid = -1
			end if
			return
		end if
		! Regime modeled by macroscopic thermodynamics
		do imol = 1,n_mol_types
			xmols(imol) = real(indices(imol),kind=kind(1.d0))/real(sum(indices),kind=kind(1.d0))
		end do
		if ((xmols(nmolbase) .gt. 1.2d0*xmols(nmolacid)) .or.&
				&(xmols(nmolbase) .lt. 0.8d0*xmols(nmolacid))) then
			valid = -1
			return
		end if
! -----------------------------------------------------------------------

		! If the composition is stable, see if the size is within the given range
		call get_masses_and_radii(m,r,indices)
		if (r .gt. r_max) then
			valid = 0
			return
		end if

		if ((valid .eq. 1) .and. any(indices .gt. ij_ind_max)) then
			write(*,*) 'Error: Composition (',indices,') classified as valid?'
			stop
		end if

	end if

end subroutine check_validity


end module restriction_criteria
