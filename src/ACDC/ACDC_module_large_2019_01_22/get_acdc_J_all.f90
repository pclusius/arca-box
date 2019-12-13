MODULE ACDC_LARGE
IMPLICIT NONE
private
public :: get_acdc_J
CONTAINS

subroutine get_acdc_J(c_acid,c_base,c_org,cs_ref,temp,ipr,t_sim,solve_ss,j_acdc,diameter_acdc,c_bins, pout)
use acdc_size_param, only: get_system_size, get_vapor_indices, group_size_bins, nbins	! parameters related to the simulation system size
use acdc_system, only : small_set_mode, n_charges, nclust_max! simulation mode and numbers of charging states
use driver_acdc_J, only : acdc_driver

	implicit none


	! Input: ambient conditions and simulation time
	real(kind(1.d0)), intent(inout) :: c_acid, c_base, c_org	! vapor concentrations (1/m^3)
	real(kind(1.d0)), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(kind(1.d0)), intent(in) :: temp				! temperature (K)
	real(kind(1.d0)), intent(in) :: ipr					! ion production rate (1/s/m^3)
	real(kind(1.d0)), intent(in) :: t_sim				! simulation time (s)
	logical, intent(in) :: solve_ss, pout						! solve the steady state or run only for the given time
	! Output: formed particles/time and their acid content
	real(kind(1.d0)), intent(out) :: j_acdc				! simulated formation rate (1/s/m^3)
	real(kind(1.d0)), intent(out) :: diameter_acdc		! mass diameter of the formed particles (m)

	real(kind(1.d0)), save, allocatable :: c(:)			! cluster concentrations
	integer, save :: neq_syst = 0, nclust_syst = 0, nout_syst(n_charges) = 0	! parameters for the simulation system size
	real(kind(1.d0)), save :: diameter_max_syst
	integer, save :: n1acid = 0, n1base = 0, n1org = 0	! cluster numbers for vapor monomers
	real(kind(1.d0)) :: c_bins(nbins)					! size bin concentrations
	integer, save :: grouped_indices(0:nbins,nclust_max)! indices of the clusters belonging to each bin
	integer, save :: nclust_per_bin(0:nbins)			! total number of clusters in each bin

	real(kind(1.d0)) :: j_out(n_charges)				! formation rate vector (for neu only, or for neu, neg and pos)
	real(kind(1.d0)) :: t_in							! input time for the driver (s)
	logical, save :: use_solver							! use the solver or the simple Euler method
	logical :: int_ok = .true.							! .false. if integration fails
	real(kind(1.d0)), save :: t_iter = 1.d-16			! iteration time step (s) for the Euler method
	integer, save :: ipar(4)							! parameters for re-calling the monomer settings and rate constants
	logical, save :: firstcall = .true.
	integer :: i


	if (firstcall) then

		firstcall = .false.

		ipar = 0

		! Determine the system size and vapor indices
		call get_system_size(neq_syst,nclust_syst,nout_syst,diameter_max_syst)
		call get_vapor_indices(n1acid,n1base,n1org)

		! Initialize the concentrations
		allocate(c(neq_syst))
		c = 0.d0

		if (.not. small_set_mode) then
			call group_size_bins(grouped_indices,nclust_per_bin)
		end if
		c_bins = 0.d0

		write(*,*) 'Solving a set of ', neq_syst, ' equations'

		!if (neq_syst .le. 2000) then
			use_solver = .true.
		!else
		!	use_solver = .false.
		!	write(*,*) 'Large system: using Eulerian integration instead of VODE'
		!end if

	end if

	! Initialize the rate constants etc. at every call because of the varying ambient conditions
	ipar(1) = 0
	ipar(3) = 0

	if (solve_ss) then
		! Override the input time with a maximum time that
		! the driver is allowed to try to get to the steady state
		t_in = 5.d6
	else
		t_in = t_sim
	end if

	! Set the vapor concentrations
	c(n1acid) = c_acid
	c(n1base) = c_base
	if (n1org .ne. 0) then
		c(n1org) = c_org
	end if

	! Run ACDC to the steady state / for the given time
	call acdc_driver(neq_syst,nclust_syst,nout_syst,c,cs_ref,temp,ipr,t_in,t_iter,solve_ss,use_solver,ipar,int_ok,j_out)

	j_acdc = sum(j_out)
	!write(*,*) 'Neutral, neg. and pos. J: ',j_out*1.d-6,' cm^-3 s^-1'
	if (j_acdc .lt. 0.d0) then
		if (j_acdc .lt. -1.d-12 .and. int_ok) then
			write(*,*) 'Warning: negative J = ',j_acdc*1.d-6,' cm^-3 s^-1, something wrong?'
		end if
		j_acdc = 0.d0
	end if

	! Update the vapor concentrations, as clusters may be assumed to act as a sink for vapors
	c_acid = c(n1acid)
	c_base = c(n1base)
	if (n1org .ne. 0) then
		c_org = c(n1org)
	end if

	! The first outgrown size is ~equal to the largest simulated size
	diameter_acdc = diameter_max_syst

	if (.not. small_set_mode) then
		! concentration in bins
		do i = 1,nbins
			if (nclust_per_bin(i) .ne. 0) then
				c_bins(i) = sum(c(grouped_indices(i,1:nclust_per_bin(i))))
			end if
      if (pout) print '(a,i0,a,t9,es12.3,a)', 'bin ',i, ': ',c_bins(i), ' /m3'
		end do
	end if


end subroutine get_acdc_J

END MODULE ACDC_LARGE
