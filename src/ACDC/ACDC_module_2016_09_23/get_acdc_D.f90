subroutine get_acdc_D(c_acid,c_base,c_org,cs_ref,temp,t_sim,solve_ss,j_acdc,diameter_acdc)
use acdc_systemD, only : nclust, neq					! number of clusters and equations
use acdc_systemD, only : cluster_names				! names of the clusters and fluxes
use acdc_systemD, only : n1A							! cluster numbers of acid and base monomers
!use acdc_system, only : n1base => n1N
!use acdc_system, only : n1B, n1N1P
!use acdc_system, only : nneg, npos					! cluster numbers of charger ions
use acdc_systemD, only : nout_neu					! numbers for formation fluxes
!use acdc_system, only : nout_neg, nout_pos
use acdc_systemD, only : nmols_out_neutral, nmolA	! criteria for "nucleation"
use acdc_systemD, only : diameter_max				! max. mass diameter (nm)
use driverD, only : acdc_driver						! driver

	implicit none
	INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=300)
	! Input: ambient conditions and simulation time
	real(real_x), intent(inout) :: c_acid, c_base, c_org	! vapor concentrations (1/m^3)
	real(real_x), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(real_x), intent(in) :: temp				! temperature (K)
!	real(real_x), intent(in) :: q_ion				! ion source rate (1/s/m^3)
	real(real_x), intent(in) :: t_sim				! simulation time (s)
	logical, intent(in) :: solve_ss						! solve the steady state or run only for the given time
	! Output: formed particles/time and their acid content
	real(real_x), intent(out) :: j_acdc				! simulated formation rate (1/s/m^3)
	real(real_x), intent(out) :: diameter_acdc		! mass diameter of the formed particles (m)
	!real(real_x), intent(out) :: nacid_acdc		! number of acid molecules at the "formation size"
														! NB: nacid_acdc is type real (dp), as it is like this in UHMA

	real(real_x), save :: c(neq)					! cluster concentrations
	character(len=11), dimension(neq) :: c_names		! cluster and flux names
	integer, save :: n1base = 0, n1org = 0				! cluster numbers of base and organic molecules
	integer, save :: nacid_out							! smallest number of acid molecules in the outgrown clusters
	real(real_x) :: j_out(1)						! formation rate vector (for neu, neg and pos)
	real(real_x) :: t_in							! input time for the driver (s)
	logical :: use_solver = .true.						! use the solver or the simple Euler method
	logical :: int_ok = .true.							! .false. if integration fails
	real(real_x), save :: t_iter = 1.d-10			! iteration time step (s) for the Euler method
	integer, save :: ipar(4)							! parameters for re-calling the monomer settings and rate constants
	logical, save :: firstcall = .true.
	integer :: n


	! Initialize the rate constants etc. at every call
	! because of the varying ambient conditions
	ipar(1:3) = 0

	if (firstcall) then
		firstcall = .false.
		! Initialize the concentrations
		c = 0.d0
		! See which compounds there are
		call cluster_names(c_names)
		do n = 1, nclust
			if ((trim(c_names(n)(:)).eq.'1N') .or. (trim(c_names(n)(:)).eq.'1D')) then
				if (n1base .ne. 0) then
					write (*,*) 'Found more than one base molecule'
					stop
				end if
				n1base = n
			else if (trim(c_names(n)(:)) .eq. '1O') then
				n1org = n
			end if
		end do
		! Get the acid content of the outgrown clusters
		! (use the smallest number of acids)
		!nacid_out = minval(nmols_out_neutral(1:size(nmols_out_neutral,1),nmolA))
	end if

	if (solve_ss) then
		! Override the input time with a maximum time that
		! the driver is allowed to try to get to the steady state
		t_in = 1.d6
	else
		t_in = t_sim
	end if

	! Set the vapor concentrations
	c(n1A) = c_acid
	c(n1base) = c_base
	if (n1org .ne. 0) then
		c(n1org) = c_org
	end if

	! Run ACDC to the steady state / for the given time
	call acdc_driver(c,cs_ref,temp,t_in,t_iter,solve_ss,use_solver,ipar,int_ok,j_out)
	j_acdc = j_out(1)
	if (j_acdc .lt. 0.d0) then
		if (j_acdc .lt. -1.d-12 .and. int_ok) then
			write(*,*) 'Warning: negative J = ',j_acdc*1.d-6,' cm^-3 s^-1, something wrong?'
		end if
		j_acdc = 0.d0
	end if

	! The first outgrown size is ~equal to the largest simulated size
	diameter_acdc = diameter_max*1.d-9

	!nacid_acdc = real(nacid_out,kind=real_x)


end subroutine get_acdc_D
