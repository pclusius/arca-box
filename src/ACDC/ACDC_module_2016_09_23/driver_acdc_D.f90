module driverD
use acdc_systemD, only : nclust, neqn => neq				! numbers of clusters and equations
use acdc_systemD, only : nout_neu						! numbers for formation fluxes
!use acdc_system, only : nout_neg, nout_pos
use solution_settings
use monomer_settingsD , only : sources_and_constants		! for determining fitted concentrations

implicit none

contains

subroutine acdc_driver(c,cs_ref,temperature,t_max,t_iter,solve_ss,use_solver,ipar,ok,j_out)

	implicit none

	! Input and output
	! Cluster distribution
	real(kind(1.d0)), intent(inout) :: c(neqn)			! initial concentrations -> final concentrations (1/m^3)
	! Ambient conditions
	real(kind(1.d0)), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(kind(1.d0)), intent(in) :: temperature			! temperature (K)
!	real(kind(1.d0)), intent(in) :: q_ion				! ion source rate (1/s/m^3)
	! Simulation settings and outcome
	real(kind(1.d0)), intent(in) :: t_max				! simulation time (s)
	real(kind(1.d0)), intent(inout) :: t_iter       	! iteration time step for the Euler method (s)
    logical, intent(in) :: solve_ss						! solve the steady state or run only for t_max
	logical, intent(in) :: use_solver					! use the solver or the simple Euler method
	integer, intent(inout) :: ipar(4)					! parameters for re-calling the monomer settings and rate constants
	logical, intent(out) :: ok							! .false. if integration failed
	real(kind(1.d0)), intent(out) :: j_out(1)			! simulated formation rates (neutral, neg and pos) (1/s/m^3)

	! Variables and parameters for the simulation
	real(kind(1.d0)) :: c00(neqn), ci(neqn), ci_temp(neqn), f_out(neqn)	! concentrations and dc/dt
	real(kind(1.d0)) :: c_out							! total concentration of outgrown particles at the beginning
	integer, parameter :: nt = 4, nt_ss = 12			! number of timeranges for integration by the solver
	real(kind(1.d0)) :: t(nt), t_ss(nt_ss)				! start and end times
	real(kind(1.d0)) :: t00, t0, tfin, ti, tt(2)
	real(kind(1.d0)) :: t_iter_tot			            ! the cumulative iterated time (s)
	real(kind(1.d0)) :: t_iter_temp			            ! temporary iterated time for monitoring the convergence to steady state
	logical :: save_ci_temp					        	! logical used in monitoring the convergence
	real(kind(1.d0)) :: parameters(2)					! one or more of the following parameters (max. 4 param.):
														! temperature (1 param.), coagulation sink (1) and ion source rate (2; neg and pos)
	integer :: i, j, k, n, nind, nt_tot

	! Parameters for the solver
	integer, parameter :: itol = 1						! same absolute tolerance for all components
	integer, parameter :: itask = 4						! stop at or beyond t=TOUT, but not beyond t_max, and return
	integer, parameter :: iopt = 1						! allow using optional input in the solver
	integer, parameter :: mf = 22						! full jacobian, computed numerically
	real(kind(1.d0)), parameter :: rtol = rtol_solver	! relative tolerance in the solver
	real(kind(1.d0)), parameter :: atol = atol_solver	! absolute tolerance in the solver
	integer :: lwork, liwork
	integer :: istate, iwork(neqn+30)
	real(kind(1.d0)) :: work(22+9*neqn+2*neqn**2)

	real(kind(1.d0)) :: j_tot, j_by_charge(4), j_by_cluster(neqn), j_all(neqn,4) ! output of the formation rate subroutine
	external fevalD, jevalD, formationD					! subroutines for equations, jacobian and formation rate


	! Initialize
	t0 = 0.d0
	t00 = t0
	call fix_fitted_c(c,c,ipar)							! Fix the concentrations of the fitted species
	c00 = c
	j_out=0.d0

	parameters = (/temperature, cs_ref/)				! Values for the varied parameters
	!parameters = (/temperature, cs_ref, q_ion, q_ion/)

	if (.not. solve_ss) then
		! Save the concentrations of outgrown particles and set them to zero in c,
		! so that only particles grown out in this time interval are used for calculating the formation rate
		c_out = c(nout_neu)
		c(nout_neu) = 0.d0
		!c_out = c((/nout_neu, nout_neg, nout_pos/))
		!c((/nout_neu, nout_neg, nout_pos/)) = 0.d0
	end if

	! Use either the solver or the Eulerian iteration
	if (use_solver) then

		if (solve_ss) then
            t_ss(1:nt_ss) = (/t00+(/1.d-4, 1.d-2, 1.d0, 1.d2, 1.d3, 2.d4, 3.d4, 4.d4, 5.d4, 1.d5, 1.d6, 1.d7/)/)
            nt_tot = nt_ss
        else
            t(1:nt) = (/t00+(/1.d-5, 1.d-3, 1.d-1/), t_max/)
            nt_tot = nt
        end if

		lwork = 22+9*neqn+2*neqn**2
		liwork = neqn+30
		work(5:10) = 0.d0						! use default parameter values in the solver
		work(1) = t_max							! don't go beyond the end time
		iwork(5:10) = 0							! use default parameter values in the solver
		iwork(6) = 1000000						! allow more steps in the solver (probably not needed)
		istate = 1								! tells the solver that this is the first call
												! -> 2 in the solver after succesful integration

		! break the simulation into several time intervals to improve convergence
		do i=1,nt_tot

			ci = c												! save the previous concentrations
			ti = t0

			if (solve_ss) then
				tfin = t_ss(i)
            else
				tfin = t(i)
            end if
			tfin = min(max(tfin,t0*10.d0),t_max)

			ok = .true.
			call DVODE (fevalD,neqn,c,t0,tfin,itol,rtol,atol,&	! calling the solver
				& itask,istate,iopt,work,lwork,iwork,liwork,&
				& jevalD,mf,parameters,ipar)
			if (istate .ne. 2) then								! checking everything went ok - if not, print error
				!write(*,*) 'ERROR: returned istate =', istate
				ok = .false.
			elseif (minval(c) .lt. negtol) then					! checking all concentrations > 0 - if not, print error
				!write(*,*) 'Negative concentrations: c_min = ', minval(c),' at t = ',t0
				ok = .false.
			end if
			if (.not. ok) then									! if there was a problem, try a shorter time interval
				!write(*,*) 'Trying to continue integration from the previous time point with a shorter time interval.'
				istate = 1
				c = ci
				t0 = ti
				tt(1) = t0+(tfin-t0)/1.d2
				tt(2) = tfin
				do j = 1, 2
					call DVODE (fevalD,neqn,c,t0,tt(j),itol,rtol,atol,&	! calling the solver
						& itask,istate,iopt,work,lwork,iwork,liwork,&
						& jevalD,mf,parameters,ipar)
					if ((istate .ne. 2) .or. (minval(c) .lt. negtol)) then
						write(*,*) 'ERROR: returned istate =', istate
						write(*,*) 'Exiting the driver'
						t0 = t00
						c = c00
						j_out = 0.d0
						return									! if there is still a problem, give up
					else
						ok = .true.
					end if
				end do
			end if

			call fix_fitted_c(c,ci,ipar)			! Fix the concentrations of the fitted species

			! Check the convergence for the steady state case
			if (solve_ss) then
				if ((maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol))) .le. sstol) .and. (t0-t00 .ge. sstimetot)) then
					!write(*,*) 'converged at t=',t0
					exit		! end the simulation when a steady state is reached
				end if
				if ((t0.ge.t_max) .or. (i.eq.nt_tot)) then
                !    write(*,*) "#########################"
				!	write(*,*) "### Did not converge! ###"
                !    write(*,*) "#########################"
					!write(*,*) "monomer concentrations: ", c(n_monomers)*1.d-6
					t0 = t00
					c = c00
					j_out = 0.d0
					return
				end if
			end if

		end do

	else

		t_iter_tot = t00
        tt(1) = t_iter

		if (solve_ss) then
            t_iter_temp = 0.d0
            save_ci_temp = .true.
        end if

		do while (t_iter_tot .lt. t_max)

            ci = c						! save the previous concentrations
			ti = t_iter_tot

			if (solve_ss .and. save_ci_temp) then
				ci_temp = c
				save_ci_temp = .false.
			end if

			tt(1) = min(tt(1),t_max-t_iter_tot)
            tt(2) = tt(1)

			!write(*,*) 'Trying a step of ',tt(1),' s'
			ok = .true.
            call fevalD(neqn,ti,ci,f_out,parameters,ipar)
            c = ci + f_out*tt(1)
            t_iter_tot = t_iter_tot + tt(1)
            if (minval(c) .lt. negtol) then								! checking all concentrations > 0 - if not, print error
				!write(*,*) 'Negative concentrations at t = ',t_iter_tot
				ok = .false.
            elseif (maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol))) .gt. chmax) then	! large changes - decrease the time interval
                !write(*,*) 'Large changes in concentrations at t = ',t_iter_tot
                !nind = maxloc(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol)),1)
                !write(*,*) 'Largest rel. change at ',nind,ci(nind),c(nind)
                ok = .false.
			end if
            ! if there was a problem, try a shorter time interval
			if (.not. ok) then
				!write(*,*) 'Trying to continue integration from the previous time point with a shorter interval at t = ',t_iter_tot
                do j = 1,10
                    c = ci
                    t_iter_tot = ti
                    tt(1) = tt(2)/2.d0**real(j,kind=kind(1.d0))
                    c = ci + f_out*tt(1)
                    t_iter_tot = t_iter_tot + tt(1)
                    !write(*,*) 'Done try #',j,' to decrease the time interval'
                    ! if there is still a problem, give up after a few tries, otherwise continue
				    if ((.not. minval(c) .lt. negtol) .and.&
						&(.not. maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol)))&
						&.gt. chmax)) then
					    ok = .true.
                        exit
				    end if
                end do
                if (.not. ok) then
                    write(*,*) 'Iteration still not ok at t = ',t_iter_tot
					write(*,*) 'Exiting the driver'
					t0 = t00
					c = c00
					j_out = 0.d0
					return
                !else
                    !write(*,*) 'Succeeded to continue iteration by decreasing the interval'
                end if
			end if

			call fix_fitted_c(c,ci,ipar)			! Fix the concentrations of the fitted species

            !write(*,*) 'Finished one iteration step with time interval ',tt(1),' s'

			if (solve_ss) then
				t_iter_temp = t_iter_temp + tt(1)
				! Check the convergence
				if (t_iter_temp .ge. sstimech) then
					if ((maxval(abs((c(1:nclust)-ci_temp(1:nclust))/max(ci_temp(1:nclust),chtol))) .le. sstol) .and. (t_iter_tot-t00 .ge. sstimetot)) then
                        ! If the concentrations haven't changed during the max. value given for t_iter_temp,
                        ! consider them converged...
			            !write(*,*) 'Converged at t=',t_iter_tot
			            exit		! end the simulation when a steady state is reached
				    else
					    ! ... otherwise save the current concentrations and
					    ! check again after t_iter_temp has reached the given value
					    t_iter_temp = 0.d0
					    save_ci_temp = .true.
		            end if
			    end if
			    if (t_iter_tot .ge. t_max) then
                    write(*,*) "#########################"
                    write(*,*) "### Did not converge! ###"
                    write(*,*) "#########################"
				    !write(*,*) "monomer concentrations: ", c(n_monomers)*1.d-6
				    t0 = t00
					c = c00
					j_out = 0.d0
				    return
			    end if
			end if

			! Increase the time step if the changes appear to be small
            if (maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol))) .lt. chmin) then
                !write(*,*) 'Increasing the iteration interval at t = ',t_iter_tot
                tt(1) = min(tt(1)+10.d0,2.d0*tt(1))
            end if

		end do

		t0 = t_iter_tot
        t_iter = tt(1)

	end if

	! Formation rate
	if (solve_ss) then
		! Steady state: flux out of the system
		call formationD(neqn,c,j_tot,j_by_charge,j_by_cluster,j_all,ipar,parameters)
		j_out = (/j_tot/)
		!j_out = (/j_by_charge(1)+j_by_charge(4), j_by_charge(2), j_by_charge(3)/)
	else
		! Fixed time step: clusters grown out of the system per integrated time
		j_out = c(nout_neu)/(t0-t00)
		!j_out = c((/nout_neu, nout_neg, nout_pos/))/(t0-t00)
		! Add the initial out-of-system population back to c
		c(nout_neu) = c_out+c(nout_neu)
		!c((/nout_neu, nout_neg, nout_pos/)) = c_out+c((/nout_neu, nout_neg, nout_pos/))
	end	if

end subroutine acdc_driver

subroutine fix_fitted_c(c,ci,ipar)

	implicit none

	real(kind(1.d0)), intent(inout) :: c(neqn)			! solved concentrations after intagration -> fixed concentrations (1/m^3)
	real(kind(1.d0)), intent(in) :: ci(neqn)			! initial concentrations before integration
	integer, intent(in) :: ipar(4)						! parameters for re-calling settings after they've been changed

	logical, save :: firstcall = .true.
	real(kind(1.d0)) :: source(neqn)
	logical :: isconst(neqn)
	integer, save :: fitted(nclust,0:nclust) = 0
	integer :: k, n
	real(kind(1.d0)) :: c_fit


	if (firstcall) then
		firstcall = .false.
		call sources_and_constants(source,isconst,fitted,ipar) ! get the fitting settings (variable "fitted")
	end if

	do n = 1, fitted(1,0)	! loop over the clusters whose concentrations are fitted (e.g. sum([1AnD])=const.)
		k = fitted(n+1,0)	! this concentration is fitted
		if (ipar(1) .eq. 0) then
			c_fit = ci(k)	! if the settings have changed, the sum([1AnD])=c_fit has been fed in as [1A] in the driver call;
							! feval will change ipar(1) to 1 when the integration begins...
		else
			c_fit = sum(ci((/k, fitted(n+1,2:fitted(n+1,1)+1)/))) ! ...after which c_fit is determined from the previous concentrations
		end if
		c(k) = c_fit-sum(c(fitted(n+1,2:fitted(n+1,1)+1)))
		if (c(k).lt.0.d0) then
			c(k) = 0.d0
			c(fitted(n+1,2:fitted(n+1,1)+1)) = c(fitted(n+1,2:fitted(n+1,1)+1))*&
				&c_fit/sum(c(fitted(n+1,2:fitted(n+1,1)+1)))
		end if
	end do

end subroutine fix_fitted_c


end module driverD
