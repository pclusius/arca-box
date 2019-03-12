!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!
! Cluster 1: 1A
! Cluster 2: 2A
! Cluster 3: 1D
! Cluster 4: 1A1D
! Cluster 5: 2A1D
! Cluster 6: 3A1D
! Cluster 7: 1A2D
! Cluster 8: 2A2D
! Cluster 9: 3A2D
! Cluster 10: 4A2D
! Cluster 11: 2A3D
! Cluster 12: 3A3D
! Cluster 13: 4A3D
! Cluster 14: 3A4D
! Cluster 15: 4A4D
! 16 is for source fluxes
! 17 is for coagulation losses
! 18 is for wall losses
! 19 is for dilution losses
! 20 is for recombination of positive and negative charger ions with each other
! 21 is for collisions that lead succesfully out of the system as neutrals
! 22 is for collisions that lead out of the system, but the resulting cluster is brought back
!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!

! differential equations: f = dc/dt
subroutine feval(neqn,t,c,f,coef,ipar)
	implicit none
	integer, parameter :: nclust = 15
	integer :: neqn, ipar(4), i, j, k, n
	real(kind(1.d0)) :: t,c(neqn),f(neqn),coef(2)
	logical, save :: isconst(22) = .false.
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,22)=0.d0,coef_lin(22,22,nclust)=0.d0,source(22)=0.d0
	integer, save :: ind_quad_loss(22,0:32)=0,ind_quad_form(22,0:112)=0
	integer, save :: ind_quad_loss_extra(22,0:24)=0,ind_quad_form_extra(22,0:45)=0
	integer, save :: ind_lin_loss(22,0:14)=0,ind_lin_form(22,0:30)=0,fitted(15,0:15)=0
	real(kind(1.d0)), save :: n_quad_form_extra(22,0:15)=0.d0
	real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)


	! the parameters are read at the very beginning and eg. if source terms or collision rates have changed
	if (ipar(1) .eq. 0) then
		! after this the same values are used until some other routine tells otherwise
		ipar(1) = 1
		call initialize_parameters(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
		&	ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)
		n_quad_form_extra(:,1:15) = real(ind_quad_form_extra(:,3:45:3),kind=kind(1.d0))
	end if

	f = 0.d0
	f = f + source ! add source term
	do i=1, neqn ! loop over all the clusters + other fluxes
		if (.not. isconst(i)) then
			! first calculate coefficients for all loss terms
			do n=1, 2*ind_quad_loss(i,0)-1, 2 ! loop over all collisions removing this cluster
				f(i) = f(i)-coef_quad(i,ind_quad_loss(i,n),ind_quad_loss(i,n+1))*c(ind_quad_loss(i,n))
			end do
			do n=1, 2*ind_lin_loss(i,0)-1, 2 ! loop over all evaporations + wall losses etc. removing this cluster
				f(i) = f(i)-coef_lin(ind_lin_loss(i,n),ind_lin_loss(i,n+1),i)
			end do
			f(i) = f(i)*c(i) ! multiplying loss coefficients with current concentration
			! then add all terms that form this cluster
			do n=1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions forming this cluster
				f(i) = f(i)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&
				&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))
			end do
			do n=1, 3*ind_quad_form_extra(i,0)-2, 3 ! loop over all collisions forming this cluster as an extra product
				f(i) = f(i)+coef_quad(ind_quad_form_extra(i,n),ind_quad_form_extra(i,n+1),i)*&
				&      n_quad_form_extra(i,(n+2)/3)*c(ind_quad_form_extra(i,n))*c(ind_quad_form_extra(i,n+1))
			end do
			do n=1, 2*ind_lin_form(i,0)-1, 2 ! loop over all evaporations forming this cluster
				f(i) = f(i)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))
			end do
		end if
	end do
	do n=1, fitted(1,0) ! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)
		f(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations
		do i=1, fitted(n+1,1) ! loop over the clusters used for the fitting
			f(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i+1)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)
		end do
	end do

end subroutine feval

!-----------------------------------------------------------

! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)
! not using this, since the solver works faster using a numerical jacobian...
subroutine jevalD(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)
	implicit none
	integer :: ldim,neqn,ierr,ipar(4),i,j,k,n
	real(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn),coef(2),ml,mu

end subroutine jevalD

!-----------------------------------------------------------

subroutine formation(neqn,c,j_tot,j_by_charge,j_by_cluster,j_all,ipar,coef)
	implicit none
	integer, parameter :: nclust = 15
	integer :: neqn, i, n, ipar(4)
	real(kind(1.d0)) :: c(neqn), j_add, j_tot,j_by_charge(4),j_by_cluster(neqn),j_all(neqn,4),coef(2)
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,22)=0.d0,coef_lin(22,22,nclust)=0.d0,source(22)=0.d0
	integer, save :: charge(nclust)=0,ind_quad_loss(22,0:32)=0,ind_quad_form(22,0:112)=0
	integer, save :: ind_quad_loss_extra(22,0:24)=0,ind_quad_form_extra(22,0:45)=0
	integer, save :: ind_lin_loss(22,0:14)=0,ind_lin_form(22,0:30)=0,fitted(nclust,0:nclust)=0
	logical, save :: isconst(22)=.false.

	if (ipar(3) .eq. 0) then
		ipar(3) = 1
		call initialize_parameters(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)
		! cluster charges
		charge = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	end if
	j_tot = 0.d0			! total formation rate
	j_by_charge = 0.d0		! contributions of neutrals, negatives, positives and recombinations
	j_by_cluster = 0.d0		! contributions of different clusters
	j_all = 0.d0			! contributions of different clusters and different charging states
	do n=1, 2*ind_quad_form(21,0)-1, 2 ! loop over all collisions leading out of the system
		j_add = coef_quad(ind_quad_form(21,n),ind_quad_form(21,n+1),21)*c(ind_quad_form(21,n))*c(ind_quad_form(21,n+1))
		j_tot = j_tot + j_add
		j_by_charge(1) = j_by_charge(1) + j_add
		j_by_cluster(ind_quad_form(21,n)) = j_by_cluster(ind_quad_form(21,n)) + j_add
		j_by_cluster(ind_quad_form(21,n+1)) = j_by_cluster(ind_quad_form(21,n+1)) + j_add
		j_all(ind_quad_form(21,n),1) = j_all(ind_quad_form(21,n),1) + j_add
		j_all(ind_quad_form(21,n+1),1) = j_all(ind_quad_form(21,n+1),1) + j_add
	end do

end subroutine formation

!-----------------------------------------------------------

subroutine initialize_parameters(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)

	use monomer_settings, only : sources_and_constants

	implicit none
	integer, parameter :: neqn=22, nclust=15
	logical :: isconst(22)
	real(kind(1.d0)) :: coef_quad(nclust,nclust,22),coef_lin(22,22,nclust)
	real(kind(1.d0)) :: source(22)
	integer :: ind_quad_loss(22,0:32),ind_quad_form(22,0:112),ind_quad_loss_extra(22,0:24),ind_quad_form_extra(22,0:45)
	integer :: ind_lin_loss(22,0:14),ind_lin_form(22,0:30)
	integer :: fitted(nclust,0:nclust)
	real(kind(1.d0)) :: coef(2)
	integer :: ipar(4)

	call sources_and_constants(source,isconst,fitted,ipar)

	ind_quad_loss = 0
	ind_quad_form = 0
	ind_lin_loss = 0
	ind_lin_form = 0
	ind_quad_loss_extra = 0
	ind_quad_form_extra = 0

	! ind_quad_loss(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the collisions
	!  i + j1 -> k1 etc. through which cluster i is lost

	! ind_quad_form(k,0:2n) = (/n,i1,j1,i2,j2,...,in,jn/) gives the cluster numbers for all the collisions
	!  i1 + j1 -> k etc. through which cluster k is formed

	! ind_lin_loss(k,0:2n) = (/n,i1,j1,i2,j2,...,lossn,lossn/) gives the cluster numbers for all the evaporations
	!  k -> i1 + j1 etc. and other losses k -> wall etc. through which cluster k is lost

	! ind_lin_form(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the evaporations
	!  k1 -> i + j1 etc. through which cluster i is formed

	! ind_quad_loss_extra(i,0:3n) = (/n,j1,k1,c1,...,jn,kn,cn/) gives the cluster numbers and coefficients
	!  i + j1 -> c1*k1 etc. for additional collision products k (e.g. monomers from the boundary)
	!  in collisions where cluster i is lost

	! ind_quad_form_extra(k,0:2n) = (/n,i1,j1,c1,...,in,jn,cn/) gives the cluster numbers and coefficients
	!  i1 + j1 -> c1*k etc. for additional ways of forming k (e.g. as a monomer from the boundary)


	! Cluster 1: 1A
	! Derivative not used and concentration set to be constant
	isconst(1) = .true.

	! Cluster 2: 2A
	ind_quad_loss(2,0:30) = (/ 15, 2,2, 2,2, 3,5, 4,6, 5,6, 6,6, 7,9, 8,10, 9,10&
			&, 10,10, 11,13, 12,13, 13,13, 14,21, 15,21 /)
	ind_quad_loss_extra(2,0:24) = (/ 8, 2,1,2, 2,1,2, 5,1,1, 6,1,2, 9,1,1, 10,1,2, 12,1,1, 13,1,2 /)
	ind_quad_form(2,0:4) = (/ 2, 1,1, 2,2 /)
	ind_lin_loss(2,0:4) = (/ 2, 1,1, 17,17 /)
	ind_lin_form(2,0:10) = (/ 5, 3,5, 4,6, 7,9, 8,10, 11,13 /)

	! Cluster 3: 1D
	! Derivative not used and concentration set to be constant
	isconst(3) = .true.

	! Cluster 4: 1A1D
	ind_quad_loss(4,0:32) = (/ 16, 1,5, 2,6, 3,7, 4,8, 4,8, 5,9, 6,10, 7,11, 8,12&
			&, 9,13, 10,13, 11,14, 12,15, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(4,0:3) = (/ 1, 10,1,1 /)
	ind_quad_form(4,0:2) = (/ 1, 1,3 /)
	ind_lin_loss(4,0:4) = (/ 2, 1,3, 17,17 /)
	ind_lin_form(4,0:24) = (/ 12, 1,5, 2,6, 3,7, 4,8, 4,8, 5,9, 6,10, 7,11, 8,12&
			&, 9,13, 11,14, 12,15 /)

	! Cluster 5: 2A1D
	ind_quad_loss(5,0:32) = (/ 16, 1,6, 2,6, 3,8, 4,9, 5,10, 5,10, 6,10, 7,12, 8,13&
			&, 9,13, 10,13, 11,15, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(5,0:12) = (/ 4, 2,1,1, 6,1,1, 9,1,1, 10,1,2 /)
	ind_quad_form(5,0:4) = (/ 2, 1,4, 2,3 /)
	ind_lin_loss(5,0:6) = (/ 3, 1,4, 2,3, 17,17 /)
	ind_lin_form(5,0:16) = (/ 8, 1,6, 3,8, 4,9, 5,10, 5,10, 7,12, 8,13, 11,15 /)

	! Cluster 6: 3A1D
	ind_quad_loss(6,0:30) = (/ 15, 2,6, 3,9, 4,10, 5,10, 6,10, 6,10, 7,13, 8,13, 9,13&
			&, 10,13, 11,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(6,0:21) = (/ 7, 2,1,2, 5,1,1, 6,1,2, 6,1,2, 8,1,1, 9,1,2, 10,1,3 /)
	ind_quad_form(6,0:8) = (/ 4, 1,5, 2,4, 2,5, 2,6 /)
	ind_lin_loss(6,0:6) = (/ 3, 1,5, 2,4, 17,17 /)
	ind_lin_form(6,0:6) = (/ 3, 3,9, 4,10, 7,13 /)

	! Cluster 7: 1A2D
	ind_quad_loss(7,0:30) = (/ 15, 1,8, 2,9, 4,11, 5,12, 6,13, 7,11, 7,11, 8,14, 9,15&
			&, 10,21, 11,14, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(7,0:9) = (/ 3, 7,3,1, 7,3,1, 11,3,1 /)
	ind_quad_form(7,0:2) = (/ 1, 3,4 /)
	ind_lin_loss(7,0:4) = (/ 2, 3,4, 17,17 /)
	ind_lin_form(7,0:14) = (/ 7, 1,8, 2,9, 4,11, 5,12, 6,13, 8,14, 9,15 /)

	! Cluster 8: 2A2D
	ind_quad_loss(8,0:32) = (/ 16, 1,9, 2,10, 3,11, 4,12, 5,13, 6,13, 7,14, 8,15, 8,15&
			&, 9,21, 10,21, 11,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(8,0:3) = (/ 1, 6,1,1 /)
	ind_quad_form(8,0:6) = (/ 3, 1,7, 3,5, 4,4 /)
	ind_lin_loss(8,0:8) = (/ 4, 1,7, 3,5, 4,4, 17,17 /)
	ind_lin_form(8,0:16) = (/ 8, 1,9, 2,10, 3,11, 4,12, 5,13, 7,14, 8,15, 8,15 /)

	! Cluster 9: 3A2D
	ind_quad_loss(9,0:32) = (/ 16, 1,10, 2,10, 3,12, 4,13, 5,13, 6,13, 7,15, 8,21, 9,21&
			&, 9,21, 10,21, 11,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(9,0:9) = (/ 3, 2,1,1, 5,1,1, 6,1,2 /)
	ind_quad_form(9,0:8) = (/ 4, 1,8, 2,7, 3,6, 4,5 /)
	ind_lin_loss(9,0:10) = (/ 5, 1,8, 2,7, 3,6, 4,5, 17,17 /)
	ind_lin_form(9,0:8) = (/ 4, 1,10, 3,12, 4,13, 7,15 /)

	! Cluster 10: 4A2D
	ind_quad_loss(10,0:30) = (/ 15, 2,10, 3,13, 4,13, 5,13, 6,13, 7,21, 8,21, 9,21, 10,21&
			&, 10,21, 11,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(10,0:12) = (/ 4, 2,1,2, 4,1,1, 5,1,2, 6,1,3 /)
	ind_quad_form(10,0:16) = (/ 8, 1,9, 2,8, 2,9, 2,10, 4,6, 5,5, 5,6, 6,6 /)
	ind_lin_loss(10,0:10) = (/ 5, 1,9, 2,8, 4,6, 5,5, 17,17 /)
	ind_lin_form(10,0:2) = (/ 1, 3,13 /)

	! Cluster 11: 2A3D
	ind_quad_loss(11,0:30) = (/ 15, 1,12, 2,13, 4,14, 5,15, 6,21, 7,14, 8,21, 9,21, 10,21&
			&, 11,21, 11,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(11,0:3) = (/ 1, 7,3,1 /)
	ind_quad_form(11,0:6) = (/ 3, 3,8, 4,7, 7,7 /)
	ind_lin_loss(11,0:6) = (/ 3, 3,8, 4,7, 17,17 /)
	ind_lin_form(11,0:8) = (/ 4, 1,12, 2,13, 4,14, 5,15 /)

	! Cluster 12: 3A3D
	ind_quad_loss(12,0:32) = (/ 16, 1,13, 2,13, 3,14, 4,15, 5,21, 6,21, 7,21, 8,21, 9,21&
			&, 10,21, 11,21, 12,21, 12,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(12,0:3) = (/ 1, 2,1,1 /)
	ind_quad_form(12,0:8) = (/ 4, 1,11, 3,9, 4,8, 5,7 /)
	ind_lin_loss(12,0:10) = (/ 5, 1,11, 3,9, 4,8, 5,7, 17,17 /)
	ind_lin_form(12,0:6) = (/ 3, 1,13, 3,14, 4,15 /)

	! Cluster 13: 4A3D
	ind_quad_loss(13,0:30) = (/ 15, 2,13, 3,15, 4,21, 5,21, 6,21, 7,21, 8,21, 9,21, 10,21&
			&, 11,21, 12,21, 13,21, 13,21, 14,21, 15,21 /)
	ind_quad_loss_extra(13,0:3) = (/ 1, 2,1,2 /)
	ind_quad_form(13,0:28) = (/ 14, 1,12, 2,11, 2,12, 2,13, 3,10, 4,9, 4,10, 5,8, 5,9&
			&, 5,10, 6,7, 6,8, 6,9, 6,10 /)
	ind_lin_loss(13,0:14) = (/ 7, 1,12, 2,11, 3,10, 4,9, 5,8, 6,7, 17,17 /)
	ind_lin_form(13,0:2) = (/ 1, 3,15 /)

	! Cluster 14: 3A4D
	ind_quad_loss(14,0:30) = (/ 15, 1,15, 2,21, 4,21, 5,21, 6,21, 7,21, 8,21, 9,21, 10,21&
			&, 11,21, 12,21, 13,21, 14,21, 14,21, 15,21 /)
	ind_quad_form(14,0:8) = (/ 4, 3,12, 4,11, 7,8, 7,11 /)
	ind_lin_loss(14,0:8) = (/ 4, 3,12, 4,11, 7,8, 17,17 /)
	ind_lin_form(14,0:2) = (/ 1, 1,15 /)

	! Cluster 15: 4A4D
	ind_quad_loss(15,0:32) = (/ 16, 1,21, 2,21, 3,21, 4,21, 5,21, 6,21, 7,21, 8,21, 9,21&
			&, 10,21, 11,21, 12,21, 13,21, 14,21, 15,21, 15,21 /)
	ind_quad_form(15,0:12) = (/ 6, 1,14, 3,13, 4,12, 5,11, 7,9, 8,8 /)
	ind_lin_loss(15,0:14) = (/ 7, 1,14, 3,13, 4,12, 5,11, 7,9, 8,8, 17,17 /)

	! Cluster 17: coag
	ind_lin_form(17,0:30) = (/ 15, 17,1, 17,2, 17,3, 17,4, 17,5, 17,6, 17,7, 17,8, 17,9&
			&, 17,10, 17,11, 17,12, 17,13, 17,14, 17,15 /)

	! Cluster 21: out_neu
	ind_quad_form(21,0:112) = (/ 56, 1,15, 2,14, 2,15, 3,15, 4,13, 4,14, 4,15, 5,12, 5,13&
			&, 5,14, 5,15, 6,11, 6,12, 6,13, 6,14, 6,15, 7,10, 7,12, 7,13&
			&, 7,14, 7,15, 8,9, 8,10, 8,11, 8,12, 8,13, 8,14, 8,15, 9,9&
			&, 9,10, 9,11, 9,12, 9,13, 9,14, 9,15, 10,10, 10,11, 10,12, 10,13&
			&, 10,14, 10,15, 11,11, 11,12, 11,13, 11,14, 11,15, 12,12, 12,13, 12,14&
			&, 12,15, 13,13, 13,14, 13,15, 14,14, 14,15, 15,15 /)

	call get_rate_coefs(coef_quad,coef_lin,coef)

end subroutine initialize_parameters

!-----------------------------------------------------------

subroutine get_rate_coefs(coef_quad,coef_lin,coef)
	implicit none
	real(kind(1.d0)) :: coef_quad(15,15,22),coef_lin(22,22,15)
	real(kind(1.d0)) :: K(15,15),E(15,15),cs(15)
	real(kind(1.d0)) :: coef(2)

	coef_quad = 0.d0
	coef_lin = 0.d0
	call get_coll(K,coef(1))
	call get_evap(E,K,coef(1))
	call get_losses(cs)
	cs = coef(2)*cs


	coef_quad(1,1,2) = 0.5d0*K(1,1)	! 1A + 1A -> 2A
	coef_lin(1,1,2) = E(1,1)	! 2A -> 1A + 1A

	coef_quad(2,2,2) = 0.5d0*K(2,2)	! 2A + 2A -> boundary -> 2A
	coef_quad(2,2,1) = 0.5d0*K(2,2)	! 2A + 2A -> boundary -> 1A

	coef_quad(3,1,4) = K(1,3)	! 1D + 1A -> 1A1D
	coef_quad(1,3,4) = K(1,3)
	coef_lin(3,1,4) = E(1,3)	! 1A1D -> 1D + 1A
	coef_lin(1,3,4) = E(1,3)
	coef_quad(3,2,5) = K(2,3)	! 1D + 2A -> 2A1D
	coef_quad(2,3,5) = K(2,3)
	coef_lin(3,2,5) = E(2,3)	! 2A1D -> 1D + 2A
	coef_lin(2,3,5) = E(2,3)

	coef_quad(4,1,5) = K(1,4)	! 1A1D + 1A -> 2A1D
	coef_quad(1,4,5) = K(1,4)
	coef_lin(4,1,5) = E(1,4)	! 2A1D -> 1A1D + 1A
	coef_lin(1,4,5) = E(1,4)
	coef_quad(4,2,6) = K(2,4)	! 1A1D + 2A -> 3A1D
	coef_quad(2,4,6) = K(2,4)
	coef_lin(4,2,6) = E(2,4)	! 3A1D -> 1A1D + 2A
	coef_lin(2,4,6) = E(2,4)
	coef_quad(4,3,7) = K(3,4)	! 1A1D + 1D -> 1A2D
	coef_quad(3,4,7) = K(3,4)
	coef_lin(4,3,7) = E(3,4)	! 1A2D -> 1A1D + 1D
	coef_lin(3,4,7) = E(3,4)
	coef_quad(4,4,8) = 0.5d0*K(4,4)	! 1A1D + 1A1D -> 2A2D
	coef_lin(4,4,8) = E(4,4)	! 2A2D -> 1A1D + 1A1D

	coef_quad(5,1,6) = K(1,5)	! 2A1D + 1A -> 3A1D
	coef_quad(1,5,6) = K(1,5)
	coef_lin(5,1,6) = E(1,5)	! 3A1D -> 2A1D + 1A
	coef_lin(1,5,6) = E(1,5)
	coef_quad(5,2,6) = K(2,5)	! 2A1D + 2A -> boundary -> 3A1D
	coef_quad(2,5,6) = K(2,5)
	coef_quad(5,2,1) = K(2,5)	! 2A1D + 2A -> boundary -> 1A
	coef_quad(2,5,1) = K(2,5)
	coef_quad(5,3,8) = K(3,5)	! 2A1D + 1D -> 2A2D
	coef_quad(3,5,8) = K(3,5)
	coef_lin(5,3,8) = E(3,5)	! 2A2D -> 2A1D + 1D
	coef_lin(3,5,8) = E(3,5)
	coef_quad(5,4,9) = K(4,5)	! 2A1D + 1A1D -> 3A2D
	coef_quad(4,5,9) = K(4,5)
	coef_lin(5,4,9) = E(4,5)	! 3A2D -> 2A1D + 1A1D
	coef_lin(4,5,9) = E(4,5)
	coef_quad(5,5,10) = 0.5d0*K(5,5)	! 2A1D + 2A1D -> 4A2D
	coef_lin(5,5,10) = E(5,5)	! 4A2D -> 2A1D + 2A1D

	coef_quad(6,2,6) = K(2,6)	! 3A1D + 2A -> boundary -> 3A1D
	coef_quad(2,6,6) = K(2,6)
	coef_quad(6,2,1) = K(2,6)	! 3A1D + 2A -> boundary -> 1A
	coef_quad(2,6,1) = K(2,6)
	coef_quad(6,3,9) = K(3,6)	! 3A1D + 1D -> 3A2D
	coef_quad(3,6,9) = K(3,6)
	coef_lin(6,3,9) = E(3,6)	! 3A2D -> 3A1D + 1D
	coef_lin(3,6,9) = E(3,6)
	coef_quad(6,4,10) = K(4,6)	! 3A1D + 1A1D -> 4A2D
	coef_quad(4,6,10) = K(4,6)
	coef_lin(6,4,10) = E(4,6)	! 4A2D -> 3A1D + 1A1D
	coef_lin(4,6,10) = E(4,6)
	coef_quad(6,5,10) = K(5,6)	! 3A1D + 2A1D -> boundary -> 4A2D
	coef_quad(5,6,10) = K(5,6)
	coef_quad(6,5,1) = K(5,6)	! 3A1D + 2A1D -> boundary -> 1A
	coef_quad(5,6,1) = K(5,6)
	coef_quad(6,6,10) = 0.5d0*K(6,6)	! 3A1D + 3A1D -> boundary -> 4A2D
	coef_quad(6,6,1) = 0.5d0*K(6,6)	! 3A1D + 3A1D -> boundary -> 1A

	coef_quad(7,1,8) = K(1,7)	! 1A2D + 1A -> 2A2D
	coef_quad(1,7,8) = K(1,7)
	coef_lin(7,1,8) = E(1,7)	! 2A2D -> 1A2D + 1A
	coef_lin(1,7,8) = E(1,7)
	coef_quad(7,2,9) = K(2,7)	! 1A2D + 2A -> 3A2D
	coef_quad(2,7,9) = K(2,7)
	coef_lin(7,2,9) = E(2,7)	! 3A2D -> 1A2D + 2A
	coef_lin(2,7,9) = E(2,7)
	coef_quad(7,4,11) = K(4,7)	! 1A2D + 1A1D -> 2A3D
	coef_quad(4,7,11) = K(4,7)
	coef_lin(7,4,11) = E(4,7)	! 2A3D -> 1A2D + 1A1D
	coef_lin(4,7,11) = E(4,7)
	coef_quad(7,5,12) = K(5,7)	! 1A2D + 2A1D -> 3A3D
	coef_quad(5,7,12) = K(5,7)
	coef_lin(7,5,12) = E(5,7)	! 3A3D -> 1A2D + 2A1D
	coef_lin(5,7,12) = E(5,7)
	coef_quad(7,6,13) = K(6,7)	! 1A2D + 3A1D -> 4A3D
	coef_quad(6,7,13) = K(6,7)
	coef_lin(7,6,13) = E(6,7)	! 4A3D -> 1A2D + 3A1D
	coef_lin(6,7,13) = E(6,7)
	coef_quad(7,7,11) = 0.5d0*K(7,7)	! 1A2D + 1A2D -> boundary -> 2A3D
	coef_quad(7,7,3) = 0.5d0*K(7,7)	! 1A2D + 1A2D -> boundary -> 1D

	coef_quad(8,1,9) = K(1,8)	! 2A2D + 1A -> 3A2D
	coef_quad(1,8,9) = K(1,8)
	coef_lin(8,1,9) = E(1,8)	! 3A2D -> 2A2D + 1A
	coef_lin(1,8,9) = E(1,8)
	coef_quad(8,2,10) = K(2,8)	! 2A2D + 2A -> 4A2D
	coef_quad(2,8,10) = K(2,8)
	coef_lin(8,2,10) = E(2,8)	! 4A2D -> 2A2D + 2A
	coef_lin(2,8,10) = E(2,8)
	coef_quad(8,3,11) = K(3,8)	! 2A2D + 1D -> 2A3D
	coef_quad(3,8,11) = K(3,8)
	coef_lin(8,3,11) = E(3,8)	! 2A3D -> 2A2D + 1D
	coef_lin(3,8,11) = E(3,8)
	coef_quad(8,4,12) = K(4,8)	! 2A2D + 1A1D -> 3A3D
	coef_quad(4,8,12) = K(4,8)
	coef_lin(8,4,12) = E(4,8)	! 3A3D -> 2A2D + 1A1D
	coef_lin(4,8,12) = E(4,8)
	coef_quad(8,5,13) = K(5,8)	! 2A2D + 2A1D -> 4A3D
	coef_quad(5,8,13) = K(5,8)
	coef_lin(8,5,13) = E(5,8)	! 4A3D -> 2A2D + 2A1D
	coef_lin(5,8,13) = E(5,8)
	coef_quad(8,6,13) = K(6,8)	! 2A2D + 3A1D -> boundary -> 4A3D
	coef_quad(6,8,13) = K(6,8)
	coef_quad(8,6,1) = K(6,8)	! 2A2D + 3A1D -> boundary -> 1A
	coef_quad(6,8,1) = K(6,8)
	coef_quad(8,7,14) = K(7,8)	! 2A2D + 1A2D -> 3A4D
	coef_quad(7,8,14) = K(7,8)
	coef_lin(8,7,14) = E(7,8)	! 3A4D -> 2A2D + 1A2D
	coef_lin(7,8,14) = E(7,8)
	coef_quad(8,8,15) = 0.5d0*K(8,8)	! 2A2D + 2A2D -> 4A4D
	coef_lin(8,8,15) = E(8,8)	! 4A4D -> 2A2D + 2A2D

	coef_quad(9,1,10) = K(1,9)	! 3A2D + 1A -> 4A2D
	coef_quad(1,9,10) = K(1,9)
	coef_lin(9,1,10) = E(1,9)	! 4A2D -> 3A2D + 1A
	coef_lin(1,9,10) = E(1,9)
	coef_quad(9,2,10) = K(2,9)	! 3A2D + 2A -> boundary -> 4A2D
	coef_quad(2,9,10) = K(2,9)
	coef_quad(9,2,1) = K(2,9)	! 3A2D + 2A -> boundary -> 1A
	coef_quad(2,9,1) = K(2,9)
	coef_quad(9,3,12) = K(3,9)	! 3A2D + 1D -> 3A3D
	coef_quad(3,9,12) = K(3,9)
	coef_lin(9,3,12) = E(3,9)	! 3A3D -> 3A2D + 1D
	coef_lin(3,9,12) = E(3,9)
	coef_quad(9,4,13) = K(4,9)	! 3A2D + 1A1D -> 4A3D
	coef_quad(4,9,13) = K(4,9)
	coef_lin(9,4,13) = E(4,9)	! 4A3D -> 3A2D + 1A1D
	coef_lin(4,9,13) = E(4,9)
	coef_quad(9,5,13) = K(5,9)	! 3A2D + 2A1D -> boundary -> 4A3D
	coef_quad(5,9,13) = K(5,9)
	coef_quad(9,5,1) = K(5,9)	! 3A2D + 2A1D -> boundary -> 1A
	coef_quad(5,9,1) = K(5,9)
	coef_quad(9,6,13) = K(6,9)	! 3A2D + 3A1D -> boundary -> 4A3D
	coef_quad(6,9,13) = K(6,9)
	coef_quad(9,6,1) = K(6,9)	! 3A2D + 3A1D -> boundary -> 1A
	coef_quad(6,9,1) = K(6,9)
	coef_quad(9,7,15) = K(7,9)	! 3A2D + 1A2D -> 4A4D
	coef_quad(7,9,15) = K(7,9)
	coef_lin(9,7,15) = E(7,9)	! 4A4D -> 3A2D + 1A2D
	coef_lin(7,9,15) = E(7,9)
	coef_quad(9,8,21) = K(8,9)	! 3A2D + 2A2D -> out_neu
	coef_quad(8,9,21) = K(8,9)
	coef_quad(9,9,21) = 0.5d0*K(9,9)	! 3A2D + 3A2D -> out_neu

	coef_quad(10,2,10) = K(2,10)	! 4A2D + 2A -> boundary -> 4A2D
	coef_quad(2,10,10) = K(2,10)
	coef_quad(10,2,1) = K(2,10)	! 4A2D + 2A -> boundary -> 1A
	coef_quad(2,10,1) = K(2,10)
	coef_quad(10,3,13) = K(3,10)	! 4A2D + 1D -> 4A3D
	coef_quad(3,10,13) = K(3,10)
	coef_lin(10,3,13) = E(3,10)	! 4A3D -> 4A2D + 1D
	coef_lin(3,10,13) = E(3,10)
	coef_quad(10,4,13) = K(4,10)	! 4A2D + 1A1D -> boundary -> 4A3D
	coef_quad(4,10,13) = K(4,10)
	coef_quad(10,4,1) = K(4,10)	! 4A2D + 1A1D -> boundary -> 1A
	coef_quad(4,10,1) = K(4,10)
	coef_quad(10,5,13) = K(5,10)	! 4A2D + 2A1D -> boundary -> 4A3D
	coef_quad(5,10,13) = K(5,10)
	coef_quad(10,5,1) = K(5,10)	! 4A2D + 2A1D -> boundary -> 1A
	coef_quad(5,10,1) = K(5,10)
	coef_quad(10,6,13) = K(6,10)	! 4A2D + 3A1D -> boundary -> 4A3D
	coef_quad(6,10,13) = K(6,10)
	coef_quad(10,6,1) = K(6,10)	! 4A2D + 3A1D -> boundary -> 1A
	coef_quad(6,10,1) = K(6,10)
	coef_quad(10,7,21) = K(7,10)	! 4A2D + 1A2D -> out_neu
	coef_quad(7,10,21) = K(7,10)
	coef_quad(10,8,21) = K(8,10)	! 4A2D + 2A2D -> out_neu
	coef_quad(8,10,21) = K(8,10)
	coef_quad(10,9,21) = K(9,10)	! 4A2D + 3A2D -> out_neu
	coef_quad(9,10,21) = K(9,10)
	coef_quad(10,10,21) = 0.5d0*K(10,10)	! 4A2D + 4A2D -> out_neu

	coef_quad(11,1,12) = K(1,11)	! 2A3D + 1A -> 3A3D
	coef_quad(1,11,12) = K(1,11)
	coef_lin(11,1,12) = E(1,11)	! 3A3D -> 2A3D + 1A
	coef_lin(1,11,12) = E(1,11)
	coef_quad(11,2,13) = K(2,11)	! 2A3D + 2A -> 4A3D
	coef_quad(2,11,13) = K(2,11)
	coef_lin(11,2,13) = E(2,11)	! 4A3D -> 2A3D + 2A
	coef_lin(2,11,13) = E(2,11)
	coef_quad(11,4,14) = K(4,11)	! 2A3D + 1A1D -> 3A4D
	coef_quad(4,11,14) = K(4,11)
	coef_lin(11,4,14) = E(4,11)	! 3A4D -> 2A3D + 1A1D
	coef_lin(4,11,14) = E(4,11)
	coef_quad(11,5,15) = K(5,11)	! 2A3D + 2A1D -> 4A4D
	coef_quad(5,11,15) = K(5,11)
	coef_lin(11,5,15) = E(5,11)	! 4A4D -> 2A3D + 2A1D
	coef_lin(5,11,15) = E(5,11)
	coef_quad(11,6,21) = K(6,11)	! 2A3D + 3A1D -> out_neu
	coef_quad(6,11,21) = K(6,11)
	coef_quad(11,7,14) = K(7,11)	! 2A3D + 1A2D -> boundary -> 3A4D
	coef_quad(7,11,14) = K(7,11)
	coef_quad(11,7,3) = K(7,11)	! 2A3D + 1A2D -> boundary -> 1D
	coef_quad(7,11,3) = K(7,11)
	coef_quad(11,8,21) = K(8,11)	! 2A3D + 2A2D -> out_neu
	coef_quad(8,11,21) = K(8,11)
	coef_quad(11,9,21) = K(9,11)	! 2A3D + 3A2D -> out_neu
	coef_quad(9,11,21) = K(9,11)
	coef_quad(11,10,21) = K(10,11)	! 2A3D + 4A2D -> out_neu
	coef_quad(10,11,21) = K(10,11)
	coef_quad(11,11,21) = 0.5d0*K(11,11)	! 2A3D + 2A3D -> out_neu

	coef_quad(12,1,13) = K(1,12)	! 3A3D + 1A -> 4A3D
	coef_quad(1,12,13) = K(1,12)
	coef_lin(12,1,13) = E(1,12)	! 4A3D -> 3A3D + 1A
	coef_lin(1,12,13) = E(1,12)
	coef_quad(12,2,13) = K(2,12)	! 3A3D + 2A -> boundary -> 4A3D
	coef_quad(2,12,13) = K(2,12)
	coef_quad(12,2,1) = K(2,12)	! 3A3D + 2A -> boundary -> 1A
	coef_quad(2,12,1) = K(2,12)
	coef_quad(12,3,14) = K(3,12)	! 3A3D + 1D -> 3A4D
	coef_quad(3,12,14) = K(3,12)
	coef_lin(12,3,14) = E(3,12)	! 3A4D -> 3A3D + 1D
	coef_lin(3,12,14) = E(3,12)
	coef_quad(12,4,15) = K(4,12)	! 3A3D + 1A1D -> 4A4D
	coef_quad(4,12,15) = K(4,12)
	coef_lin(12,4,15) = E(4,12)	! 4A4D -> 3A3D + 1A1D
	coef_lin(4,12,15) = E(4,12)
	coef_quad(12,5,21) = K(5,12)	! 3A3D + 2A1D -> out_neu
	coef_quad(5,12,21) = K(5,12)
	coef_quad(12,6,21) = K(6,12)	! 3A3D + 3A1D -> out_neu
	coef_quad(6,12,21) = K(6,12)
	coef_quad(12,7,21) = K(7,12)	! 3A3D + 1A2D -> out_neu
	coef_quad(7,12,21) = K(7,12)
	coef_quad(12,8,21) = K(8,12)	! 3A3D + 2A2D -> out_neu
	coef_quad(8,12,21) = K(8,12)
	coef_quad(12,9,21) = K(9,12)	! 3A3D + 3A2D -> out_neu
	coef_quad(9,12,21) = K(9,12)
	coef_quad(12,10,21) = K(10,12)	! 3A3D + 4A2D -> out_neu
	coef_quad(10,12,21) = K(10,12)
	coef_quad(12,11,21) = K(11,12)	! 3A3D + 2A3D -> out_neu
	coef_quad(11,12,21) = K(11,12)
	coef_quad(12,12,21) = 0.5d0*K(12,12)	! 3A3D + 3A3D -> out_neu

	coef_quad(13,2,13) = K(2,13)	! 4A3D + 2A -> boundary -> 4A3D
	coef_quad(2,13,13) = K(2,13)
	coef_quad(13,2,1) = K(2,13)	! 4A3D + 2A -> boundary -> 1A
	coef_quad(2,13,1) = K(2,13)
	coef_quad(13,3,15) = K(3,13)	! 4A3D + 1D -> 4A4D
	coef_quad(3,13,15) = K(3,13)
	coef_lin(13,3,15) = E(3,13)	! 4A4D -> 4A3D + 1D
	coef_lin(3,13,15) = E(3,13)
	coef_quad(13,4,21) = K(4,13)	! 4A3D + 1A1D -> out_neu
	coef_quad(4,13,21) = K(4,13)
	coef_quad(13,5,21) = K(5,13)	! 4A3D + 2A1D -> out_neu
	coef_quad(5,13,21) = K(5,13)
	coef_quad(13,6,21) = K(6,13)	! 4A3D + 3A1D -> out_neu
	coef_quad(6,13,21) = K(6,13)
	coef_quad(13,7,21) = K(7,13)	! 4A3D + 1A2D -> out_neu
	coef_quad(7,13,21) = K(7,13)
	coef_quad(13,8,21) = K(8,13)	! 4A3D + 2A2D -> out_neu
	coef_quad(8,13,21) = K(8,13)
	coef_quad(13,9,21) = K(9,13)	! 4A3D + 3A2D -> out_neu
	coef_quad(9,13,21) = K(9,13)
	coef_quad(13,10,21) = K(10,13)	! 4A3D + 4A2D -> out_neu
	coef_quad(10,13,21) = K(10,13)
	coef_quad(13,11,21) = K(11,13)	! 4A3D + 2A3D -> out_neu
	coef_quad(11,13,21) = K(11,13)
	coef_quad(13,12,21) = K(12,13)	! 4A3D + 3A3D -> out_neu
	coef_quad(12,13,21) = K(12,13)
	coef_quad(13,13,21) = 0.5d0*K(13,13)	! 4A3D + 4A3D -> out_neu

	coef_quad(14,1,15) = K(1,14)	! 3A4D + 1A -> 4A4D
	coef_quad(1,14,15) = K(1,14)
	coef_lin(14,1,15) = E(1,14)	! 4A4D -> 3A4D + 1A
	coef_lin(1,14,15) = E(1,14)
	coef_quad(14,2,21) = K(2,14)	! 3A4D + 2A -> out_neu
	coef_quad(2,14,21) = K(2,14)
	coef_quad(14,4,21) = K(4,14)	! 3A4D + 1A1D -> out_neu
	coef_quad(4,14,21) = K(4,14)
	coef_quad(14,5,21) = K(5,14)	! 3A4D + 2A1D -> out_neu
	coef_quad(5,14,21) = K(5,14)
	coef_quad(14,6,21) = K(6,14)	! 3A4D + 3A1D -> out_neu
	coef_quad(6,14,21) = K(6,14)
	coef_quad(14,7,21) = K(7,14)	! 3A4D + 1A2D -> out_neu
	coef_quad(7,14,21) = K(7,14)
	coef_quad(14,8,21) = K(8,14)	! 3A4D + 2A2D -> out_neu
	coef_quad(8,14,21) = K(8,14)
	coef_quad(14,9,21) = K(9,14)	! 3A4D + 3A2D -> out_neu
	coef_quad(9,14,21) = K(9,14)
	coef_quad(14,10,21) = K(10,14)	! 3A4D + 4A2D -> out_neu
	coef_quad(10,14,21) = K(10,14)
	coef_quad(14,11,21) = K(11,14)	! 3A4D + 2A3D -> out_neu
	coef_quad(11,14,21) = K(11,14)
	coef_quad(14,12,21) = K(12,14)	! 3A4D + 3A3D -> out_neu
	coef_quad(12,14,21) = K(12,14)
	coef_quad(14,13,21) = K(13,14)	! 3A4D + 4A3D -> out_neu
	coef_quad(13,14,21) = K(13,14)
	coef_quad(14,14,21) = 0.5d0*K(14,14)	! 3A4D + 3A4D -> out_neu

	coef_quad(15,1,21) = K(1,15)	! 4A4D + 1A -> out_neu
	coef_quad(1,15,21) = K(1,15)
	coef_quad(15,2,21) = K(2,15)	! 4A4D + 2A -> out_neu
	coef_quad(2,15,21) = K(2,15)
	coef_quad(15,3,21) = K(3,15)	! 4A4D + 1D -> out_neu
	coef_quad(3,15,21) = K(3,15)
	coef_quad(15,4,21) = K(4,15)	! 4A4D + 1A1D -> out_neu
	coef_quad(4,15,21) = K(4,15)
	coef_quad(15,5,21) = K(5,15)	! 4A4D + 2A1D -> out_neu
	coef_quad(5,15,21) = K(5,15)
	coef_quad(15,6,21) = K(6,15)	! 4A4D + 3A1D -> out_neu
	coef_quad(6,15,21) = K(6,15)
	coef_quad(15,7,21) = K(7,15)	! 4A4D + 1A2D -> out_neu
	coef_quad(7,15,21) = K(7,15)
	coef_quad(15,8,21) = K(8,15)	! 4A4D + 2A2D -> out_neu
	coef_quad(8,15,21) = K(8,15)
	coef_quad(15,9,21) = K(9,15)	! 4A4D + 3A2D -> out_neu
	coef_quad(9,15,21) = K(9,15)
	coef_quad(15,10,21) = K(10,15)	! 4A4D + 4A2D -> out_neu
	coef_quad(10,15,21) = K(10,15)
	coef_quad(15,11,21) = K(11,15)	! 4A4D + 2A3D -> out_neu
	coef_quad(11,15,21) = K(11,15)
	coef_quad(15,12,21) = K(12,15)	! 4A4D + 3A3D -> out_neu
	coef_quad(12,15,21) = K(12,15)
	coef_quad(15,13,21) = K(13,15)	! 4A4D + 4A3D -> out_neu
	coef_quad(13,15,21) = K(13,15)
	coef_quad(15,14,21) = K(14,15)	! 4A4D + 3A4D -> out_neu
	coef_quad(14,15,21) = K(14,15)
	coef_quad(15,15,21) = 0.5d0*K(15,15)	! 4A4D + 4A4D -> out_neu

	coef_lin(17,17,15) = cs(15)	! 4A4D -> coag
	coef_lin(17,17,14) = cs(14)	! 3A4D -> coag
	coef_lin(17,17,13) = cs(13)	! 4A3D -> coag
	coef_lin(17,17,12) = cs(12)	! 3A3D -> coag
	coef_lin(17,17,11) = cs(11)	! 2A3D -> coag
	coef_lin(17,17,10) = cs(10)	! 4A2D -> coag
	coef_lin(17,17,9) = cs(9)	! 3A2D -> coag
	coef_lin(17,17,8) = cs(8)	! 2A2D -> coag
	coef_lin(17,17,7) = cs(7)	! 1A2D -> coag
	coef_lin(17,17,6) = cs(6)	! 3A1D -> coag
	coef_lin(17,17,5) = cs(5)	! 2A1D -> coag
	coef_lin(17,17,4) = cs(4)	! 1A1D -> coag
	coef_lin(17,17,3) = cs(3)	! 1D -> coag
	coef_lin(17,17,2) = cs(2)	! 2A -> coag
	coef_lin(17,17,1) = cs(1)	! 1A -> coag

end subroutine get_rate_coefs

!-----------------------------------------------------------

subroutine get_losses(cs)
	implicit none
	real(kind(1.d0)) :: cs(15)

	! coagulation sink

	! variable cs: the following values only give the size dependence of the sink,
	! and will be multiplied by the size-independent factor given as input

	cs(1) = 1.00000000000000d+000	! coagulation loss of 1A
	cs(2) = 6.90956439983888d-001	! coagulation loss of 2A
	cs(3) = 8.92789345574936d-001	! coagulation loss of 1D
	cs(4) = 6.50906449049502d-001	! coagulation loss of 1A1D
	cs(5) = 5.34476727643759d-001	! coagulation loss of 2A1D
	cs(6) = 4.62991019676722d-001	! coagulation loss of 3A1D
	cs(7) = 5.14714704412136d-001	! coagulation loss of 1A2D
	cs(8) = 4.49748002797798d-001	! coagulation loss of 2A2D
	cs(9) = 4.03871535975795d-001	! coagulation loss of 3A2D
	cs(10) = 3.69300136986970d-001	! coagulation loss of 4A2D
	cs(11) = 3.94846710007348d-001	! coagulation loss of 2A3D
	cs(12) = 3.62287963350499d-001	! coagulation loss of 3A3D
	cs(13) = 3.36418964945053d-001	! coagulation loss of 4A3D
	cs(14) = 3.31032373593705d-001	! coagulation loss of 3A4D
	cs(15) = 3.10756278903030d-001	! coagulation loss of 4A4D

end subroutine get_losses

!-----------------------------------------------------------

subroutine get_coll(K,temperature)
	implicit none
	integer, parameter :: nclust = 15
	real(kind(1.d0)) :: K(nclust,nclust), temperature

	! collision coefficients

	K = 0.d0
	K(1,1) = 2.00300435981201d-017*sqrt(temperature)	! 1A + 1A
	K(2,2) = 2.24829637648719d-017*sqrt(temperature)	! 2A + 2A
	K(3,1) = 2.71277008815771d-017*sqrt(temperature)	! 1D + 1A
	K(1,3) = K(3,1)
	K(3,2) = 3.15349593390535d-017*sqrt(temperature)	! 1D + 2A
	K(2,3) = K(3,2)
	K(4,1) = 2.44809288262135d-017*sqrt(temperature)	! 1A1D + 1A
	K(1,4) = K(4,1)
	K(4,2) = 2.54150641906904d-017*sqrt(temperature)	! 1A1D + 2A
	K(2,4) = K(4,2)
	K(4,3) = 3.39602459742711d-017*sqrt(temperature)	! 1A1D + 1D
	K(3,4) = K(4,3)
	K(4,4) = 2.83571813260610d-017*sqrt(temperature)	! 1A1D + 1A1D
	K(5,1) = 2.58125616100570d-017*sqrt(temperature)	! 2A1D + 1A
	K(1,5) = K(5,1)
	K(5,2) = 2.52957813819825d-017*sqrt(temperature)	! 2A1D + 2A
	K(2,5) = K(5,2)
	K(5,3) = 3.70774069121475d-017*sqrt(temperature)	! 2A1D + 1D
	K(3,5) = K(5,3)
	K(5,4) = 2.87378921679481d-017*sqrt(temperature)	! 2A1D + 1A1D
	K(4,5) = K(5,4)
	K(5,5) = 2.79470816379128d-017*sqrt(temperature)	! 2A1D + 2A1D
	K(6,2) = 2.60531351555787d-017*sqrt(temperature)	! 3A1D + 2A
	K(2,6) = K(6,2)
	K(6,3) = 4.02727477268666d-017*sqrt(temperature)	! 3A1D + 1D
	K(3,6) = K(6,3)
	K(6,4) = 2.99200115944640d-017*sqrt(temperature)	! 3A1D + 1A1D
	K(4,6) = K(6,4)
	K(6,5) = 2.83326354588116d-017*sqrt(temperature)	! 3A1D + 2A1D
	K(5,6) = K(6,5)
	K(6,6) = 2.81969208207520d-017*sqrt(temperature)	! 3A1D + 3A1D
	K(7,1) = 2.76109922240466d-017*sqrt(temperature)	! 1A2D + 1A
	K(1,7) = K(7,1)
	K(7,2) = 2.75406472349300d-017*sqrt(temperature)	! 1A2D + 2A
	K(2,7) = K(7,2)
	K(7,4) = 3.09758708170216d-017*sqrt(temperature)	! 1A2D + 1A1D
	K(4,7) = K(7,4)
	K(7,5) = 3.05649415958920d-017*sqrt(temperature)	! 1A2D + 2A1D
	K(5,7) = K(7,5)
	K(7,6) = 3.12749051699915d-017*sqrt(temperature)	! 1A2D + 3A1D
	K(6,7) = K(7,6)
	K(7,7) = 3.31633158250450d-017*sqrt(temperature)	! 1A2D + 1A2D
	K(8,1) = 2.87627441579623d-017*sqrt(temperature)	! 2A2D + 1A
	K(1,8) = K(8,1)
	K(8,2) = 2.74789135147876d-017*sqrt(temperature)	! 2A2D + 2A
	K(2,8) = K(8,2)
	K(8,3) = 4.16083096768037d-017*sqrt(temperature)	! 2A2D + 1D
	K(3,8) = K(8,3)
	K(8,4) = 3.13559696468930d-017*sqrt(temperature)	! 2A2D + 1A1D
	K(4,8) = K(8,4)
	K(8,5) = 2.99670606893001d-017*sqrt(temperature)	! 2A2D + 2A1D
	K(5,8) = K(8,5)
	K(8,6) = 3.00148303388612d-017*sqrt(temperature)	! 2A2D + 3A1D
	K(6,8) = K(8,6)
	K(8,7) = 3.29046951077028d-017*sqrt(temperature)	! 2A2D + 1A2D
	K(7,8) = K(8,7)
	K(8,8) = 3.18298598355307d-017*sqrt(temperature)	! 2A2D + 2A2D
	K(9,1) = 3.02705069238893d-017*sqrt(temperature)	! 3A2D + 1A
	K(1,9) = K(9,1)
	K(9,2) = 2.81058917527580d-017*sqrt(temperature)	! 3A2D + 2A
	K(2,9) = K(9,2)
	K(9,3) = 4.43959568621110d-017*sqrt(temperature)	! 3A2D + 1D
	K(3,9) = K(9,3)
	K(9,4) = 3.23641741717275d-017*sqrt(temperature)	! 3A2D + 1A1D
	K(4,9) = K(9,4)
	K(9,5) = 3.02671039469591d-017*sqrt(temperature)	! 3A2D + 2A1D
	K(5,9) = K(9,5)
	K(9,6) = 2.98505526994026d-017*sqrt(temperature)	! 3A2D + 3A1D
	K(6,9) = K(9,6)
	K(9,7) = 3.34974069815894d-017*sqrt(temperature)	! 3A2D + 1A2D
	K(7,9) = K(9,7)
	K(9,8) = 3.18341507188162d-017*sqrt(temperature)	! 3A2D + 2A2D
	K(8,9) = K(9,8)
	K(9,9) = 3.14250788514787d-017*sqrt(temperature)	! 3A2D + 3A2D
	K(10,2) = 2.89736207769365d-017*sqrt(temperature)	! 4A2D + 2A
	K(2,10) = K(10,2)
	K(10,3) = 4.71160370674057d-017*sqrt(temperature)	! 4A2D + 1D
	K(3,10) = K(10,3)
	K(10,4) = 3.35706769634303d-017*sqrt(temperature)	! 4A2D + 1A1D
	K(4,10) = K(10,4)
	K(10,5) = 3.09025016796110d-017*sqrt(temperature)	! 4A2D + 2A1D
	K(5,10) = K(10,5)
	K(10,6) = 3.01208396825086d-017*sqrt(temperature)	! 4A2D + 3A1D
	K(6,10) = K(10,6)
	K(10,7) = 3.43919941400074d-017*sqrt(temperature)	! 4A2D + 1A2D
	K(7,10) = K(10,7)
	K(10,8) = 3.22571933069290d-017*sqrt(temperature)	! 4A2D + 2A2D
	K(8,10) = K(10,8)
	K(10,9) = 3.15232262143403d-017*sqrt(temperature)	! 4A2D + 3A2D
	K(9,10) = K(10,9)
	K(10,10) = 3.13695384995609d-017*sqrt(temperature)	! 4A2D + 4A2D
	K(11,1) = 3.13194465823373d-017*sqrt(temperature)	! 2A3D + 1A
	K(1,11) = K(11,1)
	K(11,2) = 2.93359550362804d-017*sqrt(temperature)	! 2A3D + 2A
	K(2,11) = K(11,2)
	K(11,4) = 3.36010729643498d-017*sqrt(temperature)	! 2A3D + 1A1D
	K(4,11) = K(11,4)
	K(11,5) = 3.16706022115086d-017*sqrt(temperature)	! 2A3D + 2A1D
	K(5,11) = K(11,5)
	K(11,6) = 3.14101990762351d-017*sqrt(temperature)	! 2A3D + 3A1D
	K(6,11) = K(11,6)
	K(11,7) = 3.48944858931620d-017*sqrt(temperature)	! 2A3D + 1A2D
	K(7,11) = K(11,7)
	K(11,8) = 3.33884231970716d-017*sqrt(temperature)	! 2A3D + 2A2D
	K(8,11) = K(11,8)
	K(11,9) = 3.31235468974246d-017*sqrt(temperature)	! 2A3D + 3A2D
	K(9,11) = K(11,9)
	K(11,10) = 3.33527112669434d-017*sqrt(temperature)	! 2A3D + 4A2D
	K(10,11) = K(11,10)
	K(11,11) = 3.48145043488322d-017*sqrt(temperature)	! 2A3D + 2A3D
	K(12,1) = 3.26910501730354d-017*sqrt(temperature)	! 3A3D + 1A
	K(1,12) = K(12,1)
	K(12,2) = 2.99114746938185d-017*sqrt(temperature)	! 3A3D + 2A
	K(2,12) = K(12,2)
	K(12,3) = 4.80910847211470d-017*sqrt(temperature)	! 3A3D + 1D
	K(3,12) = K(12,3)
	K(12,4) = 3.45253049402141d-017*sqrt(temperature)	! 3A3D + 1A1D
	K(4,12) = K(12,4)
	K(12,5) = 3.19568882993891d-017*sqrt(temperature)	! 3A3D + 2A1D
	K(5,12) = K(12,5)
	K(12,6) = 3.12784270392478d-017*sqrt(temperature)	! 3A3D + 3A1D
	K(6,12) = K(12,6)
	K(12,7) = 3.54492926789563d-017*sqrt(temperature)	! 3A3D + 1A2D
	K(7,12) = K(12,7)
	K(12,8) = 3.34135780138988d-017*sqrt(temperature)	! 3A3D + 2A2D
	K(8,12) = K(12,8)
	K(12,9) = 3.27766124436665d-017*sqrt(temperature)	! 3A3D + 3A2D
	K(9,12) = K(12,9)
	K(12,10) = 3.27136209916630d-017*sqrt(temperature)	! 3A3D + 4A2D
	K(10,12) = K(12,10)
	K(12,11) = 3.46023050442563d-017*sqrt(temperature)	! 3A3D + 2A3D
	K(11,12) = K(12,11)
	K(12,12) = 3.40551869990935d-017*sqrt(temperature)	! 3A3D + 3A3D
	K(13,2) = 3.06994479028858d-017*sqrt(temperature)	! 4A3D + 2A
	K(2,13) = K(13,2)
	K(13,3) = 5.05871386086251d-017*sqrt(temperature)	! 4A3D + 1D
	K(3,13) = K(13,3)
	K(13,4) = 3.56268165927125d-017*sqrt(temperature)	! 4A3D + 1A1D
	K(4,13) = K(13,4)
	K(13,5) = 3.25315618247093d-017*sqrt(temperature)	! 4A3D + 2A1D
	K(5,13) = K(13,5)
	K(13,6) = 3.15159917438775d-017*sqrt(temperature)	! 4A3D + 3A1D
	K(6,13) = K(13,6)
	K(13,7) = 3.62642762902657d-017*sqrt(temperature)	! 4A3D + 1A2D
	K(7,13) = K(13,7)
	K(13,8) = 3.37939933755832d-017*sqrt(temperature)	! 4A3D + 2A2D
	K(8,13) = K(13,8)
	K(13,9) = 3.28572029416204d-017*sqrt(temperature)	! 4A3D + 3A2D
	K(9,13) = K(13,9)
	K(13,10) = 3.25617038244942d-017*sqrt(temperature)	! 4A3D + 4A2D
	K(10,13) = K(13,10)
	K(13,11) = 3.48057086067077d-017*sqrt(temperature)	! 4A3D + 2A3D
	K(11,13) = K(13,11)
	K(13,12) = 3.39890106915372d-017*sqrt(temperature)	! 4A3D + 3A3D
	K(12,13) = K(13,12)
	K(13,13) = 3.37082033742580d-017*sqrt(temperature)	! 4A3D + 4A3D
	K(14,1) = 3.49043670857906d-017*sqrt(temperature)	! 3A4D + 1A
	K(1,14) = K(14,1)
	K(14,2) = 3.15476588186036d-017*sqrt(temperature)	! 3A4D + 2A
	K(2,14) = K(14,2)
	K(14,4) = 3.64901348708458d-017*sqrt(temperature)	! 3A4D + 1A1D
	K(4,14) = K(14,4)
	K(14,5) = 3.34815796368918d-017*sqrt(temperature)	! 3A4D + 2A1D
	K(5,14) = K(14,5)
	K(14,6) = 3.25574226703014d-017*sqrt(temperature)	! 3A4D + 3A1D
	K(6,14) = K(14,6)
	K(14,7) = 3.72165162582506d-017*sqrt(temperature)	! 3A4D + 1A2D
	K(7,14) = K(14,7)
	K(14,8) = 3.48334426724349d-017*sqrt(temperature)	! 3A4D + 2A2D
	K(8,14) = K(14,8)
	K(14,9) = 3.39830001973208d-017*sqrt(temperature)	! 3A4D + 3A2D
	K(9,14) = K(14,9)
	K(14,10) = 3.37686092929410d-017*sqrt(temperature)	! 3A4D + 4A2D
	K(10,14) = K(14,10)
	K(14,11) = 3.59271247979772d-017*sqrt(temperature)	! 3A4D + 2A3D
	K(11,14) = K(14,11)
	K(14,12) = 3.51926015491360d-017*sqrt(temperature)	! 3A4D + 3A3D
	K(12,14) = K(14,12)
	K(14,13) = 3.49886518498441d-017*sqrt(temperature)	! 3A4D + 4A3D
	K(13,14) = K(14,13)
	K(14,14) = 3.62650661050272d-017*sqrt(temperature)	! 3A4D + 3A4D
	K(15,1) = 3.62596102037558d-017*sqrt(temperature)	! 4A4D + 1A
	K(1,15) = K(15,1)
	K(15,2) = 3.22850887225448d-017*sqrt(temperature)	! 4A4D + 2A
	K(2,15) = K(15,2)
	K(15,3) = 5.38075352932129d-017*sqrt(temperature)	! 4A4D + 1D
	K(3,15) = K(15,3)
	K(15,4) = 3.75206694350681d-017*sqrt(temperature)	! 4A4D + 1A1D
	K(4,15) = K(15,4)
	K(15,5) = 3.40225187661203d-017*sqrt(temperature)	! 4A4D + 2A1D
	K(5,15) = K(15,5)
	K(15,6) = 3.27852539241479d-017*sqrt(temperature)	! 4A4D + 3A1D
	K(6,15) = K(15,6)
	K(15,7) = 3.79823162049186d-017*sqrt(temperature)	! 4A4D + 1A2D
	K(7,15) = K(15,7)
	K(15,8) = 3.51958859165781d-017*sqrt(temperature)	! 4A4D + 2A2D
	K(8,15) = K(15,8)
	K(15,9) = 3.40671002718251d-017*sqrt(temperature)	! 4A4D + 3A2D
	K(9,15) = K(15,9)
	K(15,10) = 3.36368883231231d-017*sqrt(temperature)	! 4A4D + 4A2D
	K(10,15) = K(15,10)
	K(15,11) = 3.61271601610681d-017*sqrt(temperature)	! 4A4D + 2A3D
	K(11,15) = K(15,11)
	K(15,12) = 3.51424119700538d-017*sqrt(temperature)	! 4A4D + 3A3D
	K(12,15) = K(15,12)
	K(15,13) = 3.47392500725570d-017*sqrt(temperature)	! 4A4D + 4A3D
	K(13,15) = K(15,13)
	K(15,14) = 3.60888862121434d-017*sqrt(temperature)	! 4A4D + 3A4D
	K(14,15) = K(15,14)
	K(15,15) = 3.57278096683900d-017*sqrt(temperature)	! 4A4D + 4A4D

end subroutine get_coll

!-----------------------------------------------------------

subroutine get_evap(E,K,temperature)
	implicit none
	real(kind(1.d0)) :: E(15,15), K(15,15), temperature

	! evaporation coefficients

	E = 0.d0
	E(1,1) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((-17.8487481487d0/temperature-(-33.418352d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(1,1)	! 2A -> 1A + 1A
	E(3,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-24.6465441592d0/temperature-(-31.008553d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(3,1)	! 1A1D -> 1D + 1A
	E(1,3) = E(3,1)
	E(4,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(4,1)	! 2A1D -> 1A1D + 1A
	E(1,4) = E(4,1)
	E(3,2) = 7.33893243358348d+027/temperature*exp(&
			 &((-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-0.d0&
			 &-(-17.8487481487d0/temperature-(-33.418352d0)/1.d3))&
			 &*5.03218937158374d+002)*K(3,2)	! 2A1D -> 1D + 2A
	E(2,3) = E(3,2)
	E(5,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-81.4604213712d0/temperature-(-112.898920d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(5,1)	! 3A1D -> 2A1D + 1A
	E(1,5) = E(5,1)
	E(4,2) = 7.33893243358348d+027/temperature*exp(&
			 &((-81.4604213712d0/temperature-(-112.898920d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3)&
			 &-(-17.8487481487d0/temperature-(-33.418352d0)/1.d3))&
			 &*5.03218937158374d+002)*K(4,2)	! 3A1D -> 1A1D + 2A
	E(2,4) = E(4,2)
	E(4,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(4,3)	! 1A2D -> 1A1D + 1D
	E(3,4) = E(4,3)
	E(7,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(7,1)	! 2A2D -> 1A2D + 1A
	E(1,7) = E(7,1)
	E(5,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(5,3)	! 2A2D -> 2A1D + 1D
	E(3,5) = E(5,3)
	E(4,4) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(4,4)	! 2A2D -> 1A1D + 1A1D
	E(8,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(8,1)	! 3A2D -> 2A2D + 1A
	E(1,8) = E(8,1)
	E(7,2) = 7.33893243358348d+027/temperature*exp(&
			 &((-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-(-17.8487481487d0/temperature-(-33.418352d0)/1.d3))&
			 &*5.03218937158374d+002)*K(7,2)	! 3A2D -> 1A2D + 2A
	E(2,7) = E(7,2)
	E(6,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-81.4604213712d0/temperature-(-112.898920d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(6,3)	! 3A2D -> 3A1D + 1D
	E(3,6) = E(6,3)
	E(5,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(5,4)	! 3A2D -> 2A1D + 1A1D
	E(4,5) = E(5,4)
	E(9,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-140.75630866991d0/temperature-(-194.432572d0)/1.d3)&
			 &-(-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(9,1)	! 4A2D -> 3A2D + 1A
	E(1,9) = E(9,1)
	E(8,2) = 7.33893243358348d+027/temperature*exp(&
			 &((-140.75630866991d0/temperature-(-194.432572d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-17.8487481487d0/temperature-(-33.418352d0)/1.d3))&
			 &*5.03218937158374d+002)*K(8,2)	! 4A2D -> 2A2D + 2A
	E(2,8) = E(8,2)
	E(6,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-140.75630866991d0/temperature-(-194.432572d0)/1.d3)&
			 &-(-81.4604213712d0/temperature-(-112.898920d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(6,4)	! 4A2D -> 3A1D + 1A1D
	E(4,6) = E(6,4)
	E(5,5) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((-140.75630866991d0/temperature-(-194.432572d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3))&
			 &*5.03218937158374d+002)*K(5,5)	! 4A2D -> 2A1D + 2A1D
	E(8,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(8,3)	! 2A3D -> 2A2D + 1D
	E(3,8) = E(8,3)
	E(7,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(7,4)	! 2A3D -> 1A2D + 1A1D
	E(4,7) = E(7,4)
	E(11,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-(-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(11,1)	! 3A3D -> 2A3D + 1A
	E(1,11) = E(11,1)
	E(9,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-(-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(9,3)	! 3A3D -> 3A2D + 1D
	E(3,9) = E(9,3)
	E(8,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(8,4)	! 3A3D -> 2A2D + 1A1D
	E(4,8) = E(8,4)
	E(7,5) = 7.33893243358348d+027/temperature*exp(&
			 &((-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3))&
			 &*5.03218937158374d+002)*K(7,5)	! 3A3D -> 1A2D + 2A1D
	E(5,7) = E(7,5)
	E(12,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(12,1)	! 4A3D -> 3A3D + 1A
	E(1,12) = E(12,1)
	E(11,2) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-(-17.8487481487d0/temperature-(-33.418352d0)/1.d3))&
			 &*5.03218937158374d+002)*K(11,2)	! 4A3D -> 2A3D + 2A
	E(2,11) = E(11,2)
	E(10,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-140.75630866991d0/temperature-(-194.432572d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(10,3)	! 4A3D -> 4A2D + 1D
	E(3,10) = E(10,3)
	E(9,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(9,4)	! 4A3D -> 3A2D + 1A1D
	E(4,9) = E(9,4)
	E(8,5) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3))&
			 &*5.03218937158374d+002)*K(8,5)	! 4A3D -> 2A2D + 2A1D
	E(5,8) = E(8,5)
	E(7,6) = 7.33893243358348d+027/temperature*exp(&
			 &((-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3)&
			 &-(-81.4604213712d0/temperature-(-112.898920d0)/1.d3))&
			 &*5.03218937158374d+002)*K(7,6)	! 4A3D -> 1A2D + 3A1D
	E(6,7) = E(7,6)
	E(12,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-174.6951416018d0/temperature-(-231.377403d0)/1.d3)&
			 &-(-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(12,3)	! 3A4D -> 3A3D + 1D
	E(3,12) = E(12,3)
	E(11,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-174.6951416018d0/temperature-(-231.377403d0)/1.d3)&
			 &-(-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(11,4)	! 3A4D -> 2A3D + 1A1D
	E(4,11) = E(11,4)
	E(8,7) = 7.33893243358348d+027/temperature*exp(&
			 &((-174.6951416018d0/temperature-(-231.377403d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3))&
			 &*5.03218937158374d+002)*K(8,7)	! 3A4D -> 2A2D + 1A2D
	E(7,8) = E(8,7)
	E(14,1) = 7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-174.6951416018d0/temperature-(-231.377403d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(14,1)	! 4A4D -> 3A4D + 1A
	E(1,14) = E(14,1)
	E(13,3) = 7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-173.65103607928d0/temperature-(-229.951283d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+002)*K(13,3)	! 4A4D -> 4A3D + 1D
	E(3,13) = E(13,3)
	E(12,4) = 7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-146.5927558198d0/temperature-(-190.007277d0)/1.d3)&
			 &-(-24.6465441592d0/temperature-(-31.008553d0)/1.d3))&
			 &*5.03218937158374d+002)*K(12,4)	! 4A4D -> 3A3D + 1A1D
	E(4,12) = E(12,4)
	E(11,5) = 7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-108.2392632326d0/temperature-(-152.755066d0)/1.d3)&
			 &-(-57.0595364880d0/temperature-(-74.406210d0)/1.d3))&
			 &*5.03218937158374d+002)*K(11,5)	! 4A4D -> 2A3D + 2A1D
	E(5,11) = E(11,5)
	E(9,7) = 7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-113.5835639644d0/temperature-(-155.375202d0)/1.d3)&
			 &-(-40.0759739702d0/temperature-(-66.369881d0)/1.d3))&
			 &*5.03218937158374d+002)*K(9,7)	! 4A4D -> 3A2D + 1A2D
	E(7,9) = E(9,7)
	E(8,8) = 0.5*7.33893243358348d+027/temperature*exp(&
			 &((-204.2956058491d0/temperature-(-269.137936d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3)&
			 &-(-87.5692739780d0/temperature-(-111.732097d0)/1.d3))&
			 &*5.03218937158374d+002)*K(8,8)	! 4A4D -> 2A2D + 2A2D

end subroutine get_evap
