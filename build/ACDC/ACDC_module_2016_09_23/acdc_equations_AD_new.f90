!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!
! Cluster 1: 1A
! Cluster 2: 2A
! Cluster 3: 3A
! Cluster 4: 4A
! Cluster 5: 1D
! Cluster 6: 1A1D
! Cluster 7: 2A1D
! Cluster 8: 3A1D
! Cluster 9: 4A1D
! Cluster 10: 2D
! Cluster 11: 1A2D
! Cluster 12: 2A2D
! Cluster 13: 3A2D
! Cluster 14: 4A2D
! Cluster 15: 3D
! Cluster 16: 1A3D
! Cluster 17: 2A3D
! Cluster 18: 3A3D
! Cluster 19: 4A3D
! Cluster 20: 4D
! Cluster 21: 1A4D
! Cluster 22: 2A4D
! Cluster 23: 3A4D
! Cluster 24: 4A4D
! 25 is for source fluxes
! 26 is for coagulation losses 
! 27 is for wall losses 
! 28 is for dilution losses 
! 29 is for recombination of positive and negative charger ions with each other
! 30 is for collisions that lead succesfully out of the system as neutrals
! 31 is for collisions that lead out of the system, but the resulting cluster is brought back
!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!

! differential equations: f = dc/dt
subroutine fevalD(neqn,t,c,f,coef,ipar)
	implicit none
	integer, parameter :: nclust = 24
	integer :: neqn, ipar(4), i, j, k, n
	real(kind(1.d0)) :: t, c(neqn), f(neqn), coef(2)
	logical, save :: isconst(31) = .false.
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,31)=0.d0,coef_lin(31,31,nclust)=0.d0,source(31)=0.d0
	integer, save :: ind_quad_loss(31,0:50)=0,ind_quad_form(31,0:156)=0
	integer, save :: ind_quad_loss_extra(31,0:60)=0,ind_quad_form_extra(31,0:219)=0
	integer, save :: ind_lin_loss(31,0:26)=0,ind_lin_form(31,0:48)=0,fitted(24,0:24)=0
	real(kind(1.d0)), save :: n_quad_form_extra(31,0:73)=0.d0
	real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)


	! the parameters are read at the very beginning and eg. if source terms or collision rates have changed
	if (ipar(1) .eq. 0) then
		! after this the same values are used until some other routine tells otherwise
		ipar(1) = 1
		call initialize_parametersD(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
		&	ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)
		n_quad_form_extra(:,1:73) = real(ind_quad_form_extra(:,3:219:3),kind=kind(1.d0))
	end if

	f = 0.d0
	do i = 1, neqn ! loop over all the clusters + other fluxes
		if (.not. isconst(i)) then
			! first calculate coefficients for all loss terms
			do n = 1, 2*ind_quad_loss(i,0)-1, 2 ! loop over all collisions removing this cluster
				f(i) = f(i)-coef_quad(i,ind_quad_loss(i,n),ind_quad_loss(i,n+1))*c(ind_quad_loss(i,n))
			end do
			do n = 1, 2*ind_lin_loss(i,0)-1, 2 ! loop over all evaporations + wall losses etc. removing this cluster
				f(i) = f(i)-coef_lin(ind_lin_loss(i,n),ind_lin_loss(i,n+1),i)
			end do
			f(i) = f(i)*c(i) ! multiplying loss coefficients with current concentration
			! then add all terms that form this cluster
			do n = 1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions forming this cluster
				f(i) = f(i)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&
				&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))
			end do
			do n = 1, 3*ind_quad_form_extra(i,0)-2, 3 ! loop over all collisions forming this cluster as an extra product
				f(i) = f(i)+coef_quad(ind_quad_form_extra(i,n),ind_quad_form_extra(i,n+1),i)*&
				&      n_quad_form_extra(i,(n+2)/3)*c(ind_quad_form_extra(i,n))*c(ind_quad_form_extra(i,n+1))
			end do
			do n = 1, 2*ind_lin_form(i,0)-1, 2 ! loop over all evaporations forming this cluster
				f(i) = f(i)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))
			end do
			! finally, add possible external sources
			f(i) = f(i) + source(i)
		end if
	end do
	do n = 1, fitted(1,0) ! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)
		f(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations
		do i = 1, fitted(n+1,1) ! loop over the clusters used for the fitting
			f(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i+1)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)
		end do
	end do

end subroutine fevalD

!-----------------------------------------------------------

! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)
! not using this, since the solver works faster using a numerical jacobian...
subroutine jevalD(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)
	implicit none
	integer :: ldim,neqn,ierr,ipar(4),i,j,k,n
	real(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn),coef(2),ml,mu

end subroutine jevalD

!-----------------------------------------------------------

subroutine formationD(neqn,c,j_tot,j_by_charge,j_by_cluster,j_all,ipar,coef)
	implicit none
	integer, parameter :: nclust = 24
	integer :: neqn, i, n, ipar(4)
	real(kind(1.d0)) :: c(neqn), j_add, j_tot,j_by_charge(4),j_by_cluster(neqn),j_all(neqn,4),coef(2)
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,31)=0.d0,coef_lin(31,31,nclust)=0.d0,source(31)=0.d0
	integer, save :: charge(nclust)=0,ind_quad_loss(31,0:50)=0,ind_quad_form(31,0:156)=0
	integer, save :: ind_quad_loss_extra(31,0:60)=0,ind_quad_form_extra(31,0:219)=0
	integer, save :: ind_lin_loss(31,0:26)=0,ind_lin_form(31,0:48)=0,fitted(nclust,0:nclust)=0
	logical, save :: isconst(31)=.false.

	if (ipar(3) .eq. 0) then
		ipar(3) = 1
		call initialize_parametersD(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)
		! cluster charges
		charge = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	end if
	j_tot = 0.d0			! total formation rate
	j_by_charge = 0.d0		! contributions of neutrals, negatives, positives and recombinations
	j_by_cluster = 0.d0		! contributions of different clusters
	j_all = 0.d0			! contributions of different clusters and different charging states
	do n = 1, 2*ind_quad_form(30,0)-1, 2 ! loop over all collisions leading out of the system
		j_add = coef_quad(ind_quad_form(30,n),ind_quad_form(30,n+1),30)*c(ind_quad_form(30,n))*c(ind_quad_form(30,n+1))
		j_tot = j_tot + j_add
		j_by_charge(1) = j_by_charge(1) + j_add
		j_by_cluster(ind_quad_form(30,n)) = j_by_cluster(ind_quad_form(30,n)) + j_add
		j_by_cluster(ind_quad_form(30,n+1)) = j_by_cluster(ind_quad_form(30,n+1)) + j_add
		j_all(ind_quad_form(30,n),1) = j_all(ind_quad_form(30,n),1) + j_add
		j_all(ind_quad_form(30,n+1),1) = j_all(ind_quad_form(30,n+1),1) + j_add
	end do

end subroutine formationD

!-----------------------------------------------------------

subroutine initialize_parametersD(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	source,isconst,fitted,coef,ipar)

	use monomer_settingsD, only : sources_and_constants

	implicit none
	integer, parameter :: neqn=31, nclust=24
	logical :: isconst(31)
	real(kind(1.d0)) :: coef_quad(nclust,nclust,31),coef_lin(31,31,nclust)
	real(kind(1.d0)) :: source(31)
	integer :: ind_quad_loss(31,0:50),ind_quad_form(31,0:156),ind_quad_loss_extra(31,0:60),ind_quad_form_extra(31,0:219)
	integer :: ind_lin_loss(31,0:26),ind_lin_form(31,0:48)
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
	ind_quad_loss(2,0:50) = (/ 25, 1,3, 2,4, 2,4, 3,4, 4,4, 5,7, 6,8, 7,9, 8,9&
			&, 9,9, 10,12, 11,13, 12,14, 13,14, 14,14, 15,17, 16,18, 17,19, 18,19&
			&, 19,19, 20,22, 21,23, 22,24, 23,30, 24,30 /)
	ind_quad_loss_extra(2,0:24) = (/ 8, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 13,1,1, 14,1,2, 18,1,1, 19,1,2 /)
	ind_quad_form(2,0:2) = (/ 1, 1,1 /)
	ind_lin_loss(2,0:4) = (/ 2, 1,1, 26,26 /)
	ind_lin_form(2,0:30) = (/ 15, 1,3, 2,4, 2,4, 5,7, 6,8, 7,9, 10,12, 11,13, 12,14&
			&, 15,17, 16,18, 17,19, 20,22, 21,23, 22,24 /)

	! Cluster 3: 3A
	ind_quad_loss(3,0:50) = (/ 25, 1,4, 2,4, 3,4, 3,4, 4,4, 5,8, 6,9, 7,9, 8,9&
			&, 9,9, 10,13, 11,14, 12,14, 13,14, 14,14, 15,18, 16,19, 17,19, 18,19&
			&, 19,19, 20,23, 21,24, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(3,0:39) = (/ 13, 2,1,1, 3,1,2, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 9,1,3, 12,1,1, 13,1,2&
			&, 14,1,3, 17,1,1, 18,1,2, 19,1,3 /)
	ind_quad_form(3,0:2) = (/ 1, 1,2 /)
	ind_lin_loss(3,0:4) = (/ 2, 1,2, 26,26 /)
	ind_lin_form(3,0:18) = (/ 9, 1,4, 5,8, 6,9, 10,13, 11,14, 15,18, 16,19, 20,23, 21,24 /)

	! Cluster 4: 4A
	ind_quad_loss(4,0:48) = (/ 24, 2,4, 3,4, 4,4, 4,4, 5,9, 6,9, 7,9, 8,9, 9,9&
			&, 10,14, 11,14, 12,14, 13,14, 14,14, 15,19, 16,19, 17,19, 18,19, 19,19&
			&, 20,24, 21,30, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(4,0:48) = (/ 16, 2,1,2, 3,1,3, 4,1,4, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4, 11,1,1&
			&, 12,1,2, 13,1,3, 14,1,4, 16,1,1, 17,1,2, 18,1,3, 19,1,4 /)
	ind_quad_form(4,0:14) = (/ 7, 1,3, 2,2, 2,3, 2,4, 3,3, 3,4, 4,4 /)
	ind_lin_loss(4,0:6) = (/ 3, 1,3, 2,2, 26,26 /)
	ind_lin_form(4,0:8) = (/ 4, 5,9, 10,14, 15,19, 20,24 /)

	! Cluster 5: 1D
	! Derivative not used and concentration set to be constant
	isconst(5) = .true.

	! Cluster 6: 1A1D
	ind_quad_loss(6,0:50) = (/ 25, 1,7, 2,8, 3,9, 4,9, 5,11, 6,12, 6,12, 7,13, 8,14&
			&, 9,14, 10,16, 11,17, 12,18, 13,19, 14,19, 15,21, 16,22, 17,23, 18,24&
			&, 19,30, 20,21, 21,22, 22,23, 23,24, 24,30 /)
	ind_quad_loss_extra(6,0:21) = (/ 7, 4,1,1, 9,1,1, 14,1,1, 20,5,1, 21,5,1, 22,5,1, 23,5,1 /)
	ind_quad_form(6,0:2) = (/ 1, 1,5 /)
	ind_lin_loss(6,0:4) = (/ 2, 1,5, 26,26 /)
	ind_lin_form(6,0:32) = (/ 16, 1,7, 2,8, 3,9, 5,11, 6,12, 6,12, 7,13, 8,14, 10,16&
			&, 11,17, 12,18, 13,19, 15,21, 16,22, 17,23, 18,24 /)

	! Cluster 7: 2A1D
	ind_quad_loss(7,0:50) = (/ 25, 1,8, 2,9, 3,9, 4,9, 5,12, 6,13, 7,14, 7,14, 8,14&
			&, 9,14, 10,17, 11,18, 12,19, 13,19, 14,19, 15,22, 16,23, 17,24, 18,30&
			&, 19,30, 20,22, 21,23, 22,24, 23,30, 24,30 /)
	ind_quad_loss_extra(7,0:27) = (/ 9, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 13,1,1, 14,1,2, 20,5,1, 21,5,1, 22,5,1 /)
	ind_quad_form(7,0:4) = (/ 2, 1,6, 2,5 /)
	ind_lin_loss(7,0:6) = (/ 3, 1,6, 2,5, 26,26 /)
	ind_lin_form(7,0:24) = (/ 12, 1,8, 2,9, 5,12, 6,13, 7,14, 7,14, 10,17, 11,18, 12,19&
			&, 15,22, 16,23, 17,24 /)

	! Cluster 8: 3A1D
	ind_quad_loss(8,0:50) = (/ 25, 1,9, 2,9, 3,9, 4,9, 5,13, 6,14, 7,14, 8,14, 8,14&
			&, 9,14, 10,18, 11,19, 12,19, 13,19, 14,19, 15,23, 16,24, 17,30, 18,30&
			&, 19,30, 20,23, 21,24, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(8,0:36) = (/ 12, 2,1,1, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 8,1,2, 9,1,3, 12,1,1, 13,1,2&
			&, 14,1,3, 20,5,1, 21,5,1 /)
	ind_quad_form(8,0:6) = (/ 3, 1,7, 2,6, 3,5 /)
	ind_lin_loss(8,0:8) = (/ 4, 1,7, 2,6, 3,5, 26,26 /)
	ind_lin_form(8,0:14) = (/ 7, 1,9, 5,13, 6,14, 10,18, 11,19, 15,23, 16,24 /)

	! Cluster 9: 4A1D
	ind_quad_loss(9,0:48) = (/ 24, 2,9, 3,9, 4,9, 5,14, 6,14, 7,14, 8,14, 9,14, 9,14&
			&, 10,19, 11,19, 12,19, 13,19, 14,19, 15,24, 16,30, 17,30, 18,30, 19,30&
			&, 20,24, 21,30, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(9,0:39) = (/ 13, 2,1,2, 3,1,3, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4, 9,1,4, 11,1,1&
			&, 12,1,2, 13,1,3, 14,1,4, 20,5,1 /)
	ind_quad_form(9,0:26) = (/ 13, 1,8, 2,7, 2,8, 2,9, 3,6, 3,7, 3,8, 3,9, 4,5&
			&, 4,6, 4,7, 4,8, 4,9 /)
	ind_lin_loss(9,0:10) = (/ 5, 1,8, 2,7, 3,6, 4,5, 26,26 /)
	ind_lin_form(9,0:6) = (/ 3, 5,14, 10,19, 15,24 /)

	! Cluster 10: 2D
	ind_quad_loss(10,0:50) = (/ 25, 1,11, 2,12, 3,13, 4,14, 5,15, 6,16, 7,17, 8,18, 9,19&
			&, 10,20, 10,20, 11,21, 12,22, 13,23, 14,24, 15,20, 16,21, 17,22, 18,23&
			&, 19,24, 20,20, 21,21, 22,22, 23,23, 24,24 /)
	ind_quad_loss_extra(10,0:30) = (/ 10, 15,5,1, 16,5,1, 17,5,1, 18,5,1, 19,5,1, 20,5,2, 21,5,2, 22,5,2, 23,5,2&
			&, 24,5,2 /)
	ind_quad_form(10,0:2) = (/ 1, 5,5 /)
	ind_lin_loss(10,0:4) = (/ 2, 5,5, 26,26 /)
	ind_lin_form(10,0:30) = (/ 15, 1,11, 2,12, 3,13, 4,14, 5,15, 6,16, 7,17, 8,18, 9,19&
			&, 10,20, 10,20, 11,21, 12,22, 13,23, 14,24 /)

	! Cluster 11: 1A2D
	ind_quad_loss(11,0:50) = (/ 25, 1,12, 2,13, 3,14, 4,14, 5,16, 6,17, 7,18, 8,19, 9,19&
			&, 10,21, 11,22, 11,22, 12,23, 13,24, 14,30, 15,21, 16,22, 17,23, 18,24&
			&, 19,30, 20,21, 21,22, 22,23, 23,24, 24,30 /)
	ind_quad_loss_extra(11,0:30) = (/ 10, 4,1,1, 9,1,1, 15,5,1, 16,5,1, 17,5,1, 18,5,1, 20,5,2, 21,5,2, 22,5,2&
			&, 23,5,2 /)
	ind_quad_form(11,0:4) = (/ 2, 1,10, 5,6 /)
	ind_lin_loss(11,0:6) = (/ 3, 1,10, 5,6, 26,26 /)
	ind_lin_form(11,0:24) = (/ 12, 1,12, 2,13, 3,14, 5,16, 6,17, 7,18, 8,19, 10,21, 11,22&
			&, 11,22, 12,23, 13,24 /)

	! Cluster 12: 2A2D
	ind_quad_loss(12,0:50) = (/ 25, 1,13, 2,14, 3,14, 4,14, 5,17, 6,18, 7,19, 8,19, 9,19&
			&, 10,22, 11,23, 12,24, 12,24, 13,30, 14,30, 15,22, 16,23, 17,24, 18,30&
			&, 19,30, 20,22, 21,23, 22,24, 23,30, 24,30 /)
	ind_quad_loss_extra(12,0:30) = (/ 10, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 15,5,1, 16,5,1, 17,5,1, 20,5,2, 21,5,2&
			&, 22,5,2 /)
	ind_quad_form(12,0:8) = (/ 4, 1,11, 2,10, 5,7, 6,6 /)
	ind_lin_loss(12,0:10) = (/ 5, 1,11, 2,10, 5,7, 6,6, 26,26 /)
	ind_lin_form(12,0:18) = (/ 9, 1,13, 2,14, 5,17, 6,18, 7,19, 10,22, 11,23, 12,24, 12,24 /)

	! Cluster 13: 3A2D
	ind_quad_loss(13,0:50) = (/ 25, 1,14, 2,14, 3,14, 4,14, 5,18, 6,19, 7,19, 8,19, 9,19&
			&, 10,23, 11,24, 12,30, 13,30, 13,30, 14,30, 15,23, 16,24, 17,30, 18,30&
			&, 19,30, 20,23, 21,24, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(13,0:30) = (/ 10, 2,1,1, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 9,1,3, 15,5,1, 16,5,1, 20,5,2&
			&, 21,5,2 /)
	ind_quad_form(13,0:10) = (/ 5, 1,12, 2,11, 3,10, 5,8, 6,7 /)
	ind_lin_loss(13,0:12) = (/ 6, 1,12, 2,11, 3,10, 5,8, 6,7, 26,26 /)
	ind_lin_form(13,0:10) = (/ 5, 1,14, 5,18, 6,19, 10,23, 11,24 /)

	! Cluster 14: 4A2D
	ind_quad_loss(14,0:48) = (/ 24, 2,14, 3,14, 4,14, 5,19, 6,19, 7,19, 8,19, 9,19, 10,24&
			&, 11,30, 12,30, 13,30, 14,30, 14,30, 15,24, 16,30, 17,30, 18,30, 19,30&
			&, 20,24, 21,30, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(14,0:27) = (/ 9, 2,1,2, 3,1,3, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4, 15,5,1, 20,5,2 /)
	ind_quad_form(14,0:44) = (/ 22, 1,13, 2,12, 2,13, 2,14, 3,11, 3,12, 3,13, 3,14, 4,10&
			&, 4,11, 4,12, 4,13, 4,14, 5,9, 6,8, 6,9, 7,7, 7,8, 7,9&
			&, 8,8, 8,9, 9,9 /)
	ind_lin_loss(14,0:16) = (/ 8, 1,13, 2,12, 3,11, 4,10, 5,9, 6,8, 7,7, 26,26 /)
	ind_lin_form(14,0:4) = (/ 2, 5,19, 10,24 /)

	! Cluster 15: 3D
	ind_quad_loss(15,0:50) = (/ 25, 1,16, 2,17, 3,18, 4,19, 5,20, 6,21, 7,22, 8,23, 9,24&
			&, 10,20, 11,21, 12,22, 13,23, 14,24, 15,20, 15,20, 16,21, 17,22, 18,23&
			&, 19,24, 20,20, 21,21, 22,22, 23,23, 24,24 /)
	ind_quad_loss_extra(15,0:48) = (/ 16, 10,5,1, 11,5,1, 12,5,1, 13,5,1, 14,5,1, 15,5,2, 15,5,2, 16,5,2, 17,5,2&
			&, 18,5,2, 19,5,2, 20,5,3, 21,5,3, 22,5,3, 23,5,3, 24,5,3 /)
	ind_quad_form(15,0:2) = (/ 1, 5,10 /)
	ind_lin_loss(15,0:4) = (/ 2, 5,10, 26,26 /)
	ind_lin_form(15,0:18) = (/ 9, 1,16, 2,17, 3,18, 4,19, 5,20, 6,21, 7,22, 8,23, 9,24 /)

	! Cluster 16: 1A3D
	ind_quad_loss(16,0:50) = (/ 25, 1,17, 2,18, 3,19, 4,19, 5,21, 6,22, 7,23, 8,24, 9,30&
			&, 10,21, 11,22, 12,23, 13,24, 14,30, 15,21, 16,22, 16,22, 17,23, 18,24&
			&, 19,30, 20,21, 21,22, 22,23, 23,24, 24,30 /)
	ind_quad_loss_extra(16,0:42) = (/ 14, 4,1,1, 10,5,1, 11,5,1, 12,5,1, 13,5,1, 15,5,2, 16,5,2, 16,5,2, 17,5,2&
			&, 18,5,2, 20,5,3, 21,5,3, 22,5,3, 23,5,3 /)
	ind_quad_form(16,0:6) = (/ 3, 1,15, 5,11, 6,10 /)
	ind_lin_loss(16,0:8) = (/ 4, 1,15, 5,11, 6,10, 26,26 /)
	ind_lin_form(16,0:14) = (/ 7, 1,17, 2,18, 3,19, 5,21, 6,22, 7,23, 8,24 /)

	! Cluster 17: 2A3D
	ind_quad_loss(17,0:50) = (/ 25, 1,18, 2,19, 3,19, 4,19, 5,22, 6,23, 7,24, 8,30, 9,30&
			&, 10,22, 11,23, 12,24, 13,30, 14,30, 15,22, 16,23, 17,24, 17,24, 18,30&
			&, 19,30, 20,22, 21,23, 22,24, 23,30, 24,30 /)
	ind_quad_loss_extra(17,0:36) = (/ 12, 3,1,1, 4,1,2, 10,5,1, 11,5,1, 12,5,1, 15,5,2, 16,5,2, 17,5,2, 17,5,2&
			&, 20,5,3, 21,5,3, 22,5,3 /)
	ind_quad_form(17,0:10) = (/ 5, 1,16, 2,15, 5,12, 6,11, 7,10 /)
	ind_lin_loss(17,0:12) = (/ 6, 1,16, 2,15, 5,12, 6,11, 7,10, 26,26 /)
	ind_lin_form(17,0:10) = (/ 5, 1,18, 2,19, 5,22, 6,23, 7,24 /)

	! Cluster 18: 3A3D
	ind_quad_loss(18,0:50) = (/ 25, 1,19, 2,19, 3,19, 4,19, 5,23, 6,24, 7,30, 8,30, 9,30&
			&, 10,23, 11,24, 12,30, 13,30, 14,30, 15,23, 16,24, 17,30, 18,30, 18,30&
			&, 19,30, 20,23, 21,24, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(18,0:27) = (/ 9, 2,1,1, 3,1,2, 4,1,3, 10,5,1, 11,5,1, 15,5,2, 16,5,2, 20,5,3, 21,5,3 /)
	ind_quad_form(18,0:14) = (/ 7, 1,17, 2,16, 3,15, 5,13, 6,12, 7,11, 8,10 /)
	ind_lin_loss(18,0:16) = (/ 8, 1,17, 2,16, 3,15, 5,13, 6,12, 7,11, 8,10, 26,26 /)
	ind_lin_form(18,0:6) = (/ 3, 1,19, 5,23, 6,24 /)

	! Cluster 19: 4A3D
	ind_quad_loss(19,0:48) = (/ 24, 2,19, 3,19, 4,19, 5,24, 6,30, 7,30, 8,30, 9,30, 10,24&
			&, 11,30, 12,30, 13,30, 14,30, 15,24, 16,30, 17,30, 18,30, 19,30, 19,30&
			&, 20,24, 21,30, 22,30, 23,30, 24,30 /)
	ind_quad_loss_extra(19,0:18) = (/ 6, 2,1,2, 3,1,3, 4,1,4, 10,5,1, 15,5,2, 20,5,3 /)
	ind_quad_form(19,0:56) = (/ 28, 1,18, 2,17, 2,18, 2,19, 3,16, 3,17, 3,18, 3,19, 4,15&
			&, 4,16, 4,17, 4,18, 4,19, 5,14, 6,13, 6,14, 7,12, 7,13, 7,14&
			&, 8,11, 8,12, 8,13, 8,14, 9,10, 9,11, 9,12, 9,13, 9,14 /)
	ind_lin_loss(19,0:20) = (/ 10, 1,18, 2,17, 3,16, 4,15, 5,14, 6,13, 7,12, 8,11, 9,10&
			&, 26,26 /)
	ind_lin_form(19,0:2) = (/ 1, 5,24 /)

	! Cluster 20: 4D
	ind_quad_loss(20,0:48) = (/ 24, 1,21, 2,22, 3,23, 4,24, 6,21, 7,22, 8,23, 9,24, 10,20&
			&, 11,21, 12,22, 13,23, 14,24, 15,20, 16,21, 17,22, 18,23, 19,24, 20,20&
			&, 20,20, 21,21, 22,22, 23,23, 24,24 /)
	ind_quad_loss_extra(20,0:60) = (/ 20, 6,5,1, 7,5,1, 8,5,1, 9,5,1, 10,5,2, 11,5,2, 12,5,2, 13,5,2, 14,5,2&
			&, 15,5,3, 16,5,3, 17,5,3, 18,5,3, 19,5,3, 20,5,4, 20,5,4, 21,5,4, 22,5,4, 23,5,4&
			&, 24,5,4 /)
	ind_quad_form(20,0:14) = (/ 7, 5,15, 10,10, 10,15, 10,20, 15,15, 15,20, 20,20 /)
	ind_lin_loss(20,0:6) = (/ 3, 5,15, 10,10, 26,26 /)
	ind_lin_form(20,0:8) = (/ 4, 1,21, 2,22, 3,23, 4,24 /)

	! Cluster 21: 1A4D
	ind_quad_loss(21,0:48) = (/ 24, 1,22, 2,23, 3,24, 4,30, 6,22, 7,23, 8,24, 9,30, 10,21&
			&, 11,22, 12,23, 13,24, 14,30, 15,21, 16,22, 17,23, 18,24, 19,30, 20,21&
			&, 21,22, 21,22, 22,23, 23,24, 24,30 /)
	ind_quad_loss_extra(21,0:48) = (/ 16, 6,5,1, 7,5,1, 8,5,1, 10,5,2, 11,5,2, 12,5,2, 13,5,2, 15,5,3, 16,5,3&
			&, 17,5,3, 18,5,3, 20,5,4, 21,5,4, 21,5,4, 22,5,4, 23,5,4 /)
	ind_quad_form(21,0:26) = (/ 13, 1,20, 5,16, 6,15, 6,20, 10,11, 10,16, 10,21, 11,15, 11,20&
			&, 15,16, 15,21, 16,20, 20,21 /)
	ind_lin_loss(21,0:10) = (/ 5, 1,20, 5,16, 6,15, 10,11, 26,26 /)
	ind_lin_form(21,0:6) = (/ 3, 1,22, 2,23, 3,24 /)

	! Cluster 22: 2A4D
	ind_quad_loss(22,0:48) = (/ 24, 1,23, 2,24, 3,30, 4,30, 6,23, 7,24, 8,30, 9,30, 10,22&
			&, 11,23, 12,24, 13,30, 14,30, 15,22, 16,23, 17,24, 18,30, 19,30, 20,22&
			&, 21,23, 22,24, 22,24, 23,30, 24,30 /)
	ind_quad_loss_extra(22,0:36) = (/ 12, 6,5,1, 7,5,1, 10,5,2, 11,5,2, 12,5,2, 15,5,3, 16,5,3, 17,5,3, 20,5,4&
			&, 21,5,4, 22,5,4, 22,5,4 /)
	ind_quad_form(22,0:44) = (/ 22, 1,21, 2,20, 5,17, 6,16, 6,21, 7,15, 7,20, 10,12, 10,17&
			&, 10,22, 11,11, 11,16, 11,21, 12,15, 12,20, 15,17, 15,22, 16,16, 16,21&
			&, 17,20, 20,22, 21,21 /)
	ind_lin_loss(22,0:16) = (/ 8, 1,21, 2,20, 5,17, 6,16, 7,15, 10,12, 11,11, 26,26 /)
	ind_lin_form(22,0:4) = (/ 2, 1,23, 2,24 /)

	! Cluster 23: 3A4D
	ind_quad_loss(23,0:48) = (/ 24, 1,24, 2,30, 3,30, 4,30, 6,24, 7,30, 8,30, 9,30, 10,23&
			&, 11,24, 12,30, 13,30, 14,30, 15,23, 16,24, 17,30, 18,30, 19,30, 20,23&
			&, 21,24, 22,30, 23,30, 23,30, 24,30 /)
	ind_quad_loss_extra(23,0:21) = (/ 7, 6,5,1, 10,5,2, 11,5,2, 15,5,3, 16,5,3, 20,5,4, 21,5,4 /)
	ind_quad_form(23,0:56) = (/ 28, 1,22, 2,21, 3,20, 5,18, 6,17, 6,22, 7,16, 7,21, 8,15&
			&, 8,20, 10,13, 10,18, 10,23, 11,12, 11,17, 11,22, 12,16, 12,21, 13,15&
			&, 13,20, 15,18, 15,23, 16,17, 16,22, 17,21, 18,20, 20,23, 21,22 /)
	ind_lin_loss(23,0:20) = (/ 10, 1,22, 2,21, 3,20, 5,18, 6,17, 7,16, 8,15, 10,13, 11,12&
			&, 26,26 /)
	ind_lin_form(23,0:2) = (/ 1, 1,24 /)

	! Cluster 24: 4A4D
	ind_quad_loss(24,0:48) = (/ 24, 1,30, 2,30, 3,30, 4,30, 6,30, 7,30, 8,30, 9,30, 10,24&
			&, 11,30, 12,30, 13,30, 14,30, 15,24, 16,30, 17,30, 18,30, 19,30, 20,24&
			&, 21,30, 22,30, 23,30, 24,30, 24,30 /)
	ind_quad_loss_extra(24,0:9) = (/ 3, 10,5,2, 15,5,3, 20,5,4 /)
	ind_quad_form(24,0:74) = (/ 37, 1,23, 2,22, 3,21, 4,20, 5,19, 6,18, 6,23, 7,17, 7,22&
			&, 8,16, 8,21, 9,15, 9,20, 10,14, 10,19, 10,24, 11,13, 11,18, 11,23&
			&, 12,12, 12,17, 12,22, 13,16, 13,21, 14,15, 14,20, 15,19, 15,24, 16,18&
			&, 16,23, 17,17, 17,22, 18,21, 19,20, 20,24, 21,23, 22,22 /)
	ind_lin_loss(24,0:26) = (/ 13, 1,23, 2,22, 3,21, 4,20, 5,19, 6,18, 7,17, 8,16, 9,15&
			&, 10,14, 11,13, 12,12, 26,26 /)

	! Cluster 26: coag
	ind_lin_form(26,0:48) = (/ 24, 26,1, 26,2, 26,3, 26,4, 26,5, 26,6, 26,7, 26,8, 26,9&
			&, 26,10, 26,11, 26,12, 26,13, 26,14, 26,15, 26,16, 26,17, 26,18, 26,19&
			&, 26,20, 26,21, 26,22, 26,23, 26,24 /)

	! Cluster 30: out_neu
	ind_quad_form(30,0:156) = (/ 78, 1,24, 2,23, 2,24, 3,22, 3,23, 3,24, 4,21, 4,22, 4,23&
			&, 4,24, 6,19, 6,24, 7,18, 7,19, 7,23, 7,24, 8,17, 8,18, 8,19&
			&, 8,22, 8,23, 8,24, 9,16, 9,17, 9,18, 9,19, 9,21, 9,22, 9,23&
			&, 9,24, 11,14, 11,19, 11,24, 12,13, 12,14, 12,18, 12,19, 12,23, 12,24&
			&, 13,13, 13,14, 13,17, 13,18, 13,19, 13,22, 13,23, 13,24, 14,14, 14,16&
			&, 14,17, 14,18, 14,19, 14,21, 14,22, 14,23, 14,24, 16,19, 16,24, 17,18&
			&, 17,19, 17,23, 17,24, 18,18, 18,19, 18,22, 18,23, 18,24, 19,19, 19,21&
			&, 19,22, 19,23, 19,24, 21,24, 22,23, 22,24, 23,23, 23,24, 24,24 /)

	call get_rate_coefsD(coef_quad,coef_lin,coef)

end subroutine initialize_parametersD

!-----------------------------------------------------------

subroutine get_rate_coefsD(coef_quad,coef_lin,coef)
	implicit none
	real(kind(1.d0)) :: coef_quad(24,24,31),coef_lin(31,31,24)
	real(kind(1.d0)) :: K(24,24),E(24,24),cs(24)
	real(kind(1.d0)) :: coef(2)

	coef_quad = 0.d0
	coef_lin = 0.d0
	call get_collD(K,coef(1))
	call get_evapD(E,K,coef(1))
	call get_lossesD(cs)
	cs = coef(2)*cs


	coef_quad(1,1,2) = 0.5d0*K(1,1)	! 1A + 1A -> 2A
	coef_lin(1,1,2) = E(1,1)	! 2A -> 1A + 1A

	coef_quad(2,1,3) = K(1,2)	! 2A + 1A -> 3A
	coef_quad(1,2,3) = K(1,2)
	coef_lin(2,1,3) = E(1,2)	! 3A -> 2A + 1A
	coef_lin(1,2,3) = E(1,2)
	coef_quad(2,2,4) = 0.5d0*K(2,2)	! 2A + 2A -> 4A
	coef_lin(2,2,4) = E(2,2)	! 4A -> 2A + 2A

	coef_quad(3,1,4) = K(1,3)	! 3A + 1A -> 4A
	coef_quad(1,3,4) = K(1,3)
	coef_lin(3,1,4) = E(1,3)	! 4A -> 3A + 1A
	coef_lin(1,3,4) = E(1,3)
	coef_quad(3,2,4) = K(2,3)	! 3A + 2A -> boundary -> 4A
	coef_quad(2,3,4) = K(2,3)
	coef_quad(3,2,1) = K(2,3)	! 3A + 2A -> boundary -> 1A
	coef_quad(2,3,1) = K(2,3)
	coef_quad(3,3,4) = 0.5d0*K(3,3)	! 3A + 3A -> boundary -> 4A
	coef_quad(3,3,1) = 0.5d0*K(3,3)	! 3A + 3A -> boundary -> 1A

	coef_quad(4,2,4) = K(2,4)	! 4A + 2A -> boundary -> 4A
	coef_quad(2,4,4) = K(2,4)
	coef_quad(4,2,1) = K(2,4)	! 4A + 2A -> boundary -> 1A
	coef_quad(2,4,1) = K(2,4)
	coef_quad(4,3,4) = K(3,4)	! 4A + 3A -> boundary -> 4A
	coef_quad(3,4,4) = K(3,4)
	coef_quad(4,3,1) = K(3,4)	! 4A + 3A -> boundary -> 1A
	coef_quad(3,4,1) = K(3,4)
	coef_quad(4,4,4) = 0.5d0*K(4,4)	! 4A + 4A -> boundary -> 4A
	coef_quad(4,4,1) = 0.5d0*K(4,4)	! 4A + 4A -> boundary -> 1A

	coef_quad(5,1,6) = K(1,5)	! 1D + 1A -> 1A1D
	coef_quad(1,5,6) = K(1,5)
	coef_lin(5,1,6) = E(1,5)	! 1A1D -> 1D + 1A
	coef_lin(1,5,6) = E(1,5)
	coef_quad(5,2,7) = K(2,5)	! 1D + 2A -> 2A1D
	coef_quad(2,5,7) = K(2,5)
	coef_lin(5,2,7) = E(2,5)	! 2A1D -> 1D + 2A
	coef_lin(2,5,7) = E(2,5)
	coef_quad(5,3,8) = K(3,5)	! 1D + 3A -> 3A1D
	coef_quad(3,5,8) = K(3,5)
	coef_lin(5,3,8) = E(3,5)	! 3A1D -> 1D + 3A
	coef_lin(3,5,8) = E(3,5)
	coef_quad(5,4,9) = K(4,5)	! 1D + 4A -> 4A1D
	coef_quad(4,5,9) = K(4,5)
	coef_lin(5,4,9) = E(4,5)	! 4A1D -> 1D + 4A
	coef_lin(4,5,9) = E(4,5)
	coef_quad(5,5,10) = 0.5d0*K(5,5)	! 1D + 1D -> 2D
	coef_lin(5,5,10) = E(5,5)	! 2D -> 1D + 1D

	coef_quad(6,1,7) = K(1,6)	! 1A1D + 1A -> 2A1D
	coef_quad(1,6,7) = K(1,6)
	coef_lin(6,1,7) = E(1,6)	! 2A1D -> 1A1D + 1A
	coef_lin(1,6,7) = E(1,6)
	coef_quad(6,2,8) = K(2,6)	! 1A1D + 2A -> 3A1D
	coef_quad(2,6,8) = K(2,6)
	coef_lin(6,2,8) = E(2,6)	! 3A1D -> 1A1D + 2A
	coef_lin(2,6,8) = E(2,6)
	coef_quad(6,3,9) = K(3,6)	! 1A1D + 3A -> 4A1D
	coef_quad(3,6,9) = K(3,6)
	coef_lin(6,3,9) = E(3,6)	! 4A1D -> 1A1D + 3A
	coef_lin(3,6,9) = E(3,6)
	coef_quad(6,4,9) = K(4,6)	! 1A1D + 4A -> boundary -> 4A1D
	coef_quad(4,6,9) = K(4,6)
	coef_quad(6,4,1) = K(4,6)	! 1A1D + 4A -> boundary -> 1A
	coef_quad(4,6,1) = K(4,6)
	coef_quad(6,5,11) = K(5,6)	! 1A1D + 1D -> 1A2D
	coef_quad(5,6,11) = K(5,6)
	coef_lin(6,5,11) = E(5,6)	! 1A2D -> 1A1D + 1D
	coef_lin(5,6,11) = E(5,6)
	coef_quad(6,6,12) = 0.5d0*K(6,6)	! 1A1D + 1A1D -> 2A2D
	coef_lin(6,6,12) = E(6,6)	! 2A2D -> 1A1D + 1A1D

	coef_quad(7,1,8) = K(1,7)	! 2A1D + 1A -> 3A1D
	coef_quad(1,7,8) = K(1,7)
	coef_lin(7,1,8) = E(1,7)	! 3A1D -> 2A1D + 1A
	coef_lin(1,7,8) = E(1,7)
	coef_quad(7,2,9) = K(2,7)	! 2A1D + 2A -> 4A1D
	coef_quad(2,7,9) = K(2,7)
	coef_lin(7,2,9) = E(2,7)	! 4A1D -> 2A1D + 2A
	coef_lin(2,7,9) = E(2,7)
	coef_quad(7,3,9) = K(3,7)	! 2A1D + 3A -> boundary -> 4A1D
	coef_quad(3,7,9) = K(3,7)
	coef_quad(7,3,1) = K(3,7)	! 2A1D + 3A -> boundary -> 1A
	coef_quad(3,7,1) = K(3,7)
	coef_quad(7,4,9) = K(4,7)	! 2A1D + 4A -> boundary -> 4A1D
	coef_quad(4,7,9) = K(4,7)
	coef_quad(7,4,1) = K(4,7)	! 2A1D + 4A -> boundary -> 1A
	coef_quad(4,7,1) = K(4,7)
	coef_quad(7,5,12) = K(5,7)	! 2A1D + 1D -> 2A2D
	coef_quad(5,7,12) = K(5,7)
	coef_lin(7,5,12) = E(5,7)	! 2A2D -> 2A1D + 1D
	coef_lin(5,7,12) = E(5,7)
	coef_quad(7,6,13) = K(6,7)	! 2A1D + 1A1D -> 3A2D
	coef_quad(6,7,13) = K(6,7)
	coef_lin(7,6,13) = E(6,7)	! 3A2D -> 2A1D + 1A1D
	coef_lin(6,7,13) = E(6,7)
	coef_quad(7,7,14) = 0.5d0*K(7,7)	! 2A1D + 2A1D -> 4A2D
	coef_lin(7,7,14) = E(7,7)	! 4A2D -> 2A1D + 2A1D

	coef_quad(8,1,9) = K(1,8)	! 3A1D + 1A -> 4A1D
	coef_quad(1,8,9) = K(1,8)
	coef_lin(8,1,9) = E(1,8)	! 4A1D -> 3A1D + 1A
	coef_lin(1,8,9) = E(1,8)
	coef_quad(8,2,9) = K(2,8)	! 3A1D + 2A -> boundary -> 4A1D
	coef_quad(2,8,9) = K(2,8)
	coef_quad(8,2,1) = K(2,8)	! 3A1D + 2A -> boundary -> 1A
	coef_quad(2,8,1) = K(2,8)
	coef_quad(8,3,9) = K(3,8)	! 3A1D + 3A -> boundary -> 4A1D
	coef_quad(3,8,9) = K(3,8)
	coef_quad(8,3,1) = K(3,8)	! 3A1D + 3A -> boundary -> 1A
	coef_quad(3,8,1) = K(3,8)
	coef_quad(8,4,9) = K(4,8)	! 3A1D + 4A -> boundary -> 4A1D
	coef_quad(4,8,9) = K(4,8)
	coef_quad(8,4,1) = K(4,8)	! 3A1D + 4A -> boundary -> 1A
	coef_quad(4,8,1) = K(4,8)
	coef_quad(8,5,13) = K(5,8)	! 3A1D + 1D -> 3A2D
	coef_quad(5,8,13) = K(5,8)
	coef_lin(8,5,13) = E(5,8)	! 3A2D -> 3A1D + 1D
	coef_lin(5,8,13) = E(5,8)
	coef_quad(8,6,14) = K(6,8)	! 3A1D + 1A1D -> 4A2D
	coef_quad(6,8,14) = K(6,8)
	coef_lin(8,6,14) = E(6,8)	! 4A2D -> 3A1D + 1A1D
	coef_lin(6,8,14) = E(6,8)
	coef_quad(8,7,14) = K(7,8)	! 3A1D + 2A1D -> boundary -> 4A2D
	coef_quad(7,8,14) = K(7,8)
	coef_quad(8,7,1) = K(7,8)	! 3A1D + 2A1D -> boundary -> 1A
	coef_quad(7,8,1) = K(7,8)
	coef_quad(8,8,14) = 0.5d0*K(8,8)	! 3A1D + 3A1D -> boundary -> 4A2D
	coef_quad(8,8,1) = 0.5d0*K(8,8)	! 3A1D + 3A1D -> boundary -> 1A

	coef_quad(9,2,9) = K(2,9)	! 4A1D + 2A -> boundary -> 4A1D
	coef_quad(2,9,9) = K(2,9)
	coef_quad(9,2,1) = K(2,9)	! 4A1D + 2A -> boundary -> 1A
	coef_quad(2,9,1) = K(2,9)
	coef_quad(9,3,9) = K(3,9)	! 4A1D + 3A -> boundary -> 4A1D
	coef_quad(3,9,9) = K(3,9)
	coef_quad(9,3,1) = K(3,9)	! 4A1D + 3A -> boundary -> 1A
	coef_quad(3,9,1) = K(3,9)
	coef_quad(9,4,9) = K(4,9)	! 4A1D + 4A -> boundary -> 4A1D
	coef_quad(4,9,9) = K(4,9)
	coef_quad(9,4,1) = K(4,9)	! 4A1D + 4A -> boundary -> 1A
	coef_quad(4,9,1) = K(4,9)
	coef_quad(9,5,14) = K(5,9)	! 4A1D + 1D -> 4A2D
	coef_quad(5,9,14) = K(5,9)
	coef_lin(9,5,14) = E(5,9)	! 4A2D -> 4A1D + 1D
	coef_lin(5,9,14) = E(5,9)
	coef_quad(9,6,14) = K(6,9)	! 4A1D + 1A1D -> boundary -> 4A2D
	coef_quad(6,9,14) = K(6,9)
	coef_quad(9,6,1) = K(6,9)	! 4A1D + 1A1D -> boundary -> 1A
	coef_quad(6,9,1) = K(6,9)
	coef_quad(9,7,14) = K(7,9)	! 4A1D + 2A1D -> boundary -> 4A2D
	coef_quad(7,9,14) = K(7,9)
	coef_quad(9,7,1) = K(7,9)	! 4A1D + 2A1D -> boundary -> 1A
	coef_quad(7,9,1) = K(7,9)
	coef_quad(9,8,14) = K(8,9)	! 4A1D + 3A1D -> boundary -> 4A2D
	coef_quad(8,9,14) = K(8,9)
	coef_quad(9,8,1) = K(8,9)	! 4A1D + 3A1D -> boundary -> 1A
	coef_quad(8,9,1) = K(8,9)
	coef_quad(9,9,14) = 0.5d0*K(9,9)	! 4A1D + 4A1D -> boundary -> 4A2D
	coef_quad(9,9,1) = 0.5d0*K(9,9)	! 4A1D + 4A1D -> boundary -> 1A

	coef_quad(10,1,11) = K(1,10)	! 2D + 1A -> 1A2D
	coef_quad(1,10,11) = K(1,10)
	coef_lin(10,1,11) = E(1,10)	! 1A2D -> 2D + 1A
	coef_lin(1,10,11) = E(1,10)
	coef_quad(10,2,12) = K(2,10)	! 2D + 2A -> 2A2D
	coef_quad(2,10,12) = K(2,10)
	coef_lin(10,2,12) = E(2,10)	! 2A2D -> 2D + 2A
	coef_lin(2,10,12) = E(2,10)
	coef_quad(10,3,13) = K(3,10)	! 2D + 3A -> 3A2D
	coef_quad(3,10,13) = K(3,10)
	coef_lin(10,3,13) = E(3,10)	! 3A2D -> 2D + 3A
	coef_lin(3,10,13) = E(3,10)
	coef_quad(10,4,14) = K(4,10)	! 2D + 4A -> 4A2D
	coef_quad(4,10,14) = K(4,10)
	coef_lin(10,4,14) = E(4,10)	! 4A2D -> 2D + 4A
	coef_lin(4,10,14) = E(4,10)
	coef_quad(10,5,15) = K(5,10)	! 2D + 1D -> 3D
	coef_quad(5,10,15) = K(5,10)
	coef_lin(10,5,15) = E(5,10)	! 3D -> 2D + 1D
	coef_lin(5,10,15) = E(5,10)
	coef_quad(10,6,16) = K(6,10)	! 2D + 1A1D -> 1A3D
	coef_quad(6,10,16) = K(6,10)
	coef_lin(10,6,16) = E(6,10)	! 1A3D -> 2D + 1A1D
	coef_lin(6,10,16) = E(6,10)
	coef_quad(10,7,17) = K(7,10)	! 2D + 2A1D -> 2A3D
	coef_quad(7,10,17) = K(7,10)
	coef_lin(10,7,17) = E(7,10)	! 2A3D -> 2D + 2A1D
	coef_lin(7,10,17) = E(7,10)
	coef_quad(10,8,18) = K(8,10)	! 2D + 3A1D -> 3A3D
	coef_quad(8,10,18) = K(8,10)
	coef_lin(10,8,18) = E(8,10)	! 3A3D -> 2D + 3A1D
	coef_lin(8,10,18) = E(8,10)
	coef_quad(10,9,19) = K(9,10)	! 2D + 4A1D -> 4A3D
	coef_quad(9,10,19) = K(9,10)
	coef_lin(10,9,19) = E(9,10)	! 4A3D -> 2D + 4A1D
	coef_lin(9,10,19) = E(9,10)
	coef_quad(10,10,20) = 0.5d0*K(10,10)	! 2D + 2D -> 4D
	coef_lin(10,10,20) = E(10,10)	! 4D -> 2D + 2D

	coef_quad(11,1,12) = K(1,11)	! 1A2D + 1A -> 2A2D
	coef_quad(1,11,12) = K(1,11)
	coef_lin(11,1,12) = E(1,11)	! 2A2D -> 1A2D + 1A
	coef_lin(1,11,12) = E(1,11)
	coef_quad(11,2,13) = K(2,11)	! 1A2D + 2A -> 3A2D
	coef_quad(2,11,13) = K(2,11)
	coef_lin(11,2,13) = E(2,11)	! 3A2D -> 1A2D + 2A
	coef_lin(2,11,13) = E(2,11)
	coef_quad(11,3,14) = K(3,11)	! 1A2D + 3A -> 4A2D
	coef_quad(3,11,14) = K(3,11)
	coef_lin(11,3,14) = E(3,11)	! 4A2D -> 1A2D + 3A
	coef_lin(3,11,14) = E(3,11)
	coef_quad(11,4,14) = K(4,11)	! 1A2D + 4A -> boundary -> 4A2D
	coef_quad(4,11,14) = K(4,11)
	coef_quad(11,4,1) = K(4,11)	! 1A2D + 4A -> boundary -> 1A
	coef_quad(4,11,1) = K(4,11)
	coef_quad(11,5,16) = K(5,11)	! 1A2D + 1D -> 1A3D
	coef_quad(5,11,16) = K(5,11)
	coef_lin(11,5,16) = E(5,11)	! 1A3D -> 1A2D + 1D
	coef_lin(5,11,16) = E(5,11)
	coef_quad(11,6,17) = K(6,11)	! 1A2D + 1A1D -> 2A3D
	coef_quad(6,11,17) = K(6,11)
	coef_lin(11,6,17) = E(6,11)	! 2A3D -> 1A2D + 1A1D
	coef_lin(6,11,17) = E(6,11)
	coef_quad(11,7,18) = K(7,11)	! 1A2D + 2A1D -> 3A3D
	coef_quad(7,11,18) = K(7,11)
	coef_lin(11,7,18) = E(7,11)	! 3A3D -> 1A2D + 2A1D
	coef_lin(7,11,18) = E(7,11)
	coef_quad(11,8,19) = K(8,11)	! 1A2D + 3A1D -> 4A3D
	coef_quad(8,11,19) = K(8,11)
	coef_lin(11,8,19) = E(8,11)	! 4A3D -> 1A2D + 3A1D
	coef_lin(8,11,19) = E(8,11)
	coef_quad(11,9,19) = K(9,11)	! 1A2D + 4A1D -> boundary -> 4A3D
	coef_quad(9,11,19) = K(9,11)
	coef_quad(11,9,1) = K(9,11)	! 1A2D + 4A1D -> boundary -> 1A
	coef_quad(9,11,1) = K(9,11)
	coef_quad(11,10,21) = K(10,11)	! 1A2D + 2D -> 1A4D
	coef_quad(10,11,21) = K(10,11)
	coef_lin(11,10,21) = E(10,11)	! 1A4D -> 1A2D + 2D
	coef_lin(10,11,21) = E(10,11)
	coef_quad(11,11,22) = 0.5d0*K(11,11)	! 1A2D + 1A2D -> 2A4D
	coef_lin(11,11,22) = E(11,11)	! 2A4D -> 1A2D + 1A2D

	coef_quad(12,1,13) = K(1,12)	! 2A2D + 1A -> 3A2D
	coef_quad(1,12,13) = K(1,12)
	coef_lin(12,1,13) = E(1,12)	! 3A2D -> 2A2D + 1A
	coef_lin(1,12,13) = E(1,12)
	coef_quad(12,2,14) = K(2,12)	! 2A2D + 2A -> 4A2D
	coef_quad(2,12,14) = K(2,12)
	coef_lin(12,2,14) = E(2,12)	! 4A2D -> 2A2D + 2A
	coef_lin(2,12,14) = E(2,12)
	coef_quad(12,3,14) = K(3,12)	! 2A2D + 3A -> boundary -> 4A2D
	coef_quad(3,12,14) = K(3,12)
	coef_quad(12,3,1) = K(3,12)	! 2A2D + 3A -> boundary -> 1A
	coef_quad(3,12,1) = K(3,12)
	coef_quad(12,4,14) = K(4,12)	! 2A2D + 4A -> boundary -> 4A2D
	coef_quad(4,12,14) = K(4,12)
	coef_quad(12,4,1) = K(4,12)	! 2A2D + 4A -> boundary -> 1A
	coef_quad(4,12,1) = K(4,12)
	coef_quad(12,5,17) = K(5,12)	! 2A2D + 1D -> 2A3D
	coef_quad(5,12,17) = K(5,12)
	coef_lin(12,5,17) = E(5,12)	! 2A3D -> 2A2D + 1D
	coef_lin(5,12,17) = E(5,12)
	coef_quad(12,6,18) = K(6,12)	! 2A2D + 1A1D -> 3A3D
	coef_quad(6,12,18) = K(6,12)
	coef_lin(12,6,18) = E(6,12)	! 3A3D -> 2A2D + 1A1D
	coef_lin(6,12,18) = E(6,12)
	coef_quad(12,7,19) = K(7,12)	! 2A2D + 2A1D -> 4A3D
	coef_quad(7,12,19) = K(7,12)
	coef_lin(12,7,19) = E(7,12)	! 4A3D -> 2A2D + 2A1D
	coef_lin(7,12,19) = E(7,12)
	coef_quad(12,8,19) = K(8,12)	! 2A2D + 3A1D -> boundary -> 4A3D
	coef_quad(8,12,19) = K(8,12)
	coef_quad(12,8,1) = K(8,12)	! 2A2D + 3A1D -> boundary -> 1A
	coef_quad(8,12,1) = K(8,12)
	coef_quad(12,9,19) = K(9,12)	! 2A2D + 4A1D -> boundary -> 4A3D
	coef_quad(9,12,19) = K(9,12)
	coef_quad(12,9,1) = K(9,12)	! 2A2D + 4A1D -> boundary -> 1A
	coef_quad(9,12,1) = K(9,12)
	coef_quad(12,10,22) = K(10,12)	! 2A2D + 2D -> 2A4D
	coef_quad(10,12,22) = K(10,12)
	coef_lin(12,10,22) = E(10,12)	! 2A4D -> 2A2D + 2D
	coef_lin(10,12,22) = E(10,12)
	coef_quad(12,11,23) = K(11,12)	! 2A2D + 1A2D -> 3A4D
	coef_quad(11,12,23) = K(11,12)
	coef_lin(12,11,23) = E(11,12)	! 3A4D -> 2A2D + 1A2D
	coef_lin(11,12,23) = E(11,12)
	coef_quad(12,12,24) = 0.5d0*K(12,12)	! 2A2D + 2A2D -> 4A4D
	coef_lin(12,12,24) = E(12,12)	! 4A4D -> 2A2D + 2A2D

	coef_quad(13,1,14) = K(1,13)	! 3A2D + 1A -> 4A2D
	coef_quad(1,13,14) = K(1,13)
	coef_lin(13,1,14) = E(1,13)	! 4A2D -> 3A2D + 1A
	coef_lin(1,13,14) = E(1,13)
	coef_quad(13,2,14) = K(2,13)	! 3A2D + 2A -> boundary -> 4A2D
	coef_quad(2,13,14) = K(2,13)
	coef_quad(13,2,1) = K(2,13)	! 3A2D + 2A -> boundary -> 1A
	coef_quad(2,13,1) = K(2,13)
	coef_quad(13,3,14) = K(3,13)	! 3A2D + 3A -> boundary -> 4A2D
	coef_quad(3,13,14) = K(3,13)
	coef_quad(13,3,1) = K(3,13)	! 3A2D + 3A -> boundary -> 1A
	coef_quad(3,13,1) = K(3,13)
	coef_quad(13,4,14) = K(4,13)	! 3A2D + 4A -> boundary -> 4A2D
	coef_quad(4,13,14) = K(4,13)
	coef_quad(13,4,1) = K(4,13)	! 3A2D + 4A -> boundary -> 1A
	coef_quad(4,13,1) = K(4,13)
	coef_quad(13,5,18) = K(5,13)	! 3A2D + 1D -> 3A3D
	coef_quad(5,13,18) = K(5,13)
	coef_lin(13,5,18) = E(5,13)	! 3A3D -> 3A2D + 1D
	coef_lin(5,13,18) = E(5,13)
	coef_quad(13,6,19) = K(6,13)	! 3A2D + 1A1D -> 4A3D
	coef_quad(6,13,19) = K(6,13)
	coef_lin(13,6,19) = E(6,13)	! 4A3D -> 3A2D + 1A1D
	coef_lin(6,13,19) = E(6,13)
	coef_quad(13,7,19) = K(7,13)	! 3A2D + 2A1D -> boundary -> 4A3D
	coef_quad(7,13,19) = K(7,13)
	coef_quad(13,7,1) = K(7,13)	! 3A2D + 2A1D -> boundary -> 1A
	coef_quad(7,13,1) = K(7,13)
	coef_quad(13,8,19) = K(8,13)	! 3A2D + 3A1D -> boundary -> 4A3D
	coef_quad(8,13,19) = K(8,13)
	coef_quad(13,8,1) = K(8,13)	! 3A2D + 3A1D -> boundary -> 1A
	coef_quad(8,13,1) = K(8,13)
	coef_quad(13,9,19) = K(9,13)	! 3A2D + 4A1D -> boundary -> 4A3D
	coef_quad(9,13,19) = K(9,13)
	coef_quad(13,9,1) = K(9,13)	! 3A2D + 4A1D -> boundary -> 1A
	coef_quad(9,13,1) = K(9,13)
	coef_quad(13,10,23) = K(10,13)	! 3A2D + 2D -> 3A4D
	coef_quad(10,13,23) = K(10,13)
	coef_lin(13,10,23) = E(10,13)	! 3A4D -> 3A2D + 2D
	coef_lin(10,13,23) = E(10,13)
	coef_quad(13,11,24) = K(11,13)	! 3A2D + 1A2D -> 4A4D
	coef_quad(11,13,24) = K(11,13)
	coef_lin(13,11,24) = E(11,13)	! 4A4D -> 3A2D + 1A2D
	coef_lin(11,13,24) = E(11,13)
	coef_quad(13,12,30) = K(12,13)	! 3A2D + 2A2D -> out_neu
	coef_quad(12,13,30) = K(12,13)
	coef_quad(13,13,30) = 0.5d0*K(13,13)	! 3A2D + 3A2D -> out_neu

	coef_quad(14,2,14) = K(2,14)	! 4A2D + 2A -> boundary -> 4A2D
	coef_quad(2,14,14) = K(2,14)
	coef_quad(14,2,1) = K(2,14)	! 4A2D + 2A -> boundary -> 1A
	coef_quad(2,14,1) = K(2,14)
	coef_quad(14,3,14) = K(3,14)	! 4A2D + 3A -> boundary -> 4A2D
	coef_quad(3,14,14) = K(3,14)
	coef_quad(14,3,1) = K(3,14)	! 4A2D + 3A -> boundary -> 1A
	coef_quad(3,14,1) = K(3,14)
	coef_quad(14,4,14) = K(4,14)	! 4A2D + 4A -> boundary -> 4A2D
	coef_quad(4,14,14) = K(4,14)
	coef_quad(14,4,1) = K(4,14)	! 4A2D + 4A -> boundary -> 1A
	coef_quad(4,14,1) = K(4,14)
	coef_quad(14,5,19) = K(5,14)	! 4A2D + 1D -> 4A3D
	coef_quad(5,14,19) = K(5,14)
	coef_lin(14,5,19) = E(5,14)	! 4A3D -> 4A2D + 1D
	coef_lin(5,14,19) = E(5,14)
	coef_quad(14,6,19) = K(6,14)	! 4A2D + 1A1D -> boundary -> 4A3D
	coef_quad(6,14,19) = K(6,14)
	coef_quad(14,6,1) = K(6,14)	! 4A2D + 1A1D -> boundary -> 1A
	coef_quad(6,14,1) = K(6,14)
	coef_quad(14,7,19) = K(7,14)	! 4A2D + 2A1D -> boundary -> 4A3D
	coef_quad(7,14,19) = K(7,14)
	coef_quad(14,7,1) = K(7,14)	! 4A2D + 2A1D -> boundary -> 1A
	coef_quad(7,14,1) = K(7,14)
	coef_quad(14,8,19) = K(8,14)	! 4A2D + 3A1D -> boundary -> 4A3D
	coef_quad(8,14,19) = K(8,14)
	coef_quad(14,8,1) = K(8,14)	! 4A2D + 3A1D -> boundary -> 1A
	coef_quad(8,14,1) = K(8,14)
	coef_quad(14,9,19) = K(9,14)	! 4A2D + 4A1D -> boundary -> 4A3D
	coef_quad(9,14,19) = K(9,14)
	coef_quad(14,9,1) = K(9,14)	! 4A2D + 4A1D -> boundary -> 1A
	coef_quad(9,14,1) = K(9,14)
	coef_quad(14,10,24) = K(10,14)	! 4A2D + 2D -> 4A4D
	coef_quad(10,14,24) = K(10,14)
	coef_lin(14,10,24) = E(10,14)	! 4A4D -> 4A2D + 2D
	coef_lin(10,14,24) = E(10,14)
	coef_quad(14,11,30) = K(11,14)	! 4A2D + 1A2D -> out_neu
	coef_quad(11,14,30) = K(11,14)
	coef_quad(14,12,30) = K(12,14)	! 4A2D + 2A2D -> out_neu
	coef_quad(12,14,30) = K(12,14)
	coef_quad(14,13,30) = K(13,14)	! 4A2D + 3A2D -> out_neu
	coef_quad(13,14,30) = K(13,14)
	coef_quad(14,14,30) = 0.5d0*K(14,14)	! 4A2D + 4A2D -> out_neu

	coef_quad(15,1,16) = K(1,15)	! 3D + 1A -> 1A3D
	coef_quad(1,15,16) = K(1,15)
	coef_lin(15,1,16) = E(1,15)	! 1A3D -> 3D + 1A
	coef_lin(1,15,16) = E(1,15)
	coef_quad(15,2,17) = K(2,15)	! 3D + 2A -> 2A3D
	coef_quad(2,15,17) = K(2,15)
	coef_lin(15,2,17) = E(2,15)	! 2A3D -> 3D + 2A
	coef_lin(2,15,17) = E(2,15)
	coef_quad(15,3,18) = K(3,15)	! 3D + 3A -> 3A3D
	coef_quad(3,15,18) = K(3,15)
	coef_lin(15,3,18) = E(3,15)	! 3A3D -> 3D + 3A
	coef_lin(3,15,18) = E(3,15)
	coef_quad(15,4,19) = K(4,15)	! 3D + 4A -> 4A3D
	coef_quad(4,15,19) = K(4,15)
	coef_lin(15,4,19) = E(4,15)	! 4A3D -> 3D + 4A
	coef_lin(4,15,19) = E(4,15)
	coef_quad(15,5,20) = K(5,15)	! 3D + 1D -> 4D
	coef_quad(5,15,20) = K(5,15)
	coef_lin(15,5,20) = E(5,15)	! 4D -> 3D + 1D
	coef_lin(5,15,20) = E(5,15)
	coef_quad(15,6,21) = K(6,15)	! 3D + 1A1D -> 1A4D
	coef_quad(6,15,21) = K(6,15)
	coef_lin(15,6,21) = E(6,15)	! 1A4D -> 3D + 1A1D
	coef_lin(6,15,21) = E(6,15)
	coef_quad(15,7,22) = K(7,15)	! 3D + 2A1D -> 2A4D
	coef_quad(7,15,22) = K(7,15)
	coef_lin(15,7,22) = E(7,15)	! 2A4D -> 3D + 2A1D
	coef_lin(7,15,22) = E(7,15)
	coef_quad(15,8,23) = K(8,15)	! 3D + 3A1D -> 3A4D
	coef_quad(8,15,23) = K(8,15)
	coef_lin(15,8,23) = E(8,15)	! 3A4D -> 3D + 3A1D
	coef_lin(8,15,23) = E(8,15)
	coef_quad(15,9,24) = K(9,15)	! 3D + 4A1D -> 4A4D
	coef_quad(9,15,24) = K(9,15)
	coef_lin(15,9,24) = E(9,15)	! 4A4D -> 3D + 4A1D
	coef_lin(9,15,24) = E(9,15)
	coef_quad(15,10,20) = K(10,15)	! 3D + 2D -> boundary -> 4D
	coef_quad(10,15,20) = K(10,15)
	coef_quad(15,10,5) = K(10,15)	! 3D + 2D -> boundary -> 1D
	coef_quad(10,15,5) = K(10,15)
	coef_quad(15,11,21) = K(11,15)	! 3D + 1A2D -> boundary -> 1A4D
	coef_quad(11,15,21) = K(11,15)
	coef_quad(15,11,5) = K(11,15)	! 3D + 1A2D -> boundary -> 1D
	coef_quad(11,15,5) = K(11,15)
	coef_quad(15,12,22) = K(12,15)	! 3D + 2A2D -> boundary -> 2A4D
	coef_quad(12,15,22) = K(12,15)
	coef_quad(15,12,5) = K(12,15)	! 3D + 2A2D -> boundary -> 1D
	coef_quad(12,15,5) = K(12,15)
	coef_quad(15,13,23) = K(13,15)	! 3D + 3A2D -> boundary -> 3A4D
	coef_quad(13,15,23) = K(13,15)
	coef_quad(15,13,5) = K(13,15)	! 3D + 3A2D -> boundary -> 1D
	coef_quad(13,15,5) = K(13,15)
	coef_quad(15,14,24) = K(14,15)	! 3D + 4A2D -> boundary -> 4A4D
	coef_quad(14,15,24) = K(14,15)
	coef_quad(15,14,5) = K(14,15)	! 3D + 4A2D -> boundary -> 1D
	coef_quad(14,15,5) = K(14,15)
	coef_quad(15,15,20) = 0.5d0*K(15,15)	! 3D + 3D -> boundary -> 4D
	coef_quad(15,15,5) = 0.5d0*K(15,15)	! 3D + 3D -> boundary -> 1D

	coef_quad(16,1,17) = K(1,16)	! 1A3D + 1A -> 2A3D
	coef_quad(1,16,17) = K(1,16)
	coef_lin(16,1,17) = E(1,16)	! 2A3D -> 1A3D + 1A
	coef_lin(1,16,17) = E(1,16)
	coef_quad(16,2,18) = K(2,16)	! 1A3D + 2A -> 3A3D
	coef_quad(2,16,18) = K(2,16)
	coef_lin(16,2,18) = E(2,16)	! 3A3D -> 1A3D + 2A
	coef_lin(2,16,18) = E(2,16)
	coef_quad(16,3,19) = K(3,16)	! 1A3D + 3A -> 4A3D
	coef_quad(3,16,19) = K(3,16)
	coef_lin(16,3,19) = E(3,16)	! 4A3D -> 1A3D + 3A
	coef_lin(3,16,19) = E(3,16)
	coef_quad(16,4,19) = K(4,16)	! 1A3D + 4A -> boundary -> 4A3D
	coef_quad(4,16,19) = K(4,16)
	coef_quad(16,4,1) = K(4,16)	! 1A3D + 4A -> boundary -> 1A
	coef_quad(4,16,1) = K(4,16)
	coef_quad(16,5,21) = K(5,16)	! 1A3D + 1D -> 1A4D
	coef_quad(5,16,21) = K(5,16)
	coef_lin(16,5,21) = E(5,16)	! 1A4D -> 1A3D + 1D
	coef_lin(5,16,21) = E(5,16)
	coef_quad(16,6,22) = K(6,16)	! 1A3D + 1A1D -> 2A4D
	coef_quad(6,16,22) = K(6,16)
	coef_lin(16,6,22) = E(6,16)	! 2A4D -> 1A3D + 1A1D
	coef_lin(6,16,22) = E(6,16)
	coef_quad(16,7,23) = K(7,16)	! 1A3D + 2A1D -> 3A4D
	coef_quad(7,16,23) = K(7,16)
	coef_lin(16,7,23) = E(7,16)	! 3A4D -> 1A3D + 2A1D
	coef_lin(7,16,23) = E(7,16)
	coef_quad(16,8,24) = K(8,16)	! 1A3D + 3A1D -> 4A4D
	coef_quad(8,16,24) = K(8,16)
	coef_lin(16,8,24) = E(8,16)	! 4A4D -> 1A3D + 3A1D
	coef_lin(8,16,24) = E(8,16)
	coef_quad(16,9,30) = K(9,16)	! 1A3D + 4A1D -> out_neu
	coef_quad(9,16,30) = K(9,16)
	coef_quad(16,10,21) = K(10,16)	! 1A3D + 2D -> boundary -> 1A4D
	coef_quad(10,16,21) = K(10,16)
	coef_quad(16,10,5) = K(10,16)	! 1A3D + 2D -> boundary -> 1D
	coef_quad(10,16,5) = K(10,16)
	coef_quad(16,11,22) = K(11,16)	! 1A3D + 1A2D -> boundary -> 2A4D
	coef_quad(11,16,22) = K(11,16)
	coef_quad(16,11,5) = K(11,16)	! 1A3D + 1A2D -> boundary -> 1D
	coef_quad(11,16,5) = K(11,16)
	coef_quad(16,12,23) = K(12,16)	! 1A3D + 2A2D -> boundary -> 3A4D
	coef_quad(12,16,23) = K(12,16)
	coef_quad(16,12,5) = K(12,16)	! 1A3D + 2A2D -> boundary -> 1D
	coef_quad(12,16,5) = K(12,16)
	coef_quad(16,13,24) = K(13,16)	! 1A3D + 3A2D -> boundary -> 4A4D
	coef_quad(13,16,24) = K(13,16)
	coef_quad(16,13,5) = K(13,16)	! 1A3D + 3A2D -> boundary -> 1D
	coef_quad(13,16,5) = K(13,16)
	coef_quad(16,14,30) = K(14,16)	! 1A3D + 4A2D -> out_neu
	coef_quad(14,16,30) = K(14,16)
	coef_quad(16,15,21) = K(15,16)	! 1A3D + 3D -> boundary -> 1A4D
	coef_quad(15,16,21) = K(15,16)
	coef_quad(16,15,5) = K(15,16)	! 1A3D + 3D -> boundary -> 1D
	coef_quad(15,16,5) = K(15,16)
	coef_quad(16,16,22) = 0.5d0*K(16,16)	! 1A3D + 1A3D -> boundary -> 2A4D
	coef_quad(16,16,5) = 0.5d0*K(16,16)	! 1A3D + 1A3D -> boundary -> 1D

	coef_quad(17,1,18) = K(1,17)	! 2A3D + 1A -> 3A3D
	coef_quad(1,17,18) = K(1,17)
	coef_lin(17,1,18) = E(1,17)	! 3A3D -> 2A3D + 1A
	coef_lin(1,17,18) = E(1,17)
	coef_quad(17,2,19) = K(2,17)	! 2A3D + 2A -> 4A3D
	coef_quad(2,17,19) = K(2,17)
	coef_lin(17,2,19) = E(2,17)	! 4A3D -> 2A3D + 2A
	coef_lin(2,17,19) = E(2,17)
	coef_quad(17,3,19) = K(3,17)	! 2A3D + 3A -> boundary -> 4A3D
	coef_quad(3,17,19) = K(3,17)
	coef_quad(17,3,1) = K(3,17)	! 2A3D + 3A -> boundary -> 1A
	coef_quad(3,17,1) = K(3,17)
	coef_quad(17,4,19) = K(4,17)	! 2A3D + 4A -> boundary -> 4A3D
	coef_quad(4,17,19) = K(4,17)
	coef_quad(17,4,1) = K(4,17)	! 2A3D + 4A -> boundary -> 1A
	coef_quad(4,17,1) = K(4,17)
	coef_quad(17,5,22) = K(5,17)	! 2A3D + 1D -> 2A4D
	coef_quad(5,17,22) = K(5,17)
	coef_lin(17,5,22) = E(5,17)	! 2A4D -> 2A3D + 1D
	coef_lin(5,17,22) = E(5,17)
	coef_quad(17,6,23) = K(6,17)	! 2A3D + 1A1D -> 3A4D
	coef_quad(6,17,23) = K(6,17)
	coef_lin(17,6,23) = E(6,17)	! 3A4D -> 2A3D + 1A1D
	coef_lin(6,17,23) = E(6,17)
	coef_quad(17,7,24) = K(7,17)	! 2A3D + 2A1D -> 4A4D
	coef_quad(7,17,24) = K(7,17)
	coef_lin(17,7,24) = E(7,17)	! 4A4D -> 2A3D + 2A1D
	coef_lin(7,17,24) = E(7,17)
	coef_quad(17,8,30) = K(8,17)	! 2A3D + 3A1D -> out_neu
	coef_quad(8,17,30) = K(8,17)
	coef_quad(17,9,30) = K(9,17)	! 2A3D + 4A1D -> out_neu
	coef_quad(9,17,30) = K(9,17)
	coef_quad(17,10,22) = K(10,17)	! 2A3D + 2D -> boundary -> 2A4D
	coef_quad(10,17,22) = K(10,17)
	coef_quad(17,10,5) = K(10,17)	! 2A3D + 2D -> boundary -> 1D
	coef_quad(10,17,5) = K(10,17)
	coef_quad(17,11,23) = K(11,17)	! 2A3D + 1A2D -> boundary -> 3A4D
	coef_quad(11,17,23) = K(11,17)
	coef_quad(17,11,5) = K(11,17)	! 2A3D + 1A2D -> boundary -> 1D
	coef_quad(11,17,5) = K(11,17)
	coef_quad(17,12,24) = K(12,17)	! 2A3D + 2A2D -> boundary -> 4A4D
	coef_quad(12,17,24) = K(12,17)
	coef_quad(17,12,5) = K(12,17)	! 2A3D + 2A2D -> boundary -> 1D
	coef_quad(12,17,5) = K(12,17)
	coef_quad(17,13,30) = K(13,17)	! 2A3D + 3A2D -> out_neu
	coef_quad(13,17,30) = K(13,17)
	coef_quad(17,14,30) = K(14,17)	! 2A3D + 4A2D -> out_neu
	coef_quad(14,17,30) = K(14,17)
	coef_quad(17,15,22) = K(15,17)	! 2A3D + 3D -> boundary -> 2A4D
	coef_quad(15,17,22) = K(15,17)
	coef_quad(17,15,5) = K(15,17)	! 2A3D + 3D -> boundary -> 1D
	coef_quad(15,17,5) = K(15,17)
	coef_quad(17,16,23) = K(16,17)	! 2A3D + 1A3D -> boundary -> 3A4D
	coef_quad(16,17,23) = K(16,17)
	coef_quad(17,16,5) = K(16,17)	! 2A3D + 1A3D -> boundary -> 1D
	coef_quad(16,17,5) = K(16,17)
	coef_quad(17,17,24) = 0.5d0*K(17,17)	! 2A3D + 2A3D -> boundary -> 4A4D
	coef_quad(17,17,5) = 0.5d0*K(17,17)	! 2A3D + 2A3D -> boundary -> 1D

	coef_quad(18,1,19) = K(1,18)	! 3A3D + 1A -> 4A3D
	coef_quad(1,18,19) = K(1,18)
	coef_lin(18,1,19) = E(1,18)	! 4A3D -> 3A3D + 1A
	coef_lin(1,18,19) = E(1,18)
	coef_quad(18,2,19) = K(2,18)	! 3A3D + 2A -> boundary -> 4A3D
	coef_quad(2,18,19) = K(2,18)
	coef_quad(18,2,1) = K(2,18)	! 3A3D + 2A -> boundary -> 1A
	coef_quad(2,18,1) = K(2,18)
	coef_quad(18,3,19) = K(3,18)	! 3A3D + 3A -> boundary -> 4A3D
	coef_quad(3,18,19) = K(3,18)
	coef_quad(18,3,1) = K(3,18)	! 3A3D + 3A -> boundary -> 1A
	coef_quad(3,18,1) = K(3,18)
	coef_quad(18,4,19) = K(4,18)	! 3A3D + 4A -> boundary -> 4A3D
	coef_quad(4,18,19) = K(4,18)
	coef_quad(18,4,1) = K(4,18)	! 3A3D + 4A -> boundary -> 1A
	coef_quad(4,18,1) = K(4,18)
	coef_quad(18,5,23) = K(5,18)	! 3A3D + 1D -> 3A4D
	coef_quad(5,18,23) = K(5,18)
	coef_lin(18,5,23) = E(5,18)	! 3A4D -> 3A3D + 1D
	coef_lin(5,18,23) = E(5,18)
	coef_quad(18,6,24) = K(6,18)	! 3A3D + 1A1D -> 4A4D
	coef_quad(6,18,24) = K(6,18)
	coef_lin(18,6,24) = E(6,18)	! 4A4D -> 3A3D + 1A1D
	coef_lin(6,18,24) = E(6,18)
	coef_quad(18,7,30) = K(7,18)	! 3A3D + 2A1D -> out_neu
	coef_quad(7,18,30) = K(7,18)
	coef_quad(18,8,30) = K(8,18)	! 3A3D + 3A1D -> out_neu
	coef_quad(8,18,30) = K(8,18)
	coef_quad(18,9,30) = K(9,18)	! 3A3D + 4A1D -> out_neu
	coef_quad(9,18,30) = K(9,18)
	coef_quad(18,10,23) = K(10,18)	! 3A3D + 2D -> boundary -> 3A4D
	coef_quad(10,18,23) = K(10,18)
	coef_quad(18,10,5) = K(10,18)	! 3A3D + 2D -> boundary -> 1D
	coef_quad(10,18,5) = K(10,18)
	coef_quad(18,11,24) = K(11,18)	! 3A3D + 1A2D -> boundary -> 4A4D
	coef_quad(11,18,24) = K(11,18)
	coef_quad(18,11,5) = K(11,18)	! 3A3D + 1A2D -> boundary -> 1D
	coef_quad(11,18,5) = K(11,18)
	coef_quad(18,12,30) = K(12,18)	! 3A3D + 2A2D -> out_neu
	coef_quad(12,18,30) = K(12,18)
	coef_quad(18,13,30) = K(13,18)	! 3A3D + 3A2D -> out_neu
	coef_quad(13,18,30) = K(13,18)
	coef_quad(18,14,30) = K(14,18)	! 3A3D + 4A2D -> out_neu
	coef_quad(14,18,30) = K(14,18)
	coef_quad(18,15,23) = K(15,18)	! 3A3D + 3D -> boundary -> 3A4D
	coef_quad(15,18,23) = K(15,18)
	coef_quad(18,15,5) = K(15,18)	! 3A3D + 3D -> boundary -> 1D
	coef_quad(15,18,5) = K(15,18)
	coef_quad(18,16,24) = K(16,18)	! 3A3D + 1A3D -> boundary -> 4A4D
	coef_quad(16,18,24) = K(16,18)
	coef_quad(18,16,5) = K(16,18)	! 3A3D + 1A3D -> boundary -> 1D
	coef_quad(16,18,5) = K(16,18)
	coef_quad(18,17,30) = K(17,18)	! 3A3D + 2A3D -> out_neu
	coef_quad(17,18,30) = K(17,18)
	coef_quad(18,18,30) = 0.5d0*K(18,18)	! 3A3D + 3A3D -> out_neu

	coef_quad(19,2,19) = K(2,19)	! 4A3D + 2A -> boundary -> 4A3D
	coef_quad(2,19,19) = K(2,19)
	coef_quad(19,2,1) = K(2,19)	! 4A3D + 2A -> boundary -> 1A
	coef_quad(2,19,1) = K(2,19)
	coef_quad(19,3,19) = K(3,19)	! 4A3D + 3A -> boundary -> 4A3D
	coef_quad(3,19,19) = K(3,19)
	coef_quad(19,3,1) = K(3,19)	! 4A3D + 3A -> boundary -> 1A
	coef_quad(3,19,1) = K(3,19)
	coef_quad(19,4,19) = K(4,19)	! 4A3D + 4A -> boundary -> 4A3D
	coef_quad(4,19,19) = K(4,19)
	coef_quad(19,4,1) = K(4,19)	! 4A3D + 4A -> boundary -> 1A
	coef_quad(4,19,1) = K(4,19)
	coef_quad(19,5,24) = K(5,19)	! 4A3D + 1D -> 4A4D
	coef_quad(5,19,24) = K(5,19)
	coef_lin(19,5,24) = E(5,19)	! 4A4D -> 4A3D + 1D
	coef_lin(5,19,24) = E(5,19)
	coef_quad(19,6,30) = K(6,19)	! 4A3D + 1A1D -> out_neu
	coef_quad(6,19,30) = K(6,19)
	coef_quad(19,7,30) = K(7,19)	! 4A3D + 2A1D -> out_neu
	coef_quad(7,19,30) = K(7,19)
	coef_quad(19,8,30) = K(8,19)	! 4A3D + 3A1D -> out_neu
	coef_quad(8,19,30) = K(8,19)
	coef_quad(19,9,30) = K(9,19)	! 4A3D + 4A1D -> out_neu
	coef_quad(9,19,30) = K(9,19)
	coef_quad(19,10,24) = K(10,19)	! 4A3D + 2D -> boundary -> 4A4D
	coef_quad(10,19,24) = K(10,19)
	coef_quad(19,10,5) = K(10,19)	! 4A3D + 2D -> boundary -> 1D
	coef_quad(10,19,5) = K(10,19)
	coef_quad(19,11,30) = K(11,19)	! 4A3D + 1A2D -> out_neu
	coef_quad(11,19,30) = K(11,19)
	coef_quad(19,12,30) = K(12,19)	! 4A3D + 2A2D -> out_neu
	coef_quad(12,19,30) = K(12,19)
	coef_quad(19,13,30) = K(13,19)	! 4A3D + 3A2D -> out_neu
	coef_quad(13,19,30) = K(13,19)
	coef_quad(19,14,30) = K(14,19)	! 4A3D + 4A2D -> out_neu
	coef_quad(14,19,30) = K(14,19)
	coef_quad(19,15,24) = K(15,19)	! 4A3D + 3D -> boundary -> 4A4D
	coef_quad(15,19,24) = K(15,19)
	coef_quad(19,15,5) = K(15,19)	! 4A3D + 3D -> boundary -> 1D
	coef_quad(15,19,5) = K(15,19)
	coef_quad(19,16,30) = K(16,19)	! 4A3D + 1A3D -> out_neu
	coef_quad(16,19,30) = K(16,19)
	coef_quad(19,17,30) = K(17,19)	! 4A3D + 2A3D -> out_neu
	coef_quad(17,19,30) = K(17,19)
	coef_quad(19,18,30) = K(18,19)	! 4A3D + 3A3D -> out_neu
	coef_quad(18,19,30) = K(18,19)
	coef_quad(19,19,30) = 0.5d0*K(19,19)	! 4A3D + 4A3D -> out_neu

	coef_quad(20,1,21) = K(1,20)	! 4D + 1A -> 1A4D
	coef_quad(1,20,21) = K(1,20)
	coef_lin(20,1,21) = E(1,20)	! 1A4D -> 4D + 1A
	coef_lin(1,20,21) = E(1,20)
	coef_quad(20,2,22) = K(2,20)	! 4D + 2A -> 2A4D
	coef_quad(2,20,22) = K(2,20)
	coef_lin(20,2,22) = E(2,20)	! 2A4D -> 4D + 2A
	coef_lin(2,20,22) = E(2,20)
	coef_quad(20,3,23) = K(3,20)	! 4D + 3A -> 3A4D
	coef_quad(3,20,23) = K(3,20)
	coef_lin(20,3,23) = E(3,20)	! 3A4D -> 4D + 3A
	coef_lin(3,20,23) = E(3,20)
	coef_quad(20,4,24) = K(4,20)	! 4D + 4A -> 4A4D
	coef_quad(4,20,24) = K(4,20)
	coef_lin(20,4,24) = E(4,20)	! 4A4D -> 4D + 4A
	coef_lin(4,20,24) = E(4,20)
	coef_quad(20,6,21) = K(6,20)	! 4D + 1A1D -> boundary -> 1A4D
	coef_quad(6,20,21) = K(6,20)
	coef_quad(20,6,5) = K(6,20)	! 4D + 1A1D -> boundary -> 1D
	coef_quad(6,20,5) = K(6,20)
	coef_quad(20,7,22) = K(7,20)	! 4D + 2A1D -> boundary -> 2A4D
	coef_quad(7,20,22) = K(7,20)
	coef_quad(20,7,5) = K(7,20)	! 4D + 2A1D -> boundary -> 1D
	coef_quad(7,20,5) = K(7,20)
	coef_quad(20,8,23) = K(8,20)	! 4D + 3A1D -> boundary -> 3A4D
	coef_quad(8,20,23) = K(8,20)
	coef_quad(20,8,5) = K(8,20)	! 4D + 3A1D -> boundary -> 1D
	coef_quad(8,20,5) = K(8,20)
	coef_quad(20,9,24) = K(9,20)	! 4D + 4A1D -> boundary -> 4A4D
	coef_quad(9,20,24) = K(9,20)
	coef_quad(20,9,5) = K(9,20)	! 4D + 4A1D -> boundary -> 1D
	coef_quad(9,20,5) = K(9,20)
	coef_quad(20,10,20) = K(10,20)	! 4D + 2D -> boundary -> 4D
	coef_quad(10,20,20) = K(10,20)
	coef_quad(20,10,5) = K(10,20)	! 4D + 2D -> boundary -> 1D
	coef_quad(10,20,5) = K(10,20)
	coef_quad(20,11,21) = K(11,20)	! 4D + 1A2D -> boundary -> 1A4D
	coef_quad(11,20,21) = K(11,20)
	coef_quad(20,11,5) = K(11,20)	! 4D + 1A2D -> boundary -> 1D
	coef_quad(11,20,5) = K(11,20)
	coef_quad(20,12,22) = K(12,20)	! 4D + 2A2D -> boundary -> 2A4D
	coef_quad(12,20,22) = K(12,20)
	coef_quad(20,12,5) = K(12,20)	! 4D + 2A2D -> boundary -> 1D
	coef_quad(12,20,5) = K(12,20)
	coef_quad(20,13,23) = K(13,20)	! 4D + 3A2D -> boundary -> 3A4D
	coef_quad(13,20,23) = K(13,20)
	coef_quad(20,13,5) = K(13,20)	! 4D + 3A2D -> boundary -> 1D
	coef_quad(13,20,5) = K(13,20)
	coef_quad(20,14,24) = K(14,20)	! 4D + 4A2D -> boundary -> 4A4D
	coef_quad(14,20,24) = K(14,20)
	coef_quad(20,14,5) = K(14,20)	! 4D + 4A2D -> boundary -> 1D
	coef_quad(14,20,5) = K(14,20)
	coef_quad(20,15,20) = K(15,20)	! 4D + 3D -> boundary -> 4D
	coef_quad(15,20,20) = K(15,20)
	coef_quad(20,15,5) = K(15,20)	! 4D + 3D -> boundary -> 1D
	coef_quad(15,20,5) = K(15,20)
	coef_quad(20,16,21) = K(16,20)	! 4D + 1A3D -> boundary -> 1A4D
	coef_quad(16,20,21) = K(16,20)
	coef_quad(20,16,5) = K(16,20)	! 4D + 1A3D -> boundary -> 1D
	coef_quad(16,20,5) = K(16,20)
	coef_quad(20,17,22) = K(17,20)	! 4D + 2A3D -> boundary -> 2A4D
	coef_quad(17,20,22) = K(17,20)
	coef_quad(20,17,5) = K(17,20)	! 4D + 2A3D -> boundary -> 1D
	coef_quad(17,20,5) = K(17,20)
	coef_quad(20,18,23) = K(18,20)	! 4D + 3A3D -> boundary -> 3A4D
	coef_quad(18,20,23) = K(18,20)
	coef_quad(20,18,5) = K(18,20)	! 4D + 3A3D -> boundary -> 1D
	coef_quad(18,20,5) = K(18,20)
	coef_quad(20,19,24) = K(19,20)	! 4D + 4A3D -> boundary -> 4A4D
	coef_quad(19,20,24) = K(19,20)
	coef_quad(20,19,5) = K(19,20)	! 4D + 4A3D -> boundary -> 1D
	coef_quad(19,20,5) = K(19,20)
	coef_quad(20,20,20) = 0.5d0*K(20,20)	! 4D + 4D -> boundary -> 4D
	coef_quad(20,20,5) = 0.5d0*K(20,20)	! 4D + 4D -> boundary -> 1D

	coef_quad(21,1,22) = K(1,21)	! 1A4D + 1A -> 2A4D
	coef_quad(1,21,22) = K(1,21)
	coef_lin(21,1,22) = E(1,21)	! 2A4D -> 1A4D + 1A
	coef_lin(1,21,22) = E(1,21)
	coef_quad(21,2,23) = K(2,21)	! 1A4D + 2A -> 3A4D
	coef_quad(2,21,23) = K(2,21)
	coef_lin(21,2,23) = E(2,21)	! 3A4D -> 1A4D + 2A
	coef_lin(2,21,23) = E(2,21)
	coef_quad(21,3,24) = K(3,21)	! 1A4D + 3A -> 4A4D
	coef_quad(3,21,24) = K(3,21)
	coef_lin(21,3,24) = E(3,21)	! 4A4D -> 1A4D + 3A
	coef_lin(3,21,24) = E(3,21)
	coef_quad(21,4,30) = K(4,21)	! 1A4D + 4A -> out_neu
	coef_quad(4,21,30) = K(4,21)
	coef_quad(21,6,22) = K(6,21)	! 1A4D + 1A1D -> boundary -> 2A4D
	coef_quad(6,21,22) = K(6,21)
	coef_quad(21,6,5) = K(6,21)	! 1A4D + 1A1D -> boundary -> 1D
	coef_quad(6,21,5) = K(6,21)
	coef_quad(21,7,23) = K(7,21)	! 1A4D + 2A1D -> boundary -> 3A4D
	coef_quad(7,21,23) = K(7,21)
	coef_quad(21,7,5) = K(7,21)	! 1A4D + 2A1D -> boundary -> 1D
	coef_quad(7,21,5) = K(7,21)
	coef_quad(21,8,24) = K(8,21)	! 1A4D + 3A1D -> boundary -> 4A4D
	coef_quad(8,21,24) = K(8,21)
	coef_quad(21,8,5) = K(8,21)	! 1A4D + 3A1D -> boundary -> 1D
	coef_quad(8,21,5) = K(8,21)
	coef_quad(21,9,30) = K(9,21)	! 1A4D + 4A1D -> out_neu
	coef_quad(9,21,30) = K(9,21)
	coef_quad(21,10,21) = K(10,21)	! 1A4D + 2D -> boundary -> 1A4D
	coef_quad(10,21,21) = K(10,21)
	coef_quad(21,10,5) = K(10,21)	! 1A4D + 2D -> boundary -> 1D
	coef_quad(10,21,5) = K(10,21)
	coef_quad(21,11,22) = K(11,21)	! 1A4D + 1A2D -> boundary -> 2A4D
	coef_quad(11,21,22) = K(11,21)
	coef_quad(21,11,5) = K(11,21)	! 1A4D + 1A2D -> boundary -> 1D
	coef_quad(11,21,5) = K(11,21)
	coef_quad(21,12,23) = K(12,21)	! 1A4D + 2A2D -> boundary -> 3A4D
	coef_quad(12,21,23) = K(12,21)
	coef_quad(21,12,5) = K(12,21)	! 1A4D + 2A2D -> boundary -> 1D
	coef_quad(12,21,5) = K(12,21)
	coef_quad(21,13,24) = K(13,21)	! 1A4D + 3A2D -> boundary -> 4A4D
	coef_quad(13,21,24) = K(13,21)
	coef_quad(21,13,5) = K(13,21)	! 1A4D + 3A2D -> boundary -> 1D
	coef_quad(13,21,5) = K(13,21)
	coef_quad(21,14,30) = K(14,21)	! 1A4D + 4A2D -> out_neu
	coef_quad(14,21,30) = K(14,21)
	coef_quad(21,15,21) = K(15,21)	! 1A4D + 3D -> boundary -> 1A4D
	coef_quad(15,21,21) = K(15,21)
	coef_quad(21,15,5) = K(15,21)	! 1A4D + 3D -> boundary -> 1D
	coef_quad(15,21,5) = K(15,21)
	coef_quad(21,16,22) = K(16,21)	! 1A4D + 1A3D -> boundary -> 2A4D
	coef_quad(16,21,22) = K(16,21)
	coef_quad(21,16,5) = K(16,21)	! 1A4D + 1A3D -> boundary -> 1D
	coef_quad(16,21,5) = K(16,21)
	coef_quad(21,17,23) = K(17,21)	! 1A4D + 2A3D -> boundary -> 3A4D
	coef_quad(17,21,23) = K(17,21)
	coef_quad(21,17,5) = K(17,21)	! 1A4D + 2A3D -> boundary -> 1D
	coef_quad(17,21,5) = K(17,21)
	coef_quad(21,18,24) = K(18,21)	! 1A4D + 3A3D -> boundary -> 4A4D
	coef_quad(18,21,24) = K(18,21)
	coef_quad(21,18,5) = K(18,21)	! 1A4D + 3A3D -> boundary -> 1D
	coef_quad(18,21,5) = K(18,21)
	coef_quad(21,19,30) = K(19,21)	! 1A4D + 4A3D -> out_neu
	coef_quad(19,21,30) = K(19,21)
	coef_quad(21,20,21) = K(20,21)	! 1A4D + 4D -> boundary -> 1A4D
	coef_quad(20,21,21) = K(20,21)
	coef_quad(21,20,5) = K(20,21)	! 1A4D + 4D -> boundary -> 1D
	coef_quad(20,21,5) = K(20,21)
	coef_quad(21,21,22) = 0.5d0*K(21,21)	! 1A4D + 1A4D -> boundary -> 2A4D
	coef_quad(21,21,5) = 0.5d0*K(21,21)	! 1A4D + 1A4D -> boundary -> 1D

	coef_quad(22,1,23) = K(1,22)	! 2A4D + 1A -> 3A4D
	coef_quad(1,22,23) = K(1,22)
	coef_lin(22,1,23) = E(1,22)	! 3A4D -> 2A4D + 1A
	coef_lin(1,22,23) = E(1,22)
	coef_quad(22,2,24) = K(2,22)	! 2A4D + 2A -> 4A4D
	coef_quad(2,22,24) = K(2,22)
	coef_lin(22,2,24) = E(2,22)	! 4A4D -> 2A4D + 2A
	coef_lin(2,22,24) = E(2,22)
	coef_quad(22,3,30) = K(3,22)	! 2A4D + 3A -> out_neu
	coef_quad(3,22,30) = K(3,22)
	coef_quad(22,4,30) = K(4,22)	! 2A4D + 4A -> out_neu
	coef_quad(4,22,30) = K(4,22)
	coef_quad(22,6,23) = K(6,22)	! 2A4D + 1A1D -> boundary -> 3A4D
	coef_quad(6,22,23) = K(6,22)
	coef_quad(22,6,5) = K(6,22)	! 2A4D + 1A1D -> boundary -> 1D
	coef_quad(6,22,5) = K(6,22)
	coef_quad(22,7,24) = K(7,22)	! 2A4D + 2A1D -> boundary -> 4A4D
	coef_quad(7,22,24) = K(7,22)
	coef_quad(22,7,5) = K(7,22)	! 2A4D + 2A1D -> boundary -> 1D
	coef_quad(7,22,5) = K(7,22)
	coef_quad(22,8,30) = K(8,22)	! 2A4D + 3A1D -> out_neu
	coef_quad(8,22,30) = K(8,22)
	coef_quad(22,9,30) = K(9,22)	! 2A4D + 4A1D -> out_neu
	coef_quad(9,22,30) = K(9,22)
	coef_quad(22,10,22) = K(10,22)	! 2A4D + 2D -> boundary -> 2A4D
	coef_quad(10,22,22) = K(10,22)
	coef_quad(22,10,5) = K(10,22)	! 2A4D + 2D -> boundary -> 1D
	coef_quad(10,22,5) = K(10,22)
	coef_quad(22,11,23) = K(11,22)	! 2A4D + 1A2D -> boundary -> 3A4D
	coef_quad(11,22,23) = K(11,22)
	coef_quad(22,11,5) = K(11,22)	! 2A4D + 1A2D -> boundary -> 1D
	coef_quad(11,22,5) = K(11,22)
	coef_quad(22,12,24) = K(12,22)	! 2A4D + 2A2D -> boundary -> 4A4D
	coef_quad(12,22,24) = K(12,22)
	coef_quad(22,12,5) = K(12,22)	! 2A4D + 2A2D -> boundary -> 1D
	coef_quad(12,22,5) = K(12,22)
	coef_quad(22,13,30) = K(13,22)	! 2A4D + 3A2D -> out_neu
	coef_quad(13,22,30) = K(13,22)
	coef_quad(22,14,30) = K(14,22)	! 2A4D + 4A2D -> out_neu
	coef_quad(14,22,30) = K(14,22)
	coef_quad(22,15,22) = K(15,22)	! 2A4D + 3D -> boundary -> 2A4D
	coef_quad(15,22,22) = K(15,22)
	coef_quad(22,15,5) = K(15,22)	! 2A4D + 3D -> boundary -> 1D
	coef_quad(15,22,5) = K(15,22)
	coef_quad(22,16,23) = K(16,22)	! 2A4D + 1A3D -> boundary -> 3A4D
	coef_quad(16,22,23) = K(16,22)
	coef_quad(22,16,5) = K(16,22)	! 2A4D + 1A3D -> boundary -> 1D
	coef_quad(16,22,5) = K(16,22)
	coef_quad(22,17,24) = K(17,22)	! 2A4D + 2A3D -> boundary -> 4A4D
	coef_quad(17,22,24) = K(17,22)
	coef_quad(22,17,5) = K(17,22)	! 2A4D + 2A3D -> boundary -> 1D
	coef_quad(17,22,5) = K(17,22)
	coef_quad(22,18,30) = K(18,22)	! 2A4D + 3A3D -> out_neu
	coef_quad(18,22,30) = K(18,22)
	coef_quad(22,19,30) = K(19,22)	! 2A4D + 4A3D -> out_neu
	coef_quad(19,22,30) = K(19,22)
	coef_quad(22,20,22) = K(20,22)	! 2A4D + 4D -> boundary -> 2A4D
	coef_quad(20,22,22) = K(20,22)
	coef_quad(22,20,5) = K(20,22)	! 2A4D + 4D -> boundary -> 1D
	coef_quad(20,22,5) = K(20,22)
	coef_quad(22,21,23) = K(21,22)	! 2A4D + 1A4D -> boundary -> 3A4D
	coef_quad(21,22,23) = K(21,22)
	coef_quad(22,21,5) = K(21,22)	! 2A4D + 1A4D -> boundary -> 1D
	coef_quad(21,22,5) = K(21,22)
	coef_quad(22,22,24) = 0.5d0*K(22,22)	! 2A4D + 2A4D -> boundary -> 4A4D
	coef_quad(22,22,5) = 0.5d0*K(22,22)	! 2A4D + 2A4D -> boundary -> 1D

	coef_quad(23,1,24) = K(1,23)	! 3A4D + 1A -> 4A4D
	coef_quad(1,23,24) = K(1,23)
	coef_lin(23,1,24) = E(1,23)	! 4A4D -> 3A4D + 1A
	coef_lin(1,23,24) = E(1,23)
	coef_quad(23,2,30) = K(2,23)	! 3A4D + 2A -> out_neu
	coef_quad(2,23,30) = K(2,23)
	coef_quad(23,3,30) = K(3,23)	! 3A4D + 3A -> out_neu
	coef_quad(3,23,30) = K(3,23)
	coef_quad(23,4,30) = K(4,23)	! 3A4D + 4A -> out_neu
	coef_quad(4,23,30) = K(4,23)
	coef_quad(23,6,24) = K(6,23)	! 3A4D + 1A1D -> boundary -> 4A4D
	coef_quad(6,23,24) = K(6,23)
	coef_quad(23,6,5) = K(6,23)	! 3A4D + 1A1D -> boundary -> 1D
	coef_quad(6,23,5) = K(6,23)
	coef_quad(23,7,30) = K(7,23)	! 3A4D + 2A1D -> out_neu
	coef_quad(7,23,30) = K(7,23)
	coef_quad(23,8,30) = K(8,23)	! 3A4D + 3A1D -> out_neu
	coef_quad(8,23,30) = K(8,23)
	coef_quad(23,9,30) = K(9,23)	! 3A4D + 4A1D -> out_neu
	coef_quad(9,23,30) = K(9,23)
	coef_quad(23,10,23) = K(10,23)	! 3A4D + 2D -> boundary -> 3A4D
	coef_quad(10,23,23) = K(10,23)
	coef_quad(23,10,5) = K(10,23)	! 3A4D + 2D -> boundary -> 1D
	coef_quad(10,23,5) = K(10,23)
	coef_quad(23,11,24) = K(11,23)	! 3A4D + 1A2D -> boundary -> 4A4D
	coef_quad(11,23,24) = K(11,23)
	coef_quad(23,11,5) = K(11,23)	! 3A4D + 1A2D -> boundary -> 1D
	coef_quad(11,23,5) = K(11,23)
	coef_quad(23,12,30) = K(12,23)	! 3A4D + 2A2D -> out_neu
	coef_quad(12,23,30) = K(12,23)
	coef_quad(23,13,30) = K(13,23)	! 3A4D + 3A2D -> out_neu
	coef_quad(13,23,30) = K(13,23)
	coef_quad(23,14,30) = K(14,23)	! 3A4D + 4A2D -> out_neu
	coef_quad(14,23,30) = K(14,23)
	coef_quad(23,15,23) = K(15,23)	! 3A4D + 3D -> boundary -> 3A4D
	coef_quad(15,23,23) = K(15,23)
	coef_quad(23,15,5) = K(15,23)	! 3A4D + 3D -> boundary -> 1D
	coef_quad(15,23,5) = K(15,23)
	coef_quad(23,16,24) = K(16,23)	! 3A4D + 1A3D -> boundary -> 4A4D
	coef_quad(16,23,24) = K(16,23)
	coef_quad(23,16,5) = K(16,23)	! 3A4D + 1A3D -> boundary -> 1D
	coef_quad(16,23,5) = K(16,23)
	coef_quad(23,17,30) = K(17,23)	! 3A4D + 2A3D -> out_neu
	coef_quad(17,23,30) = K(17,23)
	coef_quad(23,18,30) = K(18,23)	! 3A4D + 3A3D -> out_neu
	coef_quad(18,23,30) = K(18,23)
	coef_quad(23,19,30) = K(19,23)	! 3A4D + 4A3D -> out_neu
	coef_quad(19,23,30) = K(19,23)
	coef_quad(23,20,23) = K(20,23)	! 3A4D + 4D -> boundary -> 3A4D
	coef_quad(20,23,23) = K(20,23)
	coef_quad(23,20,5) = K(20,23)	! 3A4D + 4D -> boundary -> 1D
	coef_quad(20,23,5) = K(20,23)
	coef_quad(23,21,24) = K(21,23)	! 3A4D + 1A4D -> boundary -> 4A4D
	coef_quad(21,23,24) = K(21,23)
	coef_quad(23,21,5) = K(21,23)	! 3A4D + 1A4D -> boundary -> 1D
	coef_quad(21,23,5) = K(21,23)
	coef_quad(23,22,30) = K(22,23)	! 3A4D + 2A4D -> out_neu
	coef_quad(22,23,30) = K(22,23)
	coef_quad(23,23,30) = 0.5d0*K(23,23)	! 3A4D + 3A4D -> out_neu

	coef_quad(24,1,30) = K(1,24)	! 4A4D + 1A -> out_neu
	coef_quad(1,24,30) = K(1,24)
	coef_quad(24,2,30) = K(2,24)	! 4A4D + 2A -> out_neu
	coef_quad(2,24,30) = K(2,24)
	coef_quad(24,3,30) = K(3,24)	! 4A4D + 3A -> out_neu
	coef_quad(3,24,30) = K(3,24)
	coef_quad(24,4,30) = K(4,24)	! 4A4D + 4A -> out_neu
	coef_quad(4,24,30) = K(4,24)
	coef_quad(24,6,30) = K(6,24)	! 4A4D + 1A1D -> out_neu
	coef_quad(6,24,30) = K(6,24)
	coef_quad(24,7,30) = K(7,24)	! 4A4D + 2A1D -> out_neu
	coef_quad(7,24,30) = K(7,24)
	coef_quad(24,8,30) = K(8,24)	! 4A4D + 3A1D -> out_neu
	coef_quad(8,24,30) = K(8,24)
	coef_quad(24,9,30) = K(9,24)	! 4A4D + 4A1D -> out_neu
	coef_quad(9,24,30) = K(9,24)
	coef_quad(24,10,24) = K(10,24)	! 4A4D + 2D -> boundary -> 4A4D
	coef_quad(10,24,24) = K(10,24)
	coef_quad(24,10,5) = K(10,24)	! 4A4D + 2D -> boundary -> 1D
	coef_quad(10,24,5) = K(10,24)
	coef_quad(24,11,30) = K(11,24)	! 4A4D + 1A2D -> out_neu
	coef_quad(11,24,30) = K(11,24)
	coef_quad(24,12,30) = K(12,24)	! 4A4D + 2A2D -> out_neu
	coef_quad(12,24,30) = K(12,24)
	coef_quad(24,13,30) = K(13,24)	! 4A4D + 3A2D -> out_neu
	coef_quad(13,24,30) = K(13,24)
	coef_quad(24,14,30) = K(14,24)	! 4A4D + 4A2D -> out_neu
	coef_quad(14,24,30) = K(14,24)
	coef_quad(24,15,24) = K(15,24)	! 4A4D + 3D -> boundary -> 4A4D
	coef_quad(15,24,24) = K(15,24)
	coef_quad(24,15,5) = K(15,24)	! 4A4D + 3D -> boundary -> 1D
	coef_quad(15,24,5) = K(15,24)
	coef_quad(24,16,30) = K(16,24)	! 4A4D + 1A3D -> out_neu
	coef_quad(16,24,30) = K(16,24)
	coef_quad(24,17,30) = K(17,24)	! 4A4D + 2A3D -> out_neu
	coef_quad(17,24,30) = K(17,24)
	coef_quad(24,18,30) = K(18,24)	! 4A4D + 3A3D -> out_neu
	coef_quad(18,24,30) = K(18,24)
	coef_quad(24,19,30) = K(19,24)	! 4A4D + 4A3D -> out_neu
	coef_quad(19,24,30) = K(19,24)
	coef_quad(24,20,24) = K(20,24)	! 4A4D + 4D -> boundary -> 4A4D
	coef_quad(20,24,24) = K(20,24)
	coef_quad(24,20,5) = K(20,24)	! 4A4D + 4D -> boundary -> 1D
	coef_quad(20,24,5) = K(20,24)
	coef_quad(24,21,30) = K(21,24)	! 4A4D + 1A4D -> out_neu
	coef_quad(21,24,30) = K(21,24)
	coef_quad(24,22,30) = K(22,24)	! 4A4D + 2A4D -> out_neu
	coef_quad(22,24,30) = K(22,24)
	coef_quad(24,23,30) = K(23,24)	! 4A4D + 3A4D -> out_neu
	coef_quad(23,24,30) = K(23,24)
	coef_quad(24,24,30) = 0.5d0*K(24,24)	! 4A4D + 4A4D -> out_neu

	coef_lin(26,26,24) = cs(24)	! 4A4D -> coag
	coef_lin(26,26,23) = cs(23)	! 3A4D -> coag
	coef_lin(26,26,22) = cs(22)	! 2A4D -> coag
	coef_lin(26,26,21) = cs(21)	! 1A4D -> coag
	coef_lin(26,26,20) = cs(20)	! 4D -> coag
	coef_lin(26,26,19) = cs(19)	! 4A3D -> coag
	coef_lin(26,26,18) = cs(18)	! 3A3D -> coag
	coef_lin(26,26,17) = cs(17)	! 2A3D -> coag
	coef_lin(26,26,16) = cs(16)	! 1A3D -> coag
	coef_lin(26,26,15) = cs(15)	! 3D -> coag
	coef_lin(26,26,14) = cs(14)	! 4A2D -> coag
	coef_lin(26,26,13) = cs(13)	! 3A2D -> coag
	coef_lin(26,26,12) = cs(12)	! 2A2D -> coag
	coef_lin(26,26,11) = cs(11)	! 1A2D -> coag
	coef_lin(26,26,10) = cs(10)	! 2D -> coag
	coef_lin(26,26,9) = cs(9)	! 4A1D -> coag
	coef_lin(26,26,8) = cs(8)	! 3A1D -> coag
	coef_lin(26,26,7) = cs(7)	! 2A1D -> coag
	coef_lin(26,26,6) = cs(6)	! 1A1D -> coag
	coef_lin(26,26,5) = cs(5)	! 1D -> coag
	coef_lin(26,26,4) = cs(4)	! 4A -> coag
	coef_lin(26,26,3) = cs(3)	! 3A -> coag
	coef_lin(26,26,2) = cs(2)	! 2A -> coag
	coef_lin(26,26,1) = cs(1)	! 1A -> coag

end subroutine get_rate_coefsD

!-----------------------------------------------------------

subroutine get_lossesD(cs)
	implicit none
	real(kind(1.d0)) :: cs(24)

	! coagulation sink

	! variable cs: the following values only give the size dependence of the sink,
	! and will be multiplied by the size-independent factor given as input

	cs(1) = 1.00000000000000d+00	! coagulation loss of 1A
	cs(2) = 6.90956439983888d-01	! coagulation loss of 2A
	cs(3) = 5.56589912236292d-01	! coagulation loss of 3A
	cs(4) = 4.77420801955208d-01	! coagulation loss of 4A
	cs(5) = 8.92789345574936d-01	! coagulation loss of 1D
	cs(6) = 6.50906449049502d-01	! coagulation loss of 1A1D
	cs(7) = 5.34476727643759d-01	! coagulation loss of 2A1D
	cs(8) = 4.62991019676722d-01	! coagulation loss of 3A1D
	cs(9) = 4.13516007807625d-01	! coagulation loss of 4A1D
	cs(10) = 6.16878547874003d-01	! coagulation loss of 2D
	cs(11) = 5.14714704412136d-01	! coagulation loss of 1A2D
	cs(12) = 4.49748002797798d-01	! coagulation loss of 2A2D
	cs(13) = 4.03871535975795d-01	! coagulation loss of 3A2D
	cs(14) = 3.69300136986970d-01	! coagulation loss of 4A2D
	cs(15) = 4.96917543499051d-01	! coagulation loss of 3D
	cs(16) = 4.37538685504830d-01	! coagulation loss of 1A3D
	cs(17) = 3.94846710007348d-01	! coagulation loss of 2A3D
	cs(18) = 3.62287963350499d-01	! coagulation loss of 3A3D
	cs(19) = 3.36418964945053d-01	! coagulation loss of 4A3D
	cs(20) = 4.26236205341452d-01	! coagulation loss of 4D
	cs(21) = 3.86378395643578d-01	! coagulation loss of 1A4D
	cs(22) = 3.55645439767969d-01	! coagulation loss of 2A4D
	cs(23) = 3.31032373593705d-01	! coagulation loss of 3A4D
	cs(24) = 3.10756278903030d-01	! coagulation loss of 4A4D

end subroutine get_lossesD

!-----------------------------------------------------------

subroutine get_collD(K,temperature)
	implicit none
	integer, parameter :: nclust = 24
	real(kind(1.d0)) :: K(nclust,nclust), temperature

	! collision coefficients

	K = 0.d0
	K(1,1) = 2.00300435981201d-17*sqrt(temperature)	! 1A + 1A
	K(2,1) = 2.21482322896243d-17*sqrt(temperature)	! 2A + 1A
	K(1,2) = K(2,1)
	K(2,2) = 2.24829637648719d-17*sqrt(temperature)	! 2A + 2A
	K(3,1) = 2.43868865258245d-17*sqrt(temperature)	! 3A + 1A
	K(1,3) = K(3,1)
	K(3,2) = 2.36016202277657d-17*sqrt(temperature)	! 3A + 2A
	K(2,3) = K(3,2)
	K(3,3) = 2.40548195707689d-17*sqrt(temperature)	! 3A + 3A
	K(4,2) = 2.48605501822435d-17*sqrt(temperature)	! 4A + 2A
	K(2,4) = K(4,2)
	K(4,3) = 2.48227789571021d-17*sqrt(temperature)	! 4A + 3A
	K(3,4) = K(4,3)
	K(4,4) = 2.52362735595835d-17*sqrt(temperature)	! 4A + 4A
	K(5,1) = 2.71277008815771d-17*sqrt(temperature)	! 1D + 1A
	K(1,5) = K(5,1)
	K(5,2) = 3.15349593390535d-17*sqrt(temperature)	! 1D + 2A
	K(2,5) = K(5,2)
	K(5,3) = 3.54957741166103d-17*sqrt(temperature)	! 1D + 3A
	K(3,5) = K(5,3)
	K(5,4) = 3.90450082642634d-17*sqrt(temperature)	! 1D + 4A
	K(4,5) = K(5,4)
	K(5,5) = 3.40442617531518d-17*sqrt(temperature)	! 1D + 1D
	K(6,1) = 2.44809288262135d-17*sqrt(temperature)	! 1A1D + 1A
	K(1,6) = K(6,1)
	K(6,2) = 2.54150641906904d-17*sqrt(temperature)	! 1A1D + 2A
	K(2,6) = K(6,2)
	K(6,3) = 2.70249536387104d-17*sqrt(temperature)	! 1A1D + 3A
	K(3,6) = K(6,3)
	K(6,4) = 2.87014463732710d-17*sqrt(temperature)	! 1A1D + 4A
	K(4,6) = K(6,4)
	K(6,5) = 3.39602459742711d-17*sqrt(temperature)	! 1A1D + 1D
	K(5,6) = K(6,5)
	K(6,6) = 2.83571813260610d-17*sqrt(temperature)	! 1A1D + 1A1D
	K(7,1) = 2.58125616100570d-17*sqrt(temperature)	! 2A1D + 1A
	K(1,7) = K(7,1)
	K(7,2) = 2.52957813819825d-17*sqrt(temperature)	! 2A1D + 2A
	K(2,7) = K(7,2)
	K(7,3) = 2.59960796107242d-17*sqrt(temperature)	! 2A1D + 3A
	K(3,7) = K(7,3)
	K(7,4) = 2.69822192354348d-17*sqrt(temperature)	! 2A1D + 4A
	K(4,7) = K(7,4)
	K(7,5) = 3.70774069121475d-17*sqrt(temperature)	! 2A1D + 1D
	K(5,7) = K(7,5)
	K(7,6) = 2.87378921679481d-17*sqrt(temperature)	! 2A1D + 1A1D
	K(6,7) = K(7,6)
	K(7,7) = 2.79470816379128d-17*sqrt(temperature)	! 2A1D + 2A1D
	K(8,1) = 2.75567881474783d-17*sqrt(temperature)	! 3A1D + 1A
	K(1,8) = K(8,1)
	K(8,2) = 2.60531351555787d-17*sqrt(temperature)	! 3A1D + 2A
	K(2,8) = K(8,2)
	K(8,3) = 2.61631458720229d-17*sqrt(temperature)	! 3A1D + 3A
	K(3,8) = K(8,3)
	K(8,4) = 2.67131220926838d-17*sqrt(temperature)	! 3A1D + 4A
	K(4,8) = K(8,4)
	K(8,5) = 4.02727477268666d-17*sqrt(temperature)	! 3A1D + 1D
	K(5,8) = K(8,5)
	K(8,6) = 2.99200115944640d-17*sqrt(temperature)	! 3A1D + 1A1D
	K(6,8) = K(8,6)
	K(8,7) = 2.83326354588116d-17*sqrt(temperature)	! 3A1D + 2A1D
	K(7,8) = K(8,7)
	K(8,8) = 2.81969208207520d-17*sqrt(temperature)	! 3A1D + 3A1D
	K(9,2) = 2.70565466746347d-17*sqrt(temperature)	! 4A1D + 2A
	K(2,9) = K(9,2)
	K(9,3) = 2.67169801300549d-17*sqrt(temperature)	! 4A1D + 3A
	K(3,9) = K(9,3)
	K(9,4) = 2.69412515898993d-17*sqrt(temperature)	! 4A1D + 4A
	K(4,9) = K(9,4)
	K(9,5) = 4.33136033685207d-17*sqrt(temperature)	! 4A1D + 1D
	K(5,9) = K(9,5)
	K(9,6) = 3.12948976926336d-17*sqrt(temperature)	! 4A1D + 1A1D
	K(6,9) = K(9,6)
	K(9,7) = 2.90829841036644d-17*sqrt(temperature)	! 4A1D + 2A1D
	K(7,9) = K(9,7)
	K(9,8) = 2.85493511684281d-17*sqrt(temperature)	! 4A1D + 3A1D
	K(8,9) = K(9,8)
	K(9,9) = 2.86033780268426d-17*sqrt(temperature)	! 4A1D + 4A1D
	K(10,1) = 2.83140047857577d-17*sqrt(temperature)	! 2D + 1A
	K(1,10) = K(10,1)
	K(10,2) = 3.04498146974590d-17*sqrt(temperature)	! 2D + 2A
	K(2,10) = K(10,2)
	K(10,3) = 3.29689513431534d-17*sqrt(temperature)	! 2D + 3A
	K(3,10) = K(10,3)
	K(10,4) = 3.53967950530668d-17*sqrt(temperature)	! 2D + 4A
	K(4,10) = K(10,4)
	K(10,5) = 3.76444621173140d-17*sqrt(temperature)	! 2D + 1D
	K(5,10) = K(10,5)
	K(10,6) = 3.33667324616888d-17*sqrt(temperature)	! 2D + 1A1D
	K(6,10) = K(10,6)
	K(10,7) = 3.47090814121944d-17*sqrt(temperature)	! 2D + 2A1D
	K(7,10) = K(10,7)
	K(10,8) = 3.66646153713648d-17*sqrt(temperature)	! 2D + 3A1D
	K(8,10) = K(10,8)
	K(10,9) = 3.87018440610044d-17*sqrt(temperature)	! 2D + 4A1D
	K(9,10) = K(10,9)
	K(10,10) = 3.82133917806232d-17*sqrt(temperature)	! 2D + 2D
	K(11,1) = 2.76109922240466d-17*sqrt(temperature)	! 1A2D + 1A
	K(1,11) = K(11,1)
	K(11,2) = 2.75406472349300d-17*sqrt(temperature)	! 1A2D + 2A
	K(2,11) = K(11,2)
	K(11,3) = 2.86128234794644d-17*sqrt(temperature)	! 1A2D + 3A
	K(3,11) = K(11,3)
	K(11,4) = 2.99158924119007d-17*sqrt(temperature)	! 1A2D + 4A
	K(4,11) = K(11,4)
	K(11,5) = 3.89438517595401d-17*sqrt(temperature)	! 1A2D + 1D
	K(5,11) = K(11,5)
	K(11,6) = 3.09758708170216d-17*sqrt(temperature)	! 1A2D + 1A1D
	K(6,11) = K(11,6)
	K(11,7) = 3.05649415958920d-17*sqrt(temperature)	! 1A2D + 2A1D
	K(7,11) = K(11,7)
	K(11,8) = 3.12749051699915d-17*sqrt(temperature)	! 1A2D + 3A1D
	K(8,11) = K(11,8)
	K(11,9) = 3.23086035271087d-17*sqrt(temperature)	! 1A2D + 4A1D
	K(9,11) = K(11,9)
	K(11,10) = 3.69159922702319d-17*sqrt(temperature)	! 1A2D + 2D
	K(10,11) = K(11,10)
	K(11,11) = 3.31633158250450d-17*sqrt(temperature)	! 1A2D + 1A2D
	K(12,1) = 2.87627441579623d-17*sqrt(temperature)	! 2A2D + 1A
	K(1,12) = K(12,1)
	K(12,2) = 2.74789135147876d-17*sqrt(temperature)	! 2A2D + 2A
	K(2,12) = K(12,2)
	K(12,3) = 2.77948237234039d-17*sqrt(temperature)	! 2A2D + 3A
	K(3,12) = K(12,3)
	K(12,4) = 2.85274450093966d-17*sqrt(temperature)	! 2A2D + 4A
	K(4,12) = K(12,4)
	K(12,5) = 4.16083096768037d-17*sqrt(temperature)	! 2A2D + 1D
	K(5,12) = K(12,5)
	K(12,6) = 3.13559696468930d-17*sqrt(temperature)	! 2A2D + 1A1D
	K(6,12) = K(12,6)
	K(12,7) = 2.99670606893001d-17*sqrt(temperature)	! 2A2D + 2A1D
	K(7,12) = K(12,7)
	K(12,8) = 3.00148303388612d-17*sqrt(temperature)	! 2A2D + 3A1D
	K(8,12) = K(12,8)
	K(12,9) = 3.05325346378201d-17*sqrt(temperature)	! 2A2D + 4A1D
	K(9,12) = K(12,9)
	K(12,10) = 3.81190872573705d-17*sqrt(temperature)	! 2A2D + 2D
	K(10,12) = K(12,10)
	K(12,11) = 3.29046951077028d-17*sqrt(temperature)	! 2A2D + 1A2D
	K(11,12) = K(12,11)
	K(12,12) = 3.18298598355307d-17*sqrt(temperature)	! 2A2D + 2A2D
	K(13,1) = 3.02705069238893d-17*sqrt(temperature)	! 3A2D + 1A
	K(1,13) = K(13,1)
	K(13,2) = 2.81058917527580d-17*sqrt(temperature)	! 3A2D + 2A
	K(2,13) = K(13,2)
	K(13,3) = 2.78951206647012d-17*sqrt(temperature)	! 3A2D + 3A
	K(3,13) = K(13,3)
	K(13,4) = 2.82389562385031d-17*sqrt(temperature)	! 3A2D + 4A
	K(4,13) = K(13,4)
	K(13,5) = 4.43959568621110d-17*sqrt(temperature)	! 3A2D + 1D
	K(5,13) = K(13,5)
	K(13,6) = 3.23641741717275d-17*sqrt(temperature)	! 3A2D + 1A1D
	K(6,13) = K(13,6)
	K(13,7) = 3.02671039469591d-17*sqrt(temperature)	! 3A2D + 2A1D
	K(7,13) = K(13,7)
	K(13,8) = 2.98505526994026d-17*sqrt(temperature)	! 3A2D + 3A1D
	K(8,13) = K(13,8)
	K(13,9) = 3.00139394642724d-17*sqrt(temperature)	! 3A2D + 4A1D
	K(9,13) = K(13,9)
	K(13,10) = 3.98120395504104d-17*sqrt(temperature)	! 3A2D + 2D
	K(10,13) = K(13,10)
	K(13,11) = 3.34974069815894d-17*sqrt(temperature)	! 3A2D + 1A2D
	K(11,13) = K(13,11)
	K(13,12) = 3.18341507188162d-17*sqrt(temperature)	! 3A2D + 2A2D
	K(12,13) = K(13,12)
	K(13,13) = 3.14250788514787d-17*sqrt(temperature)	! 3A2D + 3A2D
	K(14,2) = 2.89736207769365d-17*sqrt(temperature)	! 4A2D + 2A
	K(2,14) = K(14,2)
	K(14,3) = 2.83494506682650d-17*sqrt(temperature)	! 4A2D + 3A
	K(3,14) = K(14,3)
	K(14,4) = 2.83935545836062d-17*sqrt(temperature)	! 4A2D + 4A
	K(4,14) = K(14,4)
	K(14,5) = 4.71160370674057d-17*sqrt(temperature)	! 4A2D + 1D
	K(5,14) = K(14,5)
	K(14,6) = 3.35706769634303d-17*sqrt(temperature)	! 4A2D + 1A1D
	K(6,14) = K(14,6)
	K(14,7) = 3.09025016796110d-17*sqrt(temperature)	! 4A2D + 2A1D
	K(7,14) = K(14,7)
	K(14,8) = 3.01208396825086d-17*sqrt(temperature)	! 4A2D + 3A1D
	K(8,14) = K(14,8)
	K(14,9) = 3.00100047948565d-17*sqrt(temperature)	! 4A2D + 4A1D
	K(9,14) = K(14,9)
	K(14,10) = 4.16179821086092d-17*sqrt(temperature)	! 4A2D + 2D
	K(10,14) = K(14,10)
	K(14,11) = 3.43919941400074d-17*sqrt(temperature)	! 4A2D + 1A2D
	K(11,14) = K(14,11)
	K(14,12) = 3.22571933069290d-17*sqrt(temperature)	! 4A2D + 2A2D
	K(12,14) = K(14,12)
	K(14,13) = 3.15232262143403d-17*sqrt(temperature)	! 4A2D + 3A2D
	K(13,14) = K(14,13)
	K(14,14) = 3.13695384995609d-17*sqrt(temperature)	! 4A2D + 4A2D
	K(15,1) = 3.01989286415621d-17*sqrt(temperature)	! 3D + 1A
	K(1,15) = K(15,1)
	K(15,2) = 3.09059409105638d-17*sqrt(temperature)	! 3D + 2A
	K(2,15) = K(15,2)
	K(15,3) = 3.25786584976466d-17*sqrt(temperature)	! 3D + 3A
	K(3,15) = K(15,3)
	K(15,4) = 3.43793096764031d-17*sqrt(temperature)	! 3D + 4A
	K(4,15) = K(15,4)
	K(15,5) = 4.14494129362504d-17*sqrt(temperature)	! 3D + 1D
	K(5,15) = K(15,5)
	K(15,6) = 3.42995413502484d-17*sqrt(temperature)	! 3D + 1A1D
	K(6,15) = K(15,6)
	K(15,7) = 3.45257232044349d-17*sqrt(temperature)	! 3D + 2A1D
	K(7,15) = K(15,7)
	K(15,8) = 3.57517162136901d-17*sqrt(temperature)	! 3D + 3A1D
	K(8,15) = K(15,8)
	K(15,9) = 3.72264344722750d-17*sqrt(temperature)	! 3D + 4A1D
	K(9,15) = K(15,9)
	K(15,10) = 4.01147273043354d-17*sqrt(temperature)	! 3D + 2D
	K(10,15) = K(15,10)
	K(15,11) = 3.70781454324193d-17*sqrt(temperature)	! 3D + 1A2D
	K(11,15) = K(15,11)
	K(15,12) = 3.73720338109115d-17*sqrt(temperature)	! 3D + 2A2D
	K(12,15) = K(15,12)
	K(15,13) = 3.84243188126145d-17*sqrt(temperature)	! 3D + 3A2D
	K(13,15) = K(15,13)
	K(15,14) = 3.97193293405297d-17*sqrt(temperature)	! 3D + 4A2D
	K(14,15) = K(15,14)
	K(15,15) = 4.08850120510449d-17*sqrt(temperature)	! 3D + 3D
	K(16,1) = 3.02268581550821d-17*sqrt(temperature)	! 1A3D + 1A
	K(1,16) = K(16,1)
	K(16,2) = 2.92873245342293d-17*sqrt(temperature)	! 1A3D + 2A
	K(2,16) = K(16,2)
	K(16,3) = 2.98968525621429d-17*sqrt(temperature)	! 1A3D + 3A
	K(3,16) = K(16,3)
	K(16,4) = 3.08814187734470d-17*sqrt(temperature)	! 1A3D + 4A
	K(4,16) = K(16,4)
	K(16,5) = 4.31386042019245d-17*sqrt(temperature)	! 1A3D + 1D
	K(5,16) = K(16,5)
	K(16,6) = 3.31536198639159d-17*sqrt(temperature)	! 1A3D + 1A1D
	K(6,16) = K(16,6)
	K(16,7) = 3.20633359505197d-17*sqrt(temperature)	! 1A3D + 2A1D
	K(7,16) = K(16,7)
	K(16,8) = 3.23694035737798d-17*sqrt(temperature)	! 1A3D + 3A1D
	K(8,16) = K(16,8)
	K(16,9) = 3.31134646097083d-17*sqrt(temperature)	! 1A3D + 4A1D
	K(9,16) = K(16,9)
	K(16,10) = 3.98917132542620d-17*sqrt(temperature)	! 1A3D + 2D
	K(10,16) = K(16,10)
	K(16,11) = 3.49791470568189d-17*sqrt(temperature)	! 1A3D + 1A2D
	K(11,16) = K(16,11)
	K(16,12) = 3.41741442296998d-17*sqrt(temperature)	! 1A3D + 2A2D
	K(12,16) = K(16,12)
	K(16,13) = 3.44127863105541d-17*sqrt(temperature)	! 1A3D + 3A2D
	K(13,16) = K(16,13)
	K(16,14) = 3.50435906611843d-17*sqrt(temperature)	! 1A3D + 4A2D
	K(14,16) = K(16,14)
	K(16,15) = 3.94054540376438d-17*sqrt(temperature)	! 1A3D + 3D
	K(15,16) = K(16,15)
	K(16,16) = 3.64943444796786d-17*sqrt(temperature)	! 1A3D + 1A3D
	K(17,1) = 3.13194465823373d-17*sqrt(temperature)	! 2A3D + 1A
	K(1,17) = K(17,1)
	K(17,2) = 2.93359550362804d-17*sqrt(temperature)	! 2A3D + 2A
	K(2,17) = K(17,2)
	K(17,3) = 2.92991227193255d-17*sqrt(temperature)	! 2A3D + 3A
	K(3,17) = K(17,3)
	K(17,4) = 2.97984265951940d-17*sqrt(temperature)	! 2A3D + 4A
	K(4,17) = K(17,4)
	K(17,5) = 4.55620356051974d-17*sqrt(temperature)	! 2A3D + 1D
	K(5,17) = K(17,5)
	K(17,6) = 3.36010729643498d-17*sqrt(temperature)	! 2A3D + 1A1D
	K(6,17) = K(17,6)
	K(17,7) = 3.16706022115086d-17*sqrt(temperature)	! 2A3D + 2A1D
	K(7,17) = K(17,7)
	K(17,8) = 3.14101990762351d-17*sqrt(temperature)	! 2A3D + 3A1D
	K(8,17) = K(17,8)
	K(17,9) = 3.17148861253843d-17*sqrt(temperature)	! 2A3D + 4A1D
	K(9,17) = K(17,9)
	K(17,10) = 4.10654635590272d-17*sqrt(temperature)	! 2A3D + 2D
	K(10,17) = K(17,10)
	K(17,11) = 3.48944858931620d-17*sqrt(temperature)	! 2A3D + 1A2D
	K(11,17) = K(17,11)
	K(17,12) = 3.33884231970716d-17*sqrt(temperature)	! 2A3D + 2A2D
	K(12,17) = K(17,12)
	K(17,13) = 3.31235468974246d-17*sqrt(temperature)	! 2A3D + 3A2D
	K(13,17) = K(17,13)
	K(17,14) = 3.33527112669434d-17*sqrt(temperature)	! 2A3D + 4A2D
	K(14,17) = K(17,14)
	K(17,15) = 3.98119960434104d-17*sqrt(temperature)	! 2A3D + 3D
	K(15,17) = K(17,15)
	K(17,16) = 3.59554996296975d-17*sqrt(temperature)	! 2A3D + 1A3D
	K(16,17) = K(17,16)
	K(17,17) = 3.48145043488322d-17*sqrt(temperature)	! 2A3D + 2A3D
	K(18,1) = 3.26910501730354d-17*sqrt(temperature)	! 3A3D + 1A
	K(1,18) = K(18,1)
	K(18,2) = 2.99114746938185d-17*sqrt(temperature)	! 3A3D + 2A
	K(2,18) = K(18,2)
	K(18,3) = 2.94000521244333d-17*sqrt(temperature)	! 3A3D + 3A
	K(3,18) = K(18,3)
	K(18,4) = 2.95496442637427d-17*sqrt(temperature)	! 3A3D + 4A
	K(4,18) = K(18,4)
	K(18,5) = 4.80910847211470d-17*sqrt(temperature)	! 3A3D + 1D
	K(5,18) = K(18,5)
	K(18,6) = 3.45253049402141d-17*sqrt(temperature)	! 3A3D + 1A1D
	K(6,18) = K(18,6)
	K(18,7) = 3.19568882993891d-17*sqrt(temperature)	! 3A3D + 2A1D
	K(7,18) = K(18,7)
	K(18,8) = 3.12784270392478d-17*sqrt(temperature)	! 3A3D + 3A1D
	K(8,18) = K(18,8)
	K(18,9) = 3.12644654351575d-17*sqrt(temperature)	! 3A3D + 4A1D
	K(9,18) = K(18,9)
	K(18,10) = 4.26094098340396d-17*sqrt(temperature)	! 3A3D + 2D
	K(10,18) = K(18,10)
	K(18,11) = 3.54492926789563d-17*sqrt(temperature)	! 3A3D + 1A2D
	K(11,18) = K(18,11)
	K(18,12) = 3.34135780138988d-17*sqrt(temperature)	! 3A3D + 2A2D
	K(12,18) = K(18,12)
	K(18,13) = 3.27766124436665d-17*sqrt(temperature)	! 3A3D + 3A2D
	K(13,18) = K(18,13)
	K(18,14) = 3.27136209916630d-17*sqrt(temperature)	! 3A3D + 4A2D
	K(14,18) = K(18,14)
	K(18,15) = 4.07841143973692d-17*sqrt(temperature)	! 3A3D + 3D
	K(15,18) = K(18,15)
	K(18,16) = 3.61965065772395d-17*sqrt(temperature)	! 3A3D + 1A3D
	K(16,18) = K(18,16)
	K(18,17) = 3.46023050442563d-17*sqrt(temperature)	! 3A3D + 2A3D
	K(17,18) = K(18,17)
	K(18,18) = 3.40551869990935d-17*sqrt(temperature)	! 3A3D + 3A3D
	K(19,2) = 3.06994479028858d-17*sqrt(temperature)	! 4A3D + 2A
	K(2,19) = K(19,2)
	K(19,3) = 2.98057532107381d-17*sqrt(temperature)	! 4A3D + 3A
	K(3,19) = K(19,3)
	K(19,4) = 2.96783145190470d-17*sqrt(temperature)	! 4A3D + 4A
	K(4,19) = K(19,4)
	K(19,5) = 5.05871386086251d-17*sqrt(temperature)	! 4A3D + 1D
	K(5,19) = K(19,5)
	K(19,6) = 3.56268165927125d-17*sqrt(temperature)	! 4A3D + 1A1D
	K(6,19) = K(19,6)
	K(19,7) = 3.25315618247093d-17*sqrt(temperature)	! 4A3D + 2A1D
	K(7,19) = K(19,7)
	K(19,8) = 3.15159917438775d-17*sqrt(temperature)	! 4A3D + 3A1D
	K(8,19) = K(19,8)
	K(19,9) = 3.12487693463346d-17*sqrt(temperature)	! 4A3D + 4A1D
	K(9,19) = K(19,9)
	K(19,10) = 4.42622848625614d-17*sqrt(temperature)	! 4A3D + 2D
	K(10,19) = K(19,10)
	K(19,11) = 3.62642762902657d-17*sqrt(temperature)	! 4A3D + 1A2D
	K(11,19) = K(19,11)
	K(19,12) = 3.37939933755832d-17*sqrt(temperature)	! 4A3D + 2A2D
	K(12,19) = K(19,12)
	K(19,13) = 3.28572029416204d-17*sqrt(temperature)	! 4A3D + 3A2D
	K(13,19) = K(19,13)
	K(19,14) = 3.25617038244942d-17*sqrt(temperature)	! 4A3D + 4A2D
	K(14,19) = K(19,14)
	K(19,15) = 4.19680660673797d-17*sqrt(temperature)	! 4A3D + 3D
	K(15,19) = K(19,15)
	K(19,16) = 3.67700067222789d-17*sqrt(temperature)	! 4A3D + 1A3D
	K(16,19) = K(19,16)
	K(19,17) = 3.48057086067077d-17*sqrt(temperature)	! 4A3D + 2A3D
	K(17,19) = K(19,17)
	K(19,18) = 3.39890106915372d-17*sqrt(temperature)	! 4A3D + 3A3D
	K(18,19) = K(19,18)
	K(19,19) = 3.37082033742580d-17*sqrt(temperature)	! 4A3D + 4A3D
	K(20,1) = 3.21685600799811d-17*sqrt(temperature)	! 4D + 1A
	K(1,20) = K(20,1)
	K(20,2) = 3.17813958076630d-17*sqrt(temperature)	! 4D + 2A
	K(2,20) = K(20,2)
	K(20,3) = 3.28288064560400d-17*sqrt(temperature)	! 4D + 3A
	K(3,20) = K(20,3)
	K(20,4) = 3.41787613759507d-17*sqrt(temperature)	! 4D + 4A
	K(4,20) = K(20,4)
	K(20,6) = 3.56079873052069d-17*sqrt(temperature)	! 4D + 1A1D
	K(6,20) = K(20,6)
	K(20,7) = 3.49802822146243d-17*sqrt(temperature)	! 4D + 2A1D
	K(7,20) = K(20,7)
	K(20,8) = 3.56662582350804d-17*sqrt(temperature)	! 4D + 3A1D
	K(8,20) = K(20,8)
	K(20,9) = 3.67359591994433d-17*sqrt(temperature)	! 4D + 4A1D
	K(9,20) = K(20,9)
	K(20,10) = 4.22544800557049d-17*sqrt(temperature)	! 4D + 2D
	K(10,20) = K(20,10)
	K(20,11) = 3.78512964450544d-17*sqrt(temperature)	! 4D + 1A2D
	K(11,20) = K(20,11)
	K(20,12) = 3.74528908643381d-17*sqrt(temperature)	! 4D + 2A2D
	K(12,20) = K(20,12)
	K(20,13) = 3.80320535438534d-17*sqrt(temperature)	! 4D + 3A2D
	K(13,20) = K(20,13)
	K(20,14) = 3.89596266168685d-17*sqrt(temperature)	! 4D + 4A2D
	K(14,20) = K(20,14)
	K(20,15) = 4.21902818192331d-17*sqrt(temperature)	! 4D + 3D
	K(15,20) = K(20,15)
	K(20,16) = 3.97312985600776d-17*sqrt(temperature)	! 4D + 1A3D
	K(16,20) = K(20,16)
	K(20,17) = 3.95570581413427d-17*sqrt(temperature)	! 4D + 2A3D
	K(17,20) = K(20,17)
	K(20,18) = 4.01085960065454d-17*sqrt(temperature)	! 4D + 3A3D
	K(18,20) = K(20,18)
	K(20,19) = 4.09562384115962d-17*sqrt(temperature)	! 4D + 4A3D
	K(19,20) = K(20,19)
	K(20,20) = 4.28930820109269d-17*sqrt(temperature)	! 4D + 4D
	K(21,1) = 3.25552295528681d-17*sqrt(temperature)	! 1A4D + 1A
	K(1,21) = K(21,1)
	K(21,2) = 3.08442585795698d-17*sqrt(temperature)	! 1A4D + 2A
	K(2,21) = K(21,2)
	K(21,3) = 3.10461905317041d-17*sqrt(temperature)	! 1A4D + 3A
	K(3,21) = K(21,3)
	K(21,4) = 3.17521980187661d-17*sqrt(temperature)	! 1A4D + 4A
	K(4,21) = K(21,4)
	K(21,6) = 3.50990221881463d-17*sqrt(temperature)	! 1A4D + 1A1D
	K(6,21) = K(21,6)
	K(21,7) = 3.34093584261907d-17*sqrt(temperature)	! 1A4D + 2A1D
	K(7,21) = K(21,7)
	K(21,8) = 3.33606697988163d-17*sqrt(temperature)	! 1A4D + 3A1D
	K(8,21) = K(21,8)
	K(21,9) = 3.38519461349218d-17*sqrt(temperature)	! 1A4D + 4A1D
	K(9,21) = K(21,9)
	K(21,10) = 4.25463543997013d-17*sqrt(temperature)	! 1A4D + 2D
	K(10,21) = K(21,10)
	K(21,11) = 3.66124816612156d-17*sqrt(temperature)	! 1A4D + 1A2D
	K(11,21) = K(21,11)
	K(21,12) = 3.53266365635807d-17*sqrt(temperature)	! 1A4D + 2A2D
	K(12,21) = K(21,12)
	K(21,13) = 3.52549556674088d-17*sqrt(temperature)	! 1A4D + 3A2D
	K(13,21) = K(21,13)
	K(21,14) = 3.56559835847258d-17*sqrt(temperature)	! 1A4D + 4A2D
	K(14,21) = K(21,14)
	K(21,15) = 4.14956538736041d-17*sqrt(temperature)	! 1A4D + 3D
	K(15,21) = K(21,15)
	K(21,16) = 3.78704123427971d-17*sqrt(temperature)	! 1A4D + 1A3D
	K(16,21) = K(21,16)
	K(21,17) = 3.69325295344486d-17*sqrt(temperature)	! 1A4D + 2A3D
	K(17,21) = K(21,17)
	K(21,18) = 3.68990265697328d-17*sqrt(temperature)	! 1A4D + 3A3D
	K(18,21) = K(21,18)
	K(21,19) = 3.72625724168467d-17*sqrt(temperature)	! 1A4D + 4A3D
	K(19,21) = K(21,19)
	K(21,20) = 4.14351238216624d-17*sqrt(temperature)	! 1A4D + 4D
	K(20,21) = K(21,20)
	K(21,21) = 3.90275657598138d-17*sqrt(temperature)	! 1A4D + 1A4D
	K(22,1) = 3.36215217793290d-17*sqrt(temperature)	! 2A4D + 1A
	K(1,22) = K(22,1)
	K(22,2) = 3.09922908876576d-17*sqrt(temperature)	! 2A4D + 2A
	K(2,22) = K(22,2)
	K(22,3) = 3.06295847145279d-17*sqrt(temperature)	! 2A4D + 3A
	K(3,22) = K(22,3)
	K(22,4) = 3.09133313070854d-17*sqrt(temperature)	! 2A4D + 4A
	K(4,22) = K(22,4)
	K(22,6) = 3.56118881601145d-17*sqrt(temperature)	! 2A4D + 1A1D
	K(6,22) = K(22,6)
	K(22,7) = 3.31844631680973d-17*sqrt(temperature)	! 2A4D + 2A1D
	K(7,22) = K(22,7)
	K(22,8) = 3.26404520309648d-17*sqrt(temperature)	! 2A4D + 3A1D
	K(8,22) = K(22,8)
	K(22,9) = 3.27489296842112d-17*sqrt(temperature)	! 2A4D + 4A1D
	K(9,22) = K(22,9)
	K(22,10) = 4.37129956150700d-17*sqrt(temperature)	! 2A4D + 2D
	K(10,22) = K(22,10)
	K(22,11) = 3.66703556479633d-17*sqrt(temperature)	! 2A4D + 1A2D
	K(11,22) = K(22,11)
	K(22,12) = 3.47692394054407d-17*sqrt(temperature)	! 2A4D + 2A2D
	K(12,22) = K(22,12)
	K(22,13) = 3.42570835541366d-17*sqrt(temperature)	! 2A4D + 3A2D
	K(13,22) = K(22,13)
	K(22,14) = 3.43079869501813d-17*sqrt(temperature)	! 2A4D + 4A2D
	K(14,22) = K(22,14)
	K(22,15) = 4.19974418904366d-17*sqrt(temperature)	! 2A4D + 3D
	K(15,22) = K(22,15)
	K(22,16) = 3.75407469949633d-17*sqrt(temperature)	! 2A4D + 1A3D
	K(16,22) = K(22,16)
	K(22,17) = 3.60744533264302d-17*sqrt(temperature)	! 2A4D + 2A3D
	K(17,22) = K(22,17)
	K(22,18) = 3.56444154628882d-17*sqrt(temperature)	! 2A4D + 3A3D
	K(18,22) = K(22,18)
	K(22,19) = 3.56852699615051d-17*sqrt(temperature)	! 2A4D + 4A3D
	K(19,22) = K(22,19)
	K(22,20) = 4.14368002990175d-17*sqrt(temperature)	! 2A4D + 4D
	K(20,22) = K(22,20)
	K(22,21) = 3.83579100229179d-17*sqrt(temperature)	! 2A4D + 1A4D
	K(21,22) = K(22,21)
	K(22,22) = 3.72245634097107d-17*sqrt(temperature)	! 2A4D + 2A4D
	K(23,1) = 3.49043670857906d-17*sqrt(temperature)	! 3A4D + 1A
	K(1,23) = K(23,1)
	K(23,2) = 3.15476588186036d-17*sqrt(temperature)	! 3A4D + 2A
	K(2,23) = K(23,2)
	K(23,3) = 3.07531638996627d-17*sqrt(temperature)	! 3A4D + 3A
	K(3,23) = K(23,3)
	K(23,4) = 3.07194206215132d-17*sqrt(temperature)	! 3A4D + 4A
	K(4,23) = K(23,4)
	K(23,6) = 3.64901348708458d-17*sqrt(temperature)	! 3A4D + 1A1D
	K(6,23) = K(23,6)
	K(23,7) = 3.34815796368918d-17*sqrt(temperature)	! 3A4D + 2A1D
	K(7,23) = K(23,7)
	K(23,8) = 3.25574226703014d-17*sqrt(temperature)	! 3A4D + 3A1D
	K(8,23) = K(23,8)
	K(23,9) = 3.23765740777518d-17*sqrt(temperature)	! 3A4D + 4A1D
	K(9,23) = K(23,9)
	K(23,10) = 4.51603387880226d-17*sqrt(temperature)	! 3A4D + 2D
	K(10,23) = K(23,10)
	K(23,11) = 3.72165162582506d-17*sqrt(temperature)	! 3A4D + 1A2D
	K(11,23) = K(23,11)
	K(23,12) = 3.48334426724349d-17*sqrt(temperature)	! 3A4D + 2A2D
	K(12,23) = K(23,12)
	K(23,13) = 3.39830001973208d-17*sqrt(temperature)	! 3A4D + 3A2D
	K(13,23) = K(23,13)
	K(23,14) = 3.37686092929410d-17*sqrt(temperature)	! 3A4D + 4A2D
	K(14,23) = K(23,14)
	K(23,15) = 4.29274749301986d-17*sqrt(temperature)	! 3A4D + 3D
	K(15,23) = K(23,15)
	K(23,16) = 3.78056795829335d-17*sqrt(temperature)	! 3A4D + 1A3D
	K(16,23) = K(23,16)
	K(23,17) = 3.59271247979772d-17*sqrt(temperature)	! 3A4D + 2A3D
	K(17,23) = K(23,17)
	K(23,18) = 3.51926015491360d-17*sqrt(temperature)	! 3A4D + 3A3D
	K(18,23) = K(23,18)
	K(23,19) = 3.49886518498441d-17*sqrt(temperature)	! 3A4D + 4A3D
	K(19,23) = K(23,19)
	K(23,20) = 4.19881855303853d-17*sqrt(temperature)	! 3A4D + 4D
	K(20,23) = K(23,20)
	K(23,21) = 3.83772753878316d-17*sqrt(temperature)	! 3A4D + 1A4D
	K(21,23) = K(23,21)
	K(23,22) = 3.68832257746472d-17*sqrt(temperature)	! 3A4D + 2A4D
	K(22,23) = K(23,22)
	K(23,23) = 3.62650661050272d-17*sqrt(temperature)	! 3A4D + 3A4D
	K(24,1) = 3.62596102037558d-17*sqrt(temperature)	! 4A4D + 1A
	K(1,24) = K(24,1)
	K(24,2) = 3.22850887225448d-17*sqrt(temperature)	! 4A4D + 2A
	K(2,24) = K(24,2)
	K(24,3) = 3.11350996316107d-17*sqrt(temperature)	! 4A4D + 3A
	K(3,24) = K(24,3)
	K(24,4) = 3.08440375491246d-17*sqrt(temperature)	! 4A4D + 4A
	K(4,24) = K(24,4)
	K(24,6) = 3.75206694350681d-17*sqrt(temperature)	! 4A4D + 1A1D
	K(6,24) = K(24,6)
	K(24,7) = 3.40225187661203d-17*sqrt(temperature)	! 4A4D + 2A1D
	K(7,24) = K(24,7)
	K(24,8) = 3.27852539241479d-17*sqrt(temperature)	! 4A4D + 3A1D
	K(8,24) = K(24,8)
	K(24,9) = 3.23692740187744d-17*sqrt(temperature)	! 4A4D + 4A1D
	K(9,24) = K(24,9)
	K(24,10) = 4.67037485065157d-17*sqrt(temperature)	! 4A4D + 2D
	K(10,24) = K(24,10)
	K(24,11) = 3.79823162049186d-17*sqrt(temperature)	! 4A4D + 1A2D
	K(11,24) = K(24,11)
	K(24,12) = 3.51958859165781d-17*sqrt(temperature)	! 4A4D + 2A2D
	K(12,24) = K(24,12)
	K(24,13) = 3.40671002718251d-17*sqrt(temperature)	! 4A4D + 3A2D
	K(13,24) = K(24,13)
	K(24,14) = 3.36368883231231d-17*sqrt(temperature)	! 4A4D + 4A2D
	K(14,24) = K(24,14)
	K(24,15) = 4.40367706124186d-17*sqrt(temperature)	! 4A4D + 3D
	K(15,24) = K(24,15)
	K(24,16) = 3.83490880292159d-17*sqrt(temperature)	! 4A4D + 1A3D
	K(16,24) = K(24,16)
	K(24,17) = 3.61271601610681d-17*sqrt(temperature)	! 4A4D + 2A3D
	K(17,24) = K(24,17)
	K(24,18) = 3.51424119700538d-17*sqrt(temperature)	! 4A4D + 3A3D
	K(18,24) = K(24,18)
	K(24,19) = 3.47392500725570d-17*sqrt(temperature)	! 4A4D + 4A3D
	K(19,24) = K(24,19)
	K(24,20) = 4.27872287625918d-17*sqrt(temperature)	! 4A4D + 4D
	K(20,24) = K(24,20)
	K(24,21) = 3.87276959947316d-17*sqrt(temperature)	! 4A4D + 1A4D
	K(21,24) = K(24,21)
	K(24,22) = 3.69342714695874d-17*sqrt(temperature)	! 4A4D + 2A4D
	K(22,24) = K(24,22)
	K(24,23) = 3.60888862121434d-17*sqrt(temperature)	! 4A4D + 3A4D
	K(23,24) = K(24,23)
	K(24,24) = 3.57278096683900d-17*sqrt(temperature)	! 4A4D + 4A4D

end subroutine get_collD

!-----------------------------------------------------------

subroutine get_evapD(E,K,temperature)
	implicit none
	real(kind(1.d0)) :: E(24,24), K(24,24), temperature

	! evaporation coefficients

	E = 0.d0
	E(1,1) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-17.060d0/temperature-(-34.895d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(1,1)	! 2A -> 1A + 1A
	E(2,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-33.728d0/temperature-(-75.053d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(2,1)	! 3A -> 2A + 1A
	E(1,2) = E(2,1)
	E(3,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-50.317d0/temperature-(-104.448d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(3,1)	! 4A -> 3A + 1A
	E(1,3) = E(3,1)
	E(2,2) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-50.317d0/temperature-(-104.448d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(2,2)	! 4A -> 2A + 2A
	E(5,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(5,1)	! 1A1D -> 1D + 1A
	E(1,5) = E(5,1)
	E(6,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(6,1)	! 2A1D -> 1A1D + 1A
	E(1,6) = E(6,1)
	E(5,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-0.d0&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,2)	! 2A1D -> 1D + 2A
	E(2,5) = E(5,2)
	E(7,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(7,1)	! 3A1D -> 2A1D + 1A
	E(1,7) = E(7,1)
	E(6,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,2)	! 3A1D -> 1A1D + 2A
	E(2,6) = E(6,2)
	E(5,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-0.d0&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,3)	! 3A1D -> 1D + 3A
	E(3,5) = E(5,3)
	E(8,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-93.864d0/temperature-(-155.847d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(8,1)	! 4A1D -> 3A1D + 1A
	E(1,8) = E(8,1)
	E(7,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-93.864d0/temperature-(-155.847d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,2)	! 4A1D -> 2A1D + 2A
	E(2,7) = E(7,2)
	E(6,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-93.864d0/temperature-(-155.847d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,3)	! 4A1D -> 1A1D + 3A
	E(3,6) = E(6,3)
	E(5,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-93.864d0/temperature-(-155.847d0)/1.d3)&
			 &-0.d0&
			 &-(-50.317d0/temperature-(-104.448d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,4)	! 4A1D -> 1D + 4A
	E(4,5) = E(5,4)
	E(5,5) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(5,5)	! 2D -> 1D + 1D
	E(10,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(10,1)	! 1A2D -> 2D + 1A
	E(1,10) = E(10,1)
	E(6,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(6,5)	! 1A2D -> 1A1D + 1D
	E(5,6) = E(6,5)
	E(11,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(11,1)	! 2A2D -> 1A2D + 1A
	E(1,11) = E(11,1)
	E(10,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,2)	! 2A2D -> 2D + 2A
	E(2,10) = E(10,2)
	E(7,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(7,5)	! 2A2D -> 2A1D + 1D
	E(5,7) = E(7,5)
	E(6,6) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,6)	! 2A2D -> 1A1D + 1A1D
	E(12,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(12,1)	! 3A2D -> 2A2D + 1A
	E(1,12) = E(12,1)
	E(11,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,2)	! 3A2D -> 1A2D + 2A
	E(2,11) = E(11,2)
	E(10,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,3)	! 3A2D -> 2D + 3A
	E(3,10) = E(10,3)
	E(8,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(8,5)	! 3A2D -> 3A1D + 1D
	E(5,8) = E(8,5)
	E(7,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,6)	! 3A2D -> 2A1D + 1A1D
	E(6,7) = E(7,6)
	E(13,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(13,1)	! 4A2D -> 3A2D + 1A
	E(1,13) = E(13,1)
	E(12,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,2)	! 4A2D -> 2A2D + 2A
	E(2,12) = E(12,2)
	E(11,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,3)	! 4A2D -> 1A2D + 3A
	E(3,11) = E(11,3)
	E(10,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-50.317d0/temperature-(-104.448d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,4)	! 4A2D -> 2D + 4A
	E(4,10) = E(10,4)
	E(9,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-93.864d0/temperature-(-155.847d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(9,5)	! 4A2D -> 4A1D + 1D
	E(5,9) = E(9,5)
	E(8,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(8,6)	! 4A2D -> 3A1D + 1A1D
	E(6,8) = E(8,6)
	E(7,7) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,7)	! 4A2D -> 2A1D + 2A1D
	E(10,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(10,5)	! 3D -> 2D + 1D
	E(5,10) = E(10,5)
	E(15,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(15,1)	! 1A3D -> 3D + 1A
	E(1,15) = E(15,1)
	E(11,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(11,5)	! 1A3D -> 1A2D + 1D
	E(5,11) = E(11,5)
	E(10,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,6)	! 1A3D -> 2D + 1A1D
	E(6,10) = E(10,6)
	E(16,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(16,1)	! 2A3D -> 1A3D + 1A
	E(1,16) = E(16,1)
	E(15,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,2)	! 2A3D -> 3D + 2A
	E(2,15) = E(15,2)
	E(12,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(12,5)	! 2A3D -> 2A2D + 1D
	E(5,12) = E(12,5)
	E(11,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,6)	! 2A3D -> 1A2D + 1A1D
	E(6,11) = E(11,6)
	E(10,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,7)	! 2A3D -> 2D + 2A1D
	E(7,10) = E(10,7)
	E(17,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(17,1)	! 3A3D -> 2A3D + 1A
	E(1,17) = E(17,1)
	E(16,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,2)	! 3A3D -> 1A3D + 2A
	E(2,16) = E(16,2)
	E(15,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,3)	! 3A3D -> 3D + 3A
	E(3,15) = E(15,3)
	E(13,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(13,5)	! 3A3D -> 3A2D + 1D
	E(5,13) = E(13,5)
	E(12,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,6)	! 3A3D -> 2A2D + 1A1D
	E(6,12) = E(12,6)
	E(11,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,7)	! 3A3D -> 1A2D + 2A1D
	E(7,11) = E(11,7)
	E(10,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,8)	! 3A3D -> 2D + 3A1D
	E(8,10) = E(10,8)
	E(18,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(18,1)	! 4A3D -> 3A3D + 1A
	E(1,18) = E(18,1)
	E(17,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,2)	! 4A3D -> 2A3D + 2A
	E(2,17) = E(17,2)
	E(16,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,3)	! 4A3D -> 1A3D + 3A
	E(3,16) = E(16,3)
	E(15,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-50.317d0/temperature-(-104.448d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,4)	! 4A3D -> 3D + 4A
	E(4,15) = E(15,4)
	E(14,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(14,5)	! 4A3D -> 4A2D + 1D
	E(5,14) = E(14,5)
	E(13,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,6)	! 4A3D -> 3A2D + 1A1D
	E(6,13) = E(13,6)
	E(12,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,7)	! 4A3D -> 2A2D + 2A1D
	E(7,12) = E(12,7)
	E(11,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,8)	! 4A3D -> 1A2D + 3A1D
	E(8,11) = E(11,8)
	E(10,9) = 7.33893243358348d+27/temperature*exp(&
			 &((-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-93.864d0/temperature-(-155.847d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,9)	! 4A3D -> 2D + 4A1D
	E(9,10) = E(10,9)
	E(15,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(15,5)	! 4D -> 3D + 1D
	E(5,15) = E(15,5)
	E(10,10) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,10)	! 4D -> 2D + 2D
	E(20,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(20,1)	! 1A4D -> 4D + 1A
	E(1,20) = E(20,1)
	E(16,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(16,5)	! 1A4D -> 1A3D + 1D
	E(5,16) = E(16,5)
	E(15,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,6)	! 1A4D -> 3D + 1A1D
	E(6,15) = E(15,6)
	E(11,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,10)	! 1A4D -> 1A2D + 2D
	E(10,11) = E(11,10)
	E(21,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(21,1)	! 2A4D -> 1A4D + 1A
	E(1,21) = E(21,1)
	E(20,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,2)	! 2A4D -> 4D + 2A
	E(2,20) = E(20,2)
	E(17,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(17,5)	! 2A4D -> 2A3D + 1D
	E(5,17) = E(17,5)
	E(16,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,6)	! 2A4D -> 1A3D + 1A1D
	E(6,16) = E(16,6)
	E(15,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,7)	! 2A4D -> 3D + 2A1D
	E(7,15) = E(15,7)
	E(12,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,10)	! 2A4D -> 2A2D + 2D
	E(10,12) = E(12,10)
	E(11,11) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,11)	! 2A4D -> 1A2D + 1A2D
	E(22,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(22,1)	! 3A4D -> 2A4D + 1A
	E(1,22) = E(22,1)
	E(21,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(21,2)	! 3A4D -> 1A4D + 2A
	E(2,21) = E(21,2)
	E(20,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,3)	! 3A4D -> 4D + 3A
	E(3,20) = E(20,3)
	E(18,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(18,5)	! 3A4D -> 3A3D + 1D
	E(5,18) = E(18,5)
	E(17,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,6)	! 3A4D -> 2A3D + 1A1D
	E(6,17) = E(17,6)
	E(16,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,7)	! 3A4D -> 1A3D + 2A1D
	E(7,16) = E(16,7)
	E(15,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,8)	! 3A4D -> 3D + 3A1D
	E(8,15) = E(15,8)
	E(13,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,10)	! 3A4D -> 3A2D + 2D
	E(10,13) = E(13,10)
	E(12,11) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,11)	! 3A4D -> 2A2D + 1A2D
	E(11,12) = E(12,11)
	E(23,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-158.404d0/temperature-(-241.213d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(23,1)	! 4A4D -> 3A4D + 1A
	E(1,23) = E(23,1)
	E(22,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-112.768d0/temperature-(-188.962d0)/1.d3)&
			 &-(-17.060d0/temperature-(-34.895d0)/1.d3))&
			 &*5.03218937158374d+02)*K(22,2)	! 4A4D -> 2A4D + 2A
	E(2,22) = E(22,2)
	E(21,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-62.415d0/temperature-(-150.593d0)/1.d3)&
			 &-(-33.728d0/temperature-(-75.053d0)/1.d3))&
			 &*5.03218937158374d+02)*K(21,3)	! 4A4D -> 1A4D + 3A
	E(3,21) = E(21,3)
	E(20,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-17.935d0/temperature-(-93.100d0)/1.d3)&
			 &-(-50.317d0/temperature-(-104.448d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,4)	! 4A4D -> 4D + 4A
	E(4,20) = E(20,4)
	E(19,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-161.029d0/temperature-(-240.577d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(19,5)	! 4A4D -> 4A3D + 1D
	E(5,19) = E(19,5)
	E(18,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-136.338d0/temperature-(-195.832d0)/1.d3)&
			 &-(-21.603d0/temperature-(-28.937d0)/1.d3))&
			 &*5.03218937158374d+02)*K(18,6)	! 4A4D -> 3A3D + 1A1D
	E(6,18) = E(18,6)
	E(17,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-98.720d0/temperature-(-156.132d0)/1.d3)&
			 &-(-52.957d0/temperature-(-73.226d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,7)	! 4A4D -> 2A3D + 2A1D
	E(7,17) = E(17,7)
	E(16,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-50.339d0/temperature-(-106.156d0)/1.d3)&
			 &-(-76.335d0/temperature-(-116.810d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,8)	! 4A4D -> 1A3D + 3A1D
	E(8,16) = E(16,8)
	E(15,9) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-11.864d0/temperature-(-71.068d0)/1.d3)&
			 &-(-93.864d0/temperature-(-155.847d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,9)	! 4A4D -> 3D + 4A1D
	E(9,15) = E(15,9)
	E(14,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-130.544d0/temperature-(-203.306d0)/1.d3)&
			 &-(-3.278d0/temperature-(-27.329d0)/1.d3))&
			 &*5.03218937158374d+02)*K(14,10)	! 4A4D -> 4A2D + 2D
	E(10,14) = E(14,10)
	E(13,11) = 7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-107.699d0/temperature-(-154.733d0)/1.d3)&
			 &-(-35.334d0/temperature-(-68.452d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,11)	! 4A4D -> 3A2D + 1A2D
	E(11,13) = E(13,11)
	E(12,12) = 0.5*7.33893243358348d+27/temperature*exp(&
			 &((-189.302d0/temperature-(-270.547d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3)&
			 &-(-78.988d0/temperature-(-114.033d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,12)	! 4A4D -> 2A2D + 2A2D

end subroutine get_evapD
