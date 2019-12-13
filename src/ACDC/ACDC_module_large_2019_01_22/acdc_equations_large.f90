!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!
!
! Cluster set will be restricted according to the criteria given by the user
! Available clusters from which the set is selected:
! 21041 clusters with
!   0-6 A molecules
!   0-5 N molecules
!   0-500 O molecules
! and the last element is the outgoing flux
!
!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!

! differential equations: f = dc/dt
subroutine feval(neqn,t,c,f,coef,ipar)
use acdc_system, only : monomer_indices
use monomer_settings, only : sources_and_constants
	implicit none
	integer, parameter :: nclust = 21041
	integer :: neqn, ipar(4), i, j, ij, n, ij_ind(3), imol
	real(kind(1.d0)) :: t, c(neqn), f(neqn), coef(2)
	logical, save, allocatable :: isconst(:)
	real(kind(1.d0)), save, allocatable :: K(:,:), E(:,:), loss(:), source(:)
	integer, save, allocatable :: fitted(:,:)
	integer, save :: nclust_inc=0, indices(nclust,3), clust_from_indices(0:12,0:10,0:1000), n_monomers(3)
	real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)

	! system size parameters
	if (nclust_inc .eq. 0) then
		call get_molecule_numbers(nclust_inc,indices)
		call get_cluster_numbers(clust_from_indices)
		call monomer_indices(n_monomers)
		allocate(isconst(nclust_inc))
		allocate(fitted(nclust_inc,0:nclust_inc))
		allocate(K(nclust_inc,nclust_inc))
		allocate(E(nclust_inc,nclust_inc))
		allocate(loss(nclust_inc))
		allocate(source(nclust_inc))
	end if

	! these parameters are read at the very beginning and eg. if source terms or collision rates have changed
	if (ipar(1) .eq. 0) then
		! after this the same values are used until some other routine tells otherwise
		ipar(1) = 1
		call sources_and_constants(nclust_inc,source,isconst,fitted)
		call get_coll(nclust_inc,K,coef(1))
		call get_evap(nclust_inc,E,K,coef(1))
		do i = 1,nclust_inc
			K(i,i) = 0.5d0*K(i,i)
		end do
		call get_losses(nclust_inc,loss,coef)
	end if

	f = 0.d0

	do i = 1,nclust_inc
		f(i) = f(i)+source(i)-loss(i)*c(i)
		do j = i,nclust_inc
			ij_ind = indices(i,:)+indices(j,:)
			if (clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3)) .eq. 0) then
				ij = nclust_inc+1
			else if (clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3)) .eq. -1) then
				cycle
			else
				ij = clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3))
			end if
			f(i) = f(i)-K(j,i)*c(i)*c(j)
			f(j) = f(j)-K(j,i)*c(i)*c(j)
			f(ij) = f(ij)+K(j,i)*c(i)*c(j)
			if (E(j,i) .gt. 0.d0) then
				f(i) = f(i)+E(j,i)*c(ij)
				f(j) = f(j)+E(j,i)*c(ij)
				f(ij) = f(ij)-E(j,i)*c(ij)
			end if
		end do
	end do
	do i = 1,nclust_inc
		if (isconst(i)) then
			f(i) = 0.d0
		end if
	end do
	do n = 1, fitted(1,0) ! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)
		f(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations
		do i = 1, fitted(n+1,1) ! loop over the clusters used for the fitting
			f(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i+1)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)
		end do
	end do

end subroutine feval

!-----------------------------------------------------------

! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)
! not using this, since the solver works faster using a numerical jacobian...
subroutine jeval(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)
	implicit none
	integer :: ldim,neqn,ierr,ipar(4),i,j,k,n
	real(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn),coef(2),ml,mu

end subroutine jeval

!-----------------------------------------------------------

subroutine formation(neqn,c,j_out,ipar,coef)
	implicit none
	integer, parameter :: nclust = 21041, ij_ind_max(3) = (/6,5,500/)
	integer :: neqn, i, j, ij_ind(3), imol, n, ipar(4)
	real(kind(1.d0)) :: c(neqn), j_out,coef(2)
	real(kind(1.d0)), save, allocatable :: K(:,:)
	integer, save :: nclust_inc=0, indices(nclust,3), clust_from_indices(0:12,0:10,0:1000)

	! system size parameters
	if (nclust_inc .eq. 0) then
		call get_molecule_numbers(nclust_inc,indices)
		call get_cluster_numbers(clust_from_indices)
		allocate(K(nclust_inc,nclust_inc))
	end if

	if (ipar(3) .eq. 0) then
		ipar(3) = 1
		call get_coll(nclust_inc,K,coef(1))
		do i = 1,nclust_inc
			K(i,i) = 0.5d0*K(i,i)
		end do
	end if

	j_out = 0.d0
	do i = 1,nclust_inc
		do j = i,nclust_inc
			ij_ind = indices(i,:)+indices(j,:)
			if (clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3)) .eq. 0) then
				j_out = j_out+K(j,i)*c(i)*c(j)
			end if
		end do
	end do

end subroutine formation

!-----------------------------------------------------------

subroutine get_masses_and_radii(m,r,i_ind)
	implicit none
	real(kind(1.d0)), parameter :: u2kg = 1.660538782d-27, kB=1.3806504d-23, pi=4.d0*atan(1.d0)
	real(kind(1.d0)), parameter :: m1 = 9.80800000000000d+01		! in g/mol
	real(kind(1.d0)), parameter :: m2 = 1.70400000000000d+01		! in g/mol
	real(kind(1.d0)), parameter :: m3 = 3.00000000000000d+02		! in g/mol
	real(kind(1.d0)), parameter :: monvol1 = 8.89976195915848d-29		! in m^3
	real(kind(1.d0)), parameter :: monvol2 = 4.06545702084667d-29		! in m^3
	real(kind(1.d0)), parameter :: monvol3 = 3.55829739249156d-28		! in m^3
	real(kind(1.d0)) :: r, m
	integer :: i_ind(3)

	m = (real(i_ind(1),kind=kind(1.d0))*m1+real(i_ind(2),kind=kind(1.d0))*m2+&
		&real(i_ind(3),kind=kind(1.d0))*m3)
	r = (3.d0/4.d0/pi*(real(i_ind(1),kind=kind(1.d0))*monvol1+real(i_ind(2),kind=kind(1.d0))*monvol2+&
		&real(i_ind(3),kind=kind(1.d0))*monvol3))**(1.d0/3.d0)


end subroutine get_masses_and_radii

!-----------------------------------------------------------

subroutine get_molecule_numbers(nclust_inc,indices)
	implicit none
	integer, parameter :: nclust = 21041
	integer :: n, i1, i2, i3, nclust_inc, indices(nclust,3), clust_from_indices(0:12,0:10,0:1000)

	nclust_inc = 0
	indices = 0

	call get_cluster_numbers(clust_from_indices)

	n = 0
	do i1 = 0,6
		do i2 = 0,5
			do i3 = 0,500
				if (clust_from_indices(i1, i2, i3) .gt. 0) then
					n = n+1
					indices(clust_from_indices(i1, i2, i3),:) = (/i1, i2, i3/)
				end if
			end do
		end do
	end do

	nclust_inc = n

end subroutine get_molecule_numbers

!-----------------------------------------------------------

subroutine get_cluster_numbers(clust_from_indices)
use restriction_criteria
	implicit none
	integer :: n, i1, i2, i3, ij_ind(3), clust_from_indices(0:12,0:10,0:1000), valid

	clust_from_indices = 0

	n = 0
	do i1 = 0,12
		do i2 = 0,10
			do i3 = 0,1000
				if (i1+i2+i3>0) then
					call check_validity((/i1, i2, i3/),valid)
					if (valid .eq. 1) then
						n = n+1
						clust_from_indices(i1, i2, i3) = n
					else
						clust_from_indices(i1, i2, i3) = valid
					end if
				end if
			end do
		end do
	end do

end subroutine get_cluster_numbers

!-----------------------------------------------------------

subroutine get_losses(nrates,loss,coef)
	implicit none
	real(kind(1.d0)) :: cs(nrates), loss(nrates), m, r, r_ref, coef(2)
	integer, parameter :: nclust = 21041
	integer :: nrates, i
	integer, save :: nclust_inc=0, indices(nclust,3), clust_from_indices(0:12,0:10,0:1000)

	if (nclust_inc .eq. 0) then
		call get_molecule_numbers(nclust_inc,indices)
		call get_cluster_numbers(clust_from_indices)
	end if

	! coagulation sink

	! variable cs: the following values only give the size dependence of the sink,
	! and will be multiplied by the size-independent factor given as input

	! reference radius corresponding to 1A
	call get_masses_and_radii(m,r_ref,(/1,0,0/))

	do i = 1,nclust_inc
		call get_masses_and_radii(m,r,indices(i,:))
		cs(i) = (r/r_ref)**(-1.60000000000000d+00)
	end do

	cs(clust_from_indices(1,0,0)) = 0.00000000000000d+00
	cs(clust_from_indices(0,1,0)) = 0.00000000000000d+00
	cs(clust_from_indices(0,0,1)) = 0.00000000000000d+00

	cs = coef(2)*cs

	loss = cs

end subroutine get_losses

!-----------------------------------------------------------

subroutine get_coll(nrates,K,temperature)
	implicit none
	integer, parameter :: nclust = 21041
	real(kind(1.d0)) :: K(nrates,nrates), temperature
	real(kind(1.d0)) :: ri, rj, mi, mj
	integer :: nrates, i, j
	integer, save :: nclust_inc=0, indices(nclust,3)

	if (nclust_inc .eq. 0) then
		call get_molecule_numbers(nclust_inc,indices)
	end if

	K = 0.d0

	do i = 1,nclust_inc
		call get_masses_and_radii(mi,ri,indices(i,:))
		do j = i,nclust_inc
			call get_masses_and_radii(mj,rj,indices(j,:))
			K(i,j) = sqrt(2.08965485072494d+05*temperature*(1.d0/mi+1.d0/mj))*(ri+rj)**2.d0
			K(j,i) = K(i,j)
		end do
	end do


end subroutine get_coll

!-----------------------------------------------------------

subroutine get_evap(nrates,E,K,temperature)
use acdc_system, only : monomer_indices
use acdc_hsdata, only : get_hsdata
	implicit none
	integer, parameter :: nclust = 21041
	real(kind(1.d0)) :: E(nrates,nrates), K(nrates,nrates), temperature
	integer, save :: nclust_inc=0, indices(nclust,3), clust_from_indices(0:12,0:10,0:1000)
	real(kind(1.d0)) :: m, r, xmol
	integer :: nrates, ij_ind(3), i, j, ij, imol, ind(2)
	logical, save, allocatable :: lhs(:)
	real(kind(1.d0)), save, allocatable :: deltah(:), deltas(:)
	integer, save :: n_monomers(3)
	real(kind(1.d0)), parameter :: psat(3) = (/0.00000000000000d+00,0.00000000000000d+00,1.00000000000000d-12/)
	real(kind(1.d0)), parameter :: monvol(3) = (/8.89976195915848d-29,4.06545702084667d-29,3.55829739249156d-28/)
	real(kind(1.d0)), parameter :: sigma = 5.00000000000000d-02

	if (nclust_inc .eq. 0) then
		call get_molecule_numbers(nclust_inc,indices)
		call get_cluster_numbers(clust_from_indices)
		call monomer_indices(n_monomers)
		allocate(lhs(nrates))
		allocate(deltah(nrates))
		allocate(deltas(nrates))
		do i = 1,nclust_inc
			ind(1) = indices(i,1)
			ind(2) = indices(i,2)
			call get_hsdata(ind,lhs(i),deltah(i),deltas(i))
		end do
	end if

	E = 0.d0

	do imol = 1,3
		i = n_monomers(imol)
		do j = 1,nclust_inc
			ij_ind = indices(i,:)+indices(j,:)
			call get_masses_and_radii(m,r,ij_ind)
			if (clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3)).gt.0) then
				ij = clust_from_indices(ij_ind(1),ij_ind(2),ij_ind(3))
				if (lhs(i) .and. lhs(j) .and. lhs(ij)) then
					E(i,j) = K(i,j)*1.01325d5/(1.38065040000000d-23*temperature)*&
							&exp((deltah(ij)-deltah(i)-deltah(j)-temperature*&
							&1.d-3*(deltas(ij)-deltas(i)-deltas(j)))*(5.03218937158374d+02/temperature))
				else
					xmol = real(ij_ind(imol),kind=kind(1.d0))/real(sum(ij_ind),kind=kind(1.d0))
					E(i,j) = K(i,j)*psat(imol)*xmol/(1.38065040000000d-23*temperature)*&
							&exp((2.d0*sigma*monvol(imol))/((1.38065040000000d-23*temperature)*r))
				end if
				if (i .eq. j) then
					E(i,j) = .5d0*E(i,j)
				else if (E(j,i) .gt. E(i,j)) then
					E(i,j) = E(j,i)
				end if
				E(j,i) = E(i,j)
			end if
		end do
	end do


end subroutine get_evap
