program gtritest
use variables
use seq_file_names
implicit none

call cmdline()

#ifdef fuse_namelist
 call read_initfile()
#endif

call init()
call init_gas()
call calculate_grid_parameters()
call update_sigma_z()
call update_sigma_y()

do im=1,mixstep ! loop through mixstep - meteo timestep
  ! print*, im
  do ic=1,NSPEC ! loop through components
    do iy=1,ny
      if (iy==1.and.ic==1.and.im==1) call update_terms_z(sigmaz,dt_mete)
      call gtri(a=termaz, c=termcz, b=termbz, d=Cmajor(ic,:,iy), p=Cmajor1(ic,:,iy), &
      l1=1, a1=a0Gz(ic,iy), q1=3.0d0, lm=1, am=amGz(ic,iy), qm=3.0d0, nz=nz, l2=2)
      ! l1=1, a1=0d0, q1=3.0d0, lm=2, am=amGz(ic,iy), qm=3.0d0, nz=nz, l2=2)
    Cmajor(ic,:,iy) = Cmajor1(ic,:,iy)
    end do
    do iz=2,nz
      if (ic==1.and.im==1) then
        call update_terms_y(sigmahor(iz,:),termay(:,iz),termby(:,iz),termcy(:,iz),dt_mete)
      end if
      if (ic==1.and.im==1.and.iz==5) then
        call save_sigma
      end if
      call gtri(a=termay(:,iz), c=termcy(:,iz), b=termby(:,iz), d=Cmajor(ic,iz,:), p=Cmajor1(ic,iz,:), &
                l1=2, a1=a0G(ic,iz), q1=3.0d0, lm=2, am=amG(ic,iz), qm=3.0d0, nz=ny, l2=2)
    end do
  end do ! end loop through components
  ! Cmajor1(:,1,:) = Cmajor1(:,2,:)
  where (Cmajor1<0d0) Cmajor1=0d0
  Cmajor = Cmajor1
  ! Fill corners
  ! Cmajor(:,1,1) = (Cmajor(:,2,1)+Cmajor(:,1,2))*0.5d0
  ! Cmajor(:,1,ny) = (Cmajor(:,1,ny-1)+Cmajor(:,2,ny))*0.5d0
  ! Cmajor(:,nz,1) = (Cmajor(:,nz-1,1)+Cmajor(:,nz,2))*0.5d0
  ! Cmajor(:,nz,ny) = (Cmajor(:,nz,ny-1)+Cmajor(:,nz-1,ny))*0.5d0
end do ! end loop through mixstep

call dump_gases()

contains

subroutine cmdline()
  implicit none
  ! intialize a plume if necessary
  CALL GET_COMMAND_ARGUMENT(1,f_init,status=ioi)
  if (ioi == 0) read(f_init,*) i_init

  ! 0 to use the input file as output, 1 to advance the index by 1
  CALL GET_COMMAND_ARGUMENT(2,f_init,status=ioi)
  if (ioi == 0) read(f_init,*) advance

#ifdef fuse_namelist
  ! namelist for meteo if used
  CALL GET_COMMAND_ARGUMENT(3,nmlfile,status=ioi)
  if (ioi /= 0) print*,nmlfile, 'Namelist defined, but not found, using defaults'
#endif

end subroutine cmdline

subroutine init()
  implicit none
  integer :: iz

  if (i_init==1) THEN
    print '(a)', 'Grid settings:'
    print '(a, t12,i0)', 'nz', nz
    print '(a, t12,i0)', 'ny', ny
    print '(a, t12,i0)', 'point_z', point_z
    print '(a, t12,i0)', 'point_y', point_y
    print '(a, t12,i0)', 'mixstep', mixstep
    print '(a, t12,f4.2)', 'dt_mete', dt_mete
    print '(a, t12,i0)', 'NSPEC', NSPEC
    print '(a, t12,a)', 'path', path
  end if

  call numbers()

  ! Z = [ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90,100,110,120,130,140,150,160 ,   &
  !  170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,   &
  !  350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,   &
  !  530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,   &
  !  710,720,730,740,750,760,770,780,790]
  !
  Z = [ (iz*dx, iz=0,nz-1) ]
  Y = [ (iy*dx, iy=0,ny-1) ]
  dz = Z(2:)-z(1:nz-1)
  dy = Y(2:)-Y(1:ny-1)
  daz = 1
  day = 1
  dbz = 1
  dby = 1
  dcz = 1
  dcy = 1

  i_file = findex(TRIM(path)//'/01/01/CFINAL_',4,'.r16')
  if (advance==0) i_file = max(0,i_file)
  print '(a,i0.4,a)', 'reading from '//TRIM(path)//'/*/*/CFINAL_',i_file,'.r16'
  print '(a,i0.4,a)', 'will write to '//TRIM(path)//'/*/*/CFINAL_',i_file+advance,'.r16'

end subroutine init


subroutine init_gas()
  implicit none
  integer :: i,j
  if (initW_Previous) then
    do i=1,nz
    do j=1,ny
    call readBin(ch_gas,i,j)
    Cmajor(:,i,j) = ch_gas
    end do
    end do
  else
    Cmajor = 1D0
  end if

  if (i_init==1) THEN
    print*, 'Intialization of plume'
    C = 1d0
    C(point_z,point_y) = point_multi
    do i=1,NSPEC
      if (i==2201.or.i==2206.or.i==2209.or.i==2207.or.i==1) THEN
        Cmajor(i,:,:) = Cmajor(i,:,:) * C
      end if
    end do
  ELSE
    print*, 'Mixing input'
  END IF

end subroutine init_gas

subroutine dump_gases()
  implicit none
  integer :: i,j
  do i=1,nz
    do j=1,ny
      call writeBin(Cmajor(:,i,j),i,j)
    end do
  end do

end subroutine dump_gases

subroutine calculate_grid_parameters()
  implicit none
  integer :: k

  do k=2,nz-1
     daz(k)=dz(k+1)*(dz(k)+dz(k+1))
     dbz(k)=dz(k)*dz(k+1)
     dcz(k)=dz(k)*(dz(k)+dz(k+1))
     dfz(k)=(dz(k+1)+dz(k))/dz(k+1)/4.0d0
  enddo

  dfz(1)=dfz(2)

  do k=2,ny-1
     day(k)=dy(k+1)*(dy(k)+dy(k+1))
     dby(k)=dy(k)*dy(k+1)
     dcy(k)=dy(k)*(dy(k)+dy(k+1))
     dfy(k)=(dy(k+1)+dy(k))/dy(k+1)/4.0d0
  enddo
  dfy(1)=dfy(2)

  ! sigmaz0 = [(0.05, iz=1,nz)]
  ! sigmaz0(:nz/3) = sigmaz0(:nz/3) + [(sin(1d0*iz/(nz/3-1)*3.141593)**2, iz=0,nz/3-1)]
  ! Grisogono scheme
  sigmaz0 = MAX(1d-2,0.39*friction_vel*z*exp(-0.5*(z/(0.21*pblh))**2.))

end subroutine calculate_grid_parameters

subroutine update_sigma_z()
implicit none
integer :: i
if (randomize_k==1) &
  call random_number(random_z)
sigmaz = sigmaz0 * ( 0.5 + random_z ) ! * k0 * sigmaK
where (sigmaz<0.01)
  sigmaz = 0.01
end where
! sigmaz(1:2) = 0d0
! sigmaz(:) = 3.0
! do i=1,nz
! print*, sigmaz(i),','
! end do
! stop
! sigmaz(70:) = 0d0

end subroutine update_sigma_z


subroutine update_sigma_y()
implicit none
integer :: iz
if (randomize_k==1) &
  call random_number(random_y)
do iz=1,nz
  sigmahor(iz,:) = KyKz * ( 0.5 + random_y ) * sigmaz(iz)
end do
where (sigmahor<0.01)
  sigmahor = 0.01
end where
end subroutine update_sigma_y

subroutine save_sigma()
  implicit none
  character(len=120) :: fname
  write(fname,'(a,i0.4,a)') TRIM(path)//'/../k-values/K-values_',i_file,'.txt'
  open(unit=123,file=TRIM(fname),status='replace',action='write')
  write(123,*) sigmaz
  write(123,*) sigmahor(5,:)
  write(123,*) pblh
  write(123,*) friction_vel
  write(123,*) KyKz
  write(123,*) randomize_k
  close(123)
end subroutine save_sigma

subroutine update_terms_z(kt,dt)
  implicit none
  real(dp), INTENT(in) :: kt(:), dt

  Khtz(2:nz-1) = dfz(2:nz-1)*kt(3:nz) + (1.0d0-dfz(2:nz-1))*kt(2:nz-1)
  Khbz(2:nz-1) = (1.0d0-dfz(1:nz-2))*kt(1:nz-2) + dfz(1:nz-2)*kt(2:nz-1)
  DO k=2,nz-1
    termaz(k) = dt*(2.0d0*(Khtz(k)))/daz(k)
    termcz(k) = dt*(2.0d0*(Khbz(k)))/dcz(k)
    termbz(k) = dt*(2.0d0*(dz(k)*(Khtz(k)) + dz(k+1)*(Khbz(k))) /(dz(k+1)+dz(k)) )/dbz(k)+ 1.0d0
  ENDDO

end subroutine update_terms_z

subroutine update_terms_y(kt,termay,termby,termcy,dt)
  implicit none
  real(dp), INTENT(in) :: kt(:), dt
  real(dp), INTENT(out) :: termay(:),termcy(:),termby(:)

  Khty(2:ny-1) = dfy(2:ny-1)*kt(3:ny) + (1.0d0-dfy(2:ny-1))*kt(2:ny-1)
  Khby(2:ny-1) = (1.0d0-dfy(1:ny-2))*kt(1:ny-2) + dfy(1:ny-2)*kt(2:ny-1)
  DO k=2,ny-1
    termay(k) = dt*(2.0d0*(Khty(k)))/day(k)
    termcy(k) = dt*(2.0d0*(Khby(k)))/dcy(k)
    termby(k) = dt*(2.0d0*(dy(k)*(Khty(k)) + dy(k+1)*(Khby(k))) /(dy(k+1)+dy(k)) )/dby(k)+ 1.0d0
  ENDDO

end subroutine update_terms_y
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine to solve the vertical exchange of gases - Thomas algorithm
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   a	=	array of matrix elements left of the diagonal
!   c	=	array of matrix diagonal elements
!   b	=	array of matrix elements right of the diagonal
!   d	=	input right-hand-side vector
!   p	=	output vector
!   l1	=	lower boundary condition option indicator
!   a1	=	lower boundary value
!   q1	=	lower boundary value
!   lm	=	upper boundary condition option indicator
!   am	=	upper boundary value
!   qm	=	upper boundary value
!   nz	=	upper level number
!   l2	=	lower level start number

SUBROUTINE gtri(a, c, b, d, p, l1, a1, q1, lm, am, qm, nz, l2)
  IMPLICIT NONE
  ! interface:
  INTEGER, INTENT(in) :: l1,l2,lm,nz
  REAL(kind=dp), INTENT(in) :: a1,q1,am,qm
  REAL(kind=dp), INTENT(in) :: a(nz),b(nz),c(nz),d(nz)
  REAL(kind=dp), INTENT(inout) :: p(nz)

  ! local:
  INTEGER :: m2,m,n
  REAL(kind=dp) :: den
  REAL(kind=dp) :: e(nz),f(nz)

  e=0d0
  f=0d0

  IF(l1.EQ.1) e(l2-1)=0.0d0
  IF(l1.EQ.1) f(l2-1)=a1
  IF(l1.EQ.2) e(l2-1)=1.0d0
  IF(l1.EQ.2) f(l2-1)=-a1
  IF(l1.EQ.3) e(l2-1)=a1/(a1-1.0d0)
  IF(l1.EQ.3) f(l2-1)=q1/(1.0d0-a1)

  m2=nz-1

  DO m=l2,m2
    den=b(m)-c(m)*e(m-1)
    f(m)=(d(m)+c(m)*f(m-1))/den
    e(m)=a(m)/den
  ENDDO

  IF(lm.EQ.1) p(nz)=am
  IF(lm.EQ.2) p(nz)=(f(m2)+am)/(1.0d0-e(m2))
  IF(lm.EQ.3) p(nz)=(f(m2)+qm/am)/((1.0d0+am)/am-e(m2))

  DO n=1,m2-l2+2
    m=nz-n
    p(m)=e(m)*p(m+1)+f(m)
  ENDDO
END SUBROUTINE gtri

!-------------------------------------------
subroutine readBin(X,z,y)
implicit none
integer,intent(in) :: z,y
real(dp),intent(out) :: X(:)
character(len=128) :: fname
integer            :: ioi
write(fname,'(a,i0.4,a)') TRIM(path)//'/'//numb(z)//'/'//numb(y)//'/CFINAL_',i_file,'.r16'

open(unit=666,FILE=TRIM(fname),action='read',form='unformatted', ACCESS="STREAM",iostat=ioi)
if (ioi==0) then
  read(666, POS=1) X
  close(666)
end if
end subroutine readBin

!-------------------------------------------
subroutine writeBin(X,z,y)
implicit none
integer,intent(in) :: z,y
real(dp),intent(in) :: X(:)
character(len=128) :: fname
write(fname,'(a,i0.4,a)') TRIM(path)//'/'//numb(z)//'/'//numb(y)//'/CFINAL_',i_file+advance,'.r16'
open (unit=1,file=TRIM(fname),form='unformatted', access='direct',recl=NSPEC*8)
write (1,rec=1) X
close(1)
end subroutine writeBin

subroutine numbers
  implicit none
  character(len=2),parameter :: bfr(99) = &
  ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20',&
  '21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40',&
  '41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60',&
  '61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80',&
  '81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99']
  numb = bfr(:ny)


end subroutine numbers

real(dp) function amG(ic,iz)
implicit none
integer :: ic,iz
amG = (Cmajor(ic,iz,ny-1)-Cmajor(ic,iz,ny-2))/10.0
end function amG

real(dp) function a0G(ic,iz)
implicit none
integer :: ic,iz
a0G = (Cmajor(ic,iz,3)-Cmajor(ic,iz,2))/10.0
end function a0G

real(dp) function amGz(ic,iy)
implicit none
integer :: ic,iy
! amGz = (Cmajor(ic,nz-1,iy)-Cmajor(ic,nz-2,iy))/10.0
amGz = 2d0*Cmajor(ic,nz-1,iy) - Cmajor(ic,nz-2,iy)
end function amGz

real(dp) function a0Gz(ic,iy)
implicit none
integer :: ic,iy
! a0Gz = (Cmajor(ic,2,iy)-Cmajor(ic,3,iy))/5.0
a0Gz = 2d0*Cmajor(ic,2,iy) - Cmajor(ic,3,iy)
end function a0Gz



end program gtritest
