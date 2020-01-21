MODULE aerosol_dynamics

USE second_precision, ONLY: dp
USE constants
USE INPUT
USE ParticleSizeDistribution
USE Aerosol_auxillaries

IMPLICIT NONE


Logical :: RH_DEPENDENCE=.True.
Logical :: real_code =.TRUE.
Logical :: use_raoult =.True.
Logical :: use_new_method=.true.

! integer, parameter::nr_bins=100
! integer, parameter::nr_cond=10

CONTAINS


SUBROUTINE Condensation_apc(timestep, vapour_properties, particle_properties, conc_pp,ch_gas, &
                          CH2SO4,ambient, dmass) ! Add more variables if you need it


  REAL(dp), INTENT(IN) :: timestep


  REAL(dp), DIMENSION(nr_species_p), INTENT(INOUT) :: ch_gas  ! [molec m-3], condensing vapour concentrations, which is H2SO4 and organics (ELVOC)
  REAL(dp), INTENT(INOUT) :: cH2SO4
  REAL(dp), DIMENSION(n_bins_particle) :: particle_volume_new, particle_conc_new,particle_conc_old
  REAL(dp), DIMENSION(n_bins_particle):: vp, particle_conc
  REAL(dp), DIMENSION(n_bins_particle,nr_cond) :: CR
  REAL(dp), DIMENSION(n_bins_particle) :: CR_H2SO4


  INTEGER :: j,ii,ip, a,ig
  type(ambient_properties), intent(in) :: ambient
  type(vapour_ambient), intent(in) :: vapour_properties
  type(PSD), intent(inout) :: particle_properties

  Real(dp) :: r1, r2
  Real(dp) :: sum_org
  Real(dp) :: temperature, pressure
  Real(dp), dimension(n_bins_particle,nr_species_p) :: Kelvin_Effect, kohler_effect,conc_pp_fixed
  Real(dp), dimension(n_bins_particle,nr_species_p),intent(inout) :: conc_pp

  Real(dp), dimension(n_bins_particle,nr_species_p) :: conc_pp_old
  Real(dp), dimension(n_bins_particle,nr_species_P) :: xorg                       ! mole fraction of each organic compound in each bin

  Real(dp), dimension(nr_species_p) :: conc_guess, c_total!, ch_gas_old
  Real(dp), dimension(nr_bins) :: conc_H2SO4_guess, conc_H2SO4_old, conc_H2SO4_pp, conc_H2SO4
  Real(dp), dimension(nr_bins, nr_species_P) :: composition_old
  Real(dp), dimension(nr_bins, nr_species_P), intent(inout) :: dmass
  real(dp) :: test(nr_bins)

  temperature=ambient%temperature
  pressure  = ambient%pressure
  particle_conc=particle_properties%conc_fs
  particle_conc_old=particle_properties%conc_fs
  conc_H2SO4_old = cH2SO4 ! need to replace this with the concentation of H2SO4 from chemistry
  conc_H2SO4   = cH2SO4
  composition_old = particle_properties%composition_fs
 ! write(*,*), 'ch_gas', ch_gas(ind_H2SO4)
  ! conc_H2SO4_old=conc_H2SO4
  ! Add more variabels as you need it...

if (real_code) THEN


   conc_pp_old = conc_pp
   c_total = ch_gas + sum(conc_pp,1)    ! vapour phase + particle phase


 if (use_raoult) THEN
   do ii=1,n_bins_particle
        if (particle_properties%volume_fs(ii) >0.0 ) then
         sum_org = sum(conc_pp(ii,:),DIM=1, Mask=(vapour_properties%cond_type==1))
          ! write(*,*) 'sum_org', sum_org
         if (sum_org>0) THEN
           xorg(ii,:)= conc_pp(ii,:)/sum_org
         else
           xorg(ii,:) = 1D0
         end if
        else
         xorg(ii,:) = 1D0
       end if
    end do
 end if

 do ii=1, n_bins_particle
   if (particle_properties%diameter_fs(ii) <=  1*1D-9) THEN
       ! write(*,*) 'In here',ii
     xorg(ii,:) = 1D0
   end if
 end do


 ! collision rate, kelvin_eff, kohler_effect, and apc_scheme
 do ii=1, nr_species_p

   if(vapour_properties%cond_type(ii)==1) then ! for all organics, No acids included

     CR(:,ii)  = particle_conc *collision_rate(ii,current_PSD,ambient,vapours,Diff_org)! s-1

   ! kelvinEffect
     Kelvin_Effect(:,ii) = 1D0 + 2D0*vapour_properties%surf_tension(ii)*vapour_properties%molec_volume(ii) / &
                               (kb*temperature*particle_properties%diameter_fs/2D0)     ! unit less

     kohler_effect(:,ii) = Kelvin_Effect(:,ii)*xorg(:,ii)

  !!!! apc scheme here
     conc_guess(ii) = (ch_gas(ii) + timestep*sum(CR(:,ii)*kohler_effect(:,ii)*vapour_properties%c_sat(ii))) / &
                   (1D0 + timestep*sum(CR(:,ii)))

     conc_pp(:,ii) = conc_pp_old(:,ii) + timestep*CR(:,ii)*(MIN(conc_guess(ii),c_total(ii)) - &
                                 kohler_effect(:,ii)*vapour_properties%c_sat(ii))

     WHERE( conc_pp(:,ii)<0D0) conc_pp(:,ii)=0D0 ! PREVENT OVERESTIMATION OF EVAPORATION

     ch_gas(ii) = c_total(ii) - sum(conc_pp(:,ii)) ! update ch_gas, total conc - particle phase concentration

     IF (ch_gas(ii)<0d0) then
       ch_gas(ii)=0D0
     END IF

   else ! acids mainly h2so4 for now
     CR_H2SO4  = particle_conc *collision_rate(ii,current_PSD,ambient,vapours,Diff_org)! s-1

   ! kelvinEffect
     Kelvin_Effect(:,ii) = 1D0 + 2D0*vapour_properties%surf_tension(ii)*vapour_properties%molec_volume(ii) / &
                               (kb*temperature*particle_properties%diameter_fs/2D0)     ! unit less



  !!!! apc scheme here
     conc_H2SO4_guess = (conc_H2SO4_old + timestep*sum(CR_H2SO4(:)*kelvin_effect(:,ii)*vapour_properties%c_sat(ii))) / &
                   (1D0 + timestep*sum(CR_H2SO4))

     conc_H2SO4_pp = conc_pp_old(:,ii) + timestep*CR_H2SO4(:)*(MIN(conc_H2SO4_guess,c_total(ii)) - &
                                 kohler_effect(:,ii)*vapour_properties%c_sat(ii))

     WHERE( conc_H2SO4_pp<0D0) conc_H2SO4_pp=0D0 ! PREVENT OVERESTIMATION OF EVAPORATION

     cH2SO4 = c_total(ii) - sum(conc_H2SO4_pp) ! gas phase #/m3

     conc_pp(:,ii) = conc_H2SO4_pp
   end if

 end do

! update composition
 do  ii = 1, particle_properties%nr_bins
         particle_properties%composition_fs(ii,:) = conc_pp(ii,:) * vapour_properties%molar_mass  / Na !&
                                                   ! / particle_conc(ii)
 end do

 !!! update dmass

   do ii=1,nr_bins
    dmass(ii,:) = (particle_properties%composition_fs(ii,:) - composition_old(ii,:) ) !/ vapour_properties%density &
                  ! /particle_conc(ii)
   end do

if (use_new_method) then
    CALL Mass_Number_Change('condensation')
else

do ig=1,nr_bins!_particle
   particle_volume_new(ig) = &
                             !  sum(conc_pp(ig,:)/Na*vapour_properties%molar_mass/&
                             ! vapour_properties%density/particle_conc(ig))

            sum(particle_properties%composition_fs(ig,:) / vapour_properties%density) ! another way using composition

  test(ig) =particle_properties%volume_fs(ig)/particle_conc(ig) + &
  SUM(dmass(ig,:) / vapour_properties%density(:))

end do
! write(*,*) 'particle_volume_new', sum(particle_volume_new), '', sum(test),'', &
    ! sum(particle_properties%volume_fs)!, sum(particle_properties%composition_fs)

 vp = particle_properties%volume_fs !diameter**3D0*pi/6D0
 particle_conc_new=0D0 ! Initialise a new vector with the new particle concentrations
 particle_conc_new(n_bins_particle)=particle_properties%conc_fs(n_bins_particle)
 conc_pp_fixed = 0D0


 ! Update the particle concentration in the particle_conc vector:
 ! Taken from Pontus
 DO j = 1,n_bins_particle
    a = MINLOC(vp-particle_volume_new(j),1,mask = (vp-particle_volume_new(j)) >= 0D0)
       ! write(*,*) 'a=',a,'j=',j
    IF (a > j) THEN ! growth
      r1 = (vp(a)-particle_volume_new(j))/(vp(a)-vp(a-1)) ! Fraction of particles in size bin a-1
      r2 = 1D0-r1 ! Fraction of particles in next size bin (a)
      IF (a > n_bins_particle+1) THEN
      ELSE IF (a == n_bins_particle+1) THEN
        particle_conc_new(a-1) = particle_conc_new(a-1)+r1*particle_conc(j)
        conc_pp_fixed(a-1,:) = conc_pp_fixed(a-1,:)+vp(a-1)/particle_volume_new(j)*r1*conc_pp(j,:)
      ELSE
        particle_conc_new(a-1) = particle_conc_new(a-1)+r1*particle_conc(j)
        particle_conc_new(a) = particle_conc_new(a)+r2*particle_conc(j)
        conc_pp_fixed(a-1,:) = conc_pp_fixed(a-1,:)+vp(a-1)/particle_volume_new(j)*r1*conc_pp(j,:)
        conc_pp_fixed(a,:) = conc_pp_fixed(a,:)+vp(a)/particle_volume_new(j)*r2*conc_pp(j,:)
      END IF
    ELSE ! evaporation
      IF (a == 0) THEN
      ELSEIF (a == 1) THEN
        r1 = (particle_volume_new(j))/(vp(1)) ! Fraction of particles in size bin a
        particle_conc_new(1) = particle_conc_new(1)+r1*particle_conc(j)
        conc_pp_fixed(1,:) = conc_pp_fixed(1,:)+vp(1)/particle_volume_new(j)*r1*conc_pp(j,:)
      ELSE
        r1 = (particle_volume_new(j)-vp(a-1))/(vp(a)-vp(a-1)) ! Fraction of particles in size bin i
        r2 = 1-r1 ! Fraction of particles in previous size bin (i-1)

        particle_conc_new(a) = particle_conc_new(a)+r1*particle_conc(j)
        particle_conc_new(a-1) = particle_conc_new(a-1)+r2*particle_conc(j)
        conc_pp_fixed(a,:) = conc_pp_fixed(a,:)+vp(a)/particle_volume_new(j)*r1*conc_pp(j,:)
        conc_pp_fixed(a-1,:) = conc_pp_fixed(a-1,:)+vp(a-1)/particle_volume_new(j)*r2*conc_pp(j,:)
      END IF
      END IF
 END DO

 particle_conc=particle_conc_new
 conc_pp=conc_pp_fixed

! write(*,*) 'after',sum(conc_pp,2)

  ! if particle_conc is <0 . Is possible for full stationary method
  do ip=1, n_bins_particle
    if (particle_conc(ip)<=0d0) then
      ! write(*,*) 'ip', ip
      particle_conc(ip) = 1d0
      ! particle_properties%conc_fs(ip)=particle_conc(ip)
       conc_pp(ip,:)=conc_pp_old(ip,:)*particle_properties%conc_fs(ip)/&
                                         particle_conc_old(ip)
    end if
  end do

particle_properties%conc_fs=particle_conc

! do  ii = 1, particle_properties%nr_bins
!         particle_properties%composition_fs(ii,:) = conc_pp(ii,:) * vapour_properties%molar_mass  / Na &
!                                                   / particle_conc(ii)
! end do
end if ! if use_new_method
end if



END SUBROUTINE Condensation_apc



function collision_rate(jj,particle_properties,ambient,vapour_properties,Diff_org)

    implicit none

    type(ambient_properties), intent(in) :: ambient
    type(vapour_ambient), intent(in) :: vapour_properties
    type(PSD), intent(in) :: particle_properties

    integer, intent(in) :: jj                      ! the 'number' of condensing vapour
    real(dp), dimension(nr_cond), intent(in) :: Diff_org
    real(dp),dimension(n_bins_particle) :: collision_rate

    ! integer :: ii

    real(dp) :: air_free_path,temp,rh,pres,viscosity !dif_vap, r_vap, n1, , r_hyd, r_h2o,zero
    real(dp), dimension(n_bins_particle) :: knudsen, Diff_par, slip_correction


    REAL(DP):: Diff_H2SO4_0    ! gas diffusivity for H2SO4 AT RH=0%

    real(dp) :: Diff_H2SO4,Keq1,Keq2,d_H2SO4,mH2SO4,speedH2SO4,&
                dorg,dens_air,speedorg

    real(dp),dimension(n_bins_particle) :: gasmeanfpH2SO4,speed_p,KnH2SO4,f_corH2SO4,DH2SO4eff,gasmeanfporg,&
                                        Knorg,f_cororg,Dorgeff

    temp  = ambient%temperature
    rh    = ambient%rh              ! RH unit %,values from 1-100
    pres  = ambient%pressure



    ! viscosity of air, density oif air and mean free path in air
    viscosity     = 1.8D-5*(temp/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air      = Mair*pres/(R*temp)
    air_free_path = 2D0*viscosity/(pres*SQRT(8D0*Mair/(pi*R*Temp))) ! gas mean free path in air

    ! knudsen number
    knudsen = 2D0 * air_free_path/particle_properties%diameter_fs

    ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34 pg 407)
    slip_correction = 1D0 + knudsen * &
                      (1.257D0 + 0.4D0*exp(-1.1D0/(2D0*air_free_path/particle_properties%diameter_fs)))

    ! particle diffusion coefficient (m^2 s^-1)
    Diff_par = (kb * temp * slip_correction) / (3D0 * pi * viscosity * particle_properties%diameter_fs)



    where (particle_properties%particle_mass_fs> 0)
      speed_p = SQRT(8D0*kb*temp/(pi*particle_properties%particle_mass_fs)) ! speed of particle unit(mass)=kg
    end where
   ! write(*,*) 'Carlton debug check molecvol line 276', molecvol(jj)

  IF (RH_DEPENDENCE) THEN
    IF (vapour_properties%cond_type(jj)==2) THEN
     ! For testing. In the main model we use RH dependence
     ! H2SO4 mass transfer rate (m3/s) ,RH dependent Diffusion, diameter and mass

     ! diffusivity of H2SO4 at RH=0%
      Diff_H2SO4_0 = 0.09D-4

      Keq1=0.13D0      ! Hanson and Eisele
      Keq2=0.016D0     ! Hanson and Eisele

    ! Diffusivity H2SO4 at ambient RH
       Diff_H2SO4 =  (Diff_H2SO4_0 + 0.85D0*Diff_H2SO4_0*Keq1*RH + 0.76D0*Diff_H2SO4_0*Keq1*Keq2*RH**2D0) &
        /(1D0 + Keq1*RH + Keq1*Keq2*RH**2D0)

       ! RH dependent H2SO4 diameter
       d_H2SO4 = (((vapour_properties%molar_mass(jj)/vapour_properties%density(jj)/Na)*6D0/pi)**(1D0/3D0) + &
          Keq1*RH*(((vapour_properties%molar_mass(jj) + 18D-3)/vapour_properties%density(jj)/Na)*6D0/pi)**(1D0/3D0) + &
          Keq1*Keq2*RH**2D0*(((vapour_properties%molar_mass(jj) + 36D-3)/vapour_properties%density(jj)/Na)*6D0/pi)**(1D0/3D0)) / &
          (1D0 + Keq1*RH + Keq1*Keq2*RH**2D0)
      !

      !  ! RH H2SO4 molecular mass
       mH2SO4 = (vapour_properties%molar_mass(jj)/Na + Keq1*RH*(vapour_properties%molar_mass(jj) + 18D-3)/Na +  &
                Keq1*Keq2*rh**2D0*(vapour_properties%molar_mass(jj) + 36D-3)/Na) / &
                (1D0 + Keq1*RH + Keq1*Keq2*RH**2D0)

      ! write(*,*) 'line 320 mH2so4', mH2SO4,'',d_h2so4,'',vapour_properties%molar_mass(jj),vapour_properties%molec_mass(jj)
      !
      !  ! speed of H2SO4 molecule
       speedH2SO4 = SQRT(8D0*kb*temp/(pi*mH2SO4))
      !  ! gasmeanfpH2SO4 = 0
      !
      !  ! collision of H2SO4 molecule with particle. Only for bins where mass > 0
        gasmeanfpH2SO4 = 3D0*(Diff_H2SO4 + Diff_par)/SQRT(speedH2SO4**2D0 + speed_p**2D0)
      !
      ! !  Knudsen number H2SO4
       KnH2SO4=2D0*gasmeanfpH2SO4/(particle_properties%diameter_fs + d_H2SO4)
      !
      ! ! Fuchs-Sutugin correction factor for transit
       f_corH2SO4 = (0.75*vapour_properties%alpha(jj)*(1D0 + KnH2SO4))/ &
                    (KnH2SO4**2D0 + KnH2SO4 + 0.283*KnH2SO4*vapour_properties%alpha(jj) + 0.75*vapour_properties%alpha(jj))
      !
       DH2SO4eff = (Diff_H2SO4 + Diff_par)*f_corH2SO4    !m2/s

       collision_rate = 2D0*pi*(particle_properties%diameter_fs + d_H2SO4)*DH2SO4eff  !mass transfer rate m3/s

       ! write(*,*) 'mH2SO4', mH2SO4

! Knudsen number  for organic condensable compounds
    ELSE
      !............................................................................
      !Organics mass transfer rate (m3/s)


       ! ! diamter of organic compounds
       dorg= (6D0*vapour_properties%molec_volume(jj)/pi)**(1D0/3D0)               !estimated diameter (m)
       !
       ! ! diffusivity organic compound
       ! Diff_org=5D0/(16D0*Na*dorg**2D0*dens_air)*&
       !    sqrt(R*temp*Mair/(2D0*pi)*((vapour_properties%molar_mass(jj) + Mair)/vapour_properties%molar_mass(jj)))
       !
       speedorg=SQRT(8D0*Kb*temp/(pi*vapour_properties%molec_mass(jj))) !spped of organic molecules

       ! gasmeanfporg=0
       gasmeanfporg=3D0*(Diff_org(jj) + Diff_par)/SQRT(speedorg**2D0 + speed_p**2D0)

       ! Knudsen number organic comp
       Knorg=2D0*gasmeanfporg/(particle_properties%diameter_fs + dorg)

       ! Fuchs-Sutugin correction factor for transit
       f_cororg=(0.75*vapour_properties%alpha(jj)*(1D0 + Knorg))/&
                 (Knorg**2D0 + Knorg+0.283*Knorg*vapour_properties%alpha(jj) + 0.75*vapour_properties%alpha(jj))

       Dorgeff = (Diff_org(jj) + Diff_par)*f_cororg                    ! m^2/s

       collision_rate = 2D0*pi*(particle_properties%diameter_fs + dorg)*Dorgeff        ! mass transfer coefficient s^-1
        ! write(*,*) 'line 340 collision_rate of org', sum(collision_rate)

     ENDIF



   ENDIF

end function collision_rate
!
! !
!
END MODULE aerosol_dynamics

!! for testing
! SUBROUTINE Nucleation_apc (dt, nucl_gas, particle_conc)! (Add input and output variables here)
!
!   ! Consider how kinetic H2SO4 nucleation influence the number concentrations of particles
!   ! in the fist size bin particle_conc(1) within one model time step
!   real(dp), intent(in   ) :: dt      , &  ! [s], integration time step
!                              nucl_gas     ! [molec m-3], concentration of condensed vapors
!   real(dp), intent(  out) :: particle_conc(nr_bins)  ! [molec m-3], particle number concentration
!   real(dp) :: Knucl            ! [m3 molec-1 s-1], nucleation coefficient
!   real(dp) :: nucleation_rate  ! [molec m-3 s-1], nucleation rate
!   ! real(dp),dimension(nr_bins, nr_cond), intent(inout) ::  conc_pp
!   ! real(dp),dimension(nr_cond), intent(in) :: c_p_nucl
!   ! real(dp),dimension(nr_bins),intent(inout):: diameter
!
!
!   Knucl = 1.0e-20_dp  ! [m3 molec-1 s-1]
!
!   nucleation_rate = Knucl * nucl_gas**2
!   particle_conc(1)    = particle_conc(1) + nucleation_rate*dt
!
!   ! write(*,*) 'in nuc before update',sum(conc_pp)
!   ! from Pontus.. update the conc_pp in particle phase and diameter
!
!   ! vp = sum(conc_pp(1,:)/Na*molar_mass/density/particle_conc(1)) !particle_volume(1)*particle_conc(1)*density/molecular_mass
!   ! ! write(*,*) 'in nuc sfter update',sum(conc_pp)
!   ! ! write(*,*) vp,'', sum(conc_pp)
!   ! diameter(1) = (vp*6D0/pi)**(1D0/3D0)
!
! END SUBROUTINE Nucleation_apc
