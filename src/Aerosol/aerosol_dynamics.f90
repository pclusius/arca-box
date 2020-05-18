MODULE aerosol_dynamics

USE second_precision, ONLY: dp
USE constants
USE INPUT
USE ParticleSizeDistribution
USE Aerosol_auxillaries
use omp_lib

IMPLICIT NONE

CONTAINS


SUBROUTINE Nucleation_routine (J_TOTAL,particle_conc,composition, dt_nuc)
  IMPLICIT NONE
  REAL(dp), intent(in   ) :: J_TOTAL ![m-3 s-1]
  REAL(dp), intent(inout) :: particle_conc  ! [molec m-3], particle number concentration
  REAL(dp), intent(inout) :: composition(:)  ! [molec m-3], particle number concentration
  REAL(dp), INTENT(IN) :: dt_nuc ! integration timestep for condensation [s]

  integer :: i
  ! REAL(dp) :: Knucl            ! [m3 molec-1 s-1], nucleation coefficient
  ! REAL(dp) :: nucleation_rate  ! [molec m-3 s-1], nucleation rate

  ! Knucl = 1.0e-20_dp  ! [m3 molec-1 s-1]
  ! nucleation_rate=Knucl*CH_h2so4 **2
  ! particle_conc(1) = old_par + nucleation_rate * GTIME%dt

  particle_conc = particle_conc  + J_TOTAL * dt_nuc
  composition(:) = VAPOUR_PROP%mfractions * pack(get_volume(),[(i<2,i=1,n_bins_particle)]) * VAPOUR_PROP%density
END SUBROUTINE Nucleation_routine



! ======================================================================================================================
! Main condensation subroutine. We use the Analytical predictor of condensation scheme Jacobson 1997c, 2002 and
! fundamentals of atmospheric modelling. Dmass, which is the flux onto or removed from the particles is the outcome
! of this subroutine which is fed to PSD
! ======================================================================================================================
SUBROUTINE Condensation_apc(vapour_prop, conc_vap, dmass, dt_cond, d_dpar,d_vap)
  IMPLICIT NONE
  type(vapour_ambient), INTENT(IN) :: vapour_prop ! Properties of condensing vapours
  ! type(PSD)           , INTENT(IN) :: particles ! Current PSD properties
  REAL(dp), INTENT(INOUT) :: conc_vap(:)  ! [#/m^3], condensing vapour concentrations, DIM(n_cond_tot)
  REAL(dp), INTENT(INOUT) :: dmass(:,:)   ! [kg/m^3] change of mass per particle in particle phase, DIM(n_bins_particle, n_cond_tot)
  REAL(dp) :: n_conc(n_bins_particle), diameter(n_bins_particle), volume(n_bins_particle), mass(n_bins_particle) ! [#/m^3] XXX WORKING DUMMIES
  REAL(dp) :: conc_pp(n_bins_particle,n_cond_tot) ! [#/m^3] Particle phase concentrations, DIM(n_bins_particle,n_cond_tot)
  REAL(dp), DIMENSION(n_bins_particle,n_cond_tot) :: CR ! [/s] collisions per second, collision rate * concentration

  REAL(dp), DIMENSION(n_bins_particle,n_cond_tot) :: Kelvin_Effect, kohler_effect
  REAL(dp), DIMENSION(n_bins_particle,n_cond_tot) :: conc_pp_old
  REAL(dp), DIMENSION(n_bins_particle,n_cond_tot) :: xorg ! mole fraction of each organic compound in each bin
  REAL(dp), DIMENSION(n_bins_particle, n_cond_tot) :: conc_pp_eq
  REAL(dp), DIMENSION(n_cond_tot) :: conc_tot
  REAL(dp), DIMENSION(n_cond_tot) :: conc_vap_old
  REAL(dp) :: conc_guess
  REAL(dp) :: sum_org ! Sum of organic concentration in bin, transient variable
  REAL(dp), INTENT(IN) :: dt_cond ! integration timestep for condensation [s]
  REAL(dp), INTENT(INOUT) :: d_dpar(:) ! relative change in particle diameter [1]
  REAL(dp), INTENT(INOUT) :: d_vap(:) ! relative change in vapor concentrations [1]

  INTEGER :: ii

  conc_pp = 0d0
  conc_pp = get_composition()
  n_conc = get_conc()
  diameter = get_dp()
  volume = get_volume()
  mass = sum(conc_pp, 2)

  ! Fill the number concentration composition matrix
  do ii=1, n_bins_particle
    conc_pp(ii,:) = conc_pp(ii,:) * Na / vapour_prop%molar_mass * n_conc(ii)
  end do

  ! original gas phase concentration
  conc_vap_old =  conc_vap
  ! original particle phase concentration
  conc_pp_old = conc_pp
  ! vapour phase + particle phase
  conc_tot = conc_vap + sum(conc_pp,1)

  if (use_raoult) THEN
    DO ii=1,n_bins_particle
      if (volume(ii) >0.0 .and. diameter(ii) > 1D-9) then
        sum_org = sum(conc_pp(ii,:),DIM=1, Mask=(vapour_prop%cond_type==1))
        if (sum_org>0) THEN
          xorg(ii,:)= conc_pp(ii,:)/sum_org
        else
          xorg(ii,:) = 1D0
        END if
      else
        xorg(ii,:) = 1D0
      END if
    END DO
  END if

  ! Kelvin and Kohler factors.
  DO ii = 1, n_cond_tot
    ! Kelvin factor takes into account the curvature of particles. Unitless
    Kelvin_Effect(:,ii) = 1d0

    Kelvin_Effect(:,ii) = 1D0 + 2D0*vapour_prop%surf_tension(ii)*vapour_prop%molar_mass(ii) &
                        / (R*GTEMPK*vapour_prop%density(ii)*diameter/2D0)

    ! Kohler factor the solute partial pressure effect. Unitless
    kohler_effect(:,ii) = Kelvin_Effect(:,ii)*xorg(:,ii)
  END DO
  ! Treat H2SO4 specially
  kohler_effect(:,n_cond_tot) = Kelvin_Effect(:,n_cond_tot)

  ! Approximate equilibrium concentration (#/m^3) of each compound in each size bin
  DO ii=1,n_bins_particle
    conc_pp_eq(ii,1:n_cond_tot-1) = conc_vap_old(1:n_cond_tot-1)*SUM(conc_pp_old(ii,1:n_cond_tot-1)) &
                                    / (Kelvin_Effect(ii,1:n_cond_tot-1)*vapour_prop%c_sat(1:n_cond_tot-1))
  END DO

  DO ii = 1, n_cond_tot
    ! Total number of collisions/s
    if (Use_atoms) THEN
         CR(:,ii) = n_conc * collision_rate(ii,diameter,mass,vapour_prop)
    ELSE
         CR(:,ii) = n_conc * collision_rate_uhma(ii,diameter,mass,vapour_prop)
    END IF


    if ((ii == vapour_prop%ind_H2SO4) .and. GTIME%printnow) print FMT10_CVU, 'CS SA: ',sum(CR(:,ii)), ' [/s]'
    ! apc scheme here. NOTE that for sulfuric kohler_effect = Kelvin_Effect
    conc_guess = (conc_vap(ii) + dt_cond*sum(CR(:,ii)*kohler_effect(:,ii)*vapour_prop%c_sat(ii))) &
               / (1D0 + dt_cond*sum(CR(:,ii)))


    ! apc scheme here. NOTE that for sulfuric acid vapour_prop%c_sat(ii) = 0 so the last term will vanish
    conc_pp(:,ii) = conc_pp_old(:,ii) + dt_cond*CR(:,ii)*(MIN(conc_guess,conc_tot(ii)) &
                  - kohler_effect(:,ii)*vapour_prop%c_sat(ii))

    ! Prevent overestimation of evaporation by setting negative concentrations to zero
    WHERE ( conc_pp(:,ii)<0D0 ) conc_pp(:,ii) = 0D0

    ! Prevents particles to grow over the saturation limit. For organics, no acids or HOA included
    IF (ii <= n_cond_org-1) then
      WHERE (diameter>0d0 .and. conc_pp(:,ii)>conc_pp_eq(:,ii) .AND. conc_vap_old(ii)<(Kelvin_Effect(:,ii) * vapour_prop%c_sat(ii))) conc_pp(:,ii) = conc_pp_eq(:,ii)
    END IF

    ! Update conc_vap: total conc - particle phase concentration
    conc_vap(ii) = conc_tot(ii) - sum(conc_pp(:,ii))


  END DO

  ! Calculate dmass for PSD
  DO ii = 1, n_cond_tot
    dmass(:,ii) = (conc_pp(:,ii) - conc_pp_old(:,ii))*vapour_prop%molar_mass(ii) / Na
  END DO
  ! dmass(:,n_cond_tot-1) = max(dmass(:,n_cond_tot-1), 0)
  do ii = 1, n_bins_particle
    if (n_conc(ii)>0d0) THEN
      dmass(ii,:) = dmass(ii,:) / n_conc(ii)
    else
      dmass(ii,:) = 0d0
    END if
  END DO

  ! derive diameter changes for integration time-step optimization
  DO ii = 1, n_bins_particle
    IF (SUM(conc_pp_old(ii,:)) > 0.d0 .and. n_conc(ii) > 1.d-10) THEN
      d_dpar(ii) = (SUM(conc_pp_old(ii,:)*vapour_prop%molar_mass(:)) / Na /n_conc(ii) + SUM(dmass(ii,:))) /  (SUM(conc_pp_old(ii,:)*vapour_prop%molar_mass(:)) / Na /n_conc(ii))
      d_dpar(ii) = d_dpar(ii) ** (1.d0/3.d0) - 1.d0
      !PRINT*, 'i,mass, mass change, d_dp', ii, SUM(conc_pp_old(ii,:)*vapour_prop%molar_mass(:)) / Na /n_conc(ii), SUM(dmass(ii,:)), d_dpar(ii)
    end if
  END DO

  ! derive the relative changes in the vapor phase
  DO ii = 1, n_cond_tot
    IF (conc_vap(ii) > 1.d5) THEN
      d_vap(ii) = (conc_vap_old(ii) - conc_vap(ii)) / conc_vap(ii)
      !PRINT*, 'ii, conc, dconc, d_vap', ii, conc_vap(ii), conc_vap_old(ii), conc_vap(ii) - conc_vap_old(ii), d_vap(ii)
    END IF
  END DO

END SUBROUTINE Condensation_apc




SUBROUTINE Coagulation_routine(dconc_coag, dt_coag, d_npar) ! Add more variables if you need it
  use omp_lib
  REAL(dp), DIMENSION(n_bins_particle) :: n_conc,diameter, mass
  ! type(PSD), INTENT(IN) :: particles
  ! REAL(dp), INTENT(IN) :: timestep
  REAL(dp), DIMENSION(n_bins_particle) :: volume
  ! REAL(dp) :: dp_max
  ! REAL(dp) :: temperature, pressure
  integer :: i,j,m,ii
  REAL(dp), DIMENSION(n_bins_particle,n_bins_particle) :: coagulation_coef        ! coagulation coefficients [m^3/s]
  REAL(dp), DIMENSION(n_bins_particle,n_bins_particle), intent(inout) :: dconc_coag    ! coagulation coefficients [m^3/s]
  REAL(dp), DIMENSION(n_bins_particle) :: slip_correction,diffusivity,dist, speed_p,free_path_p
  REAL(dp), DIMENSION(n_bins_particle, n_bins_particle) :: Beta_Fuchs
  ! REAL(dp), DIMENSION(n_bins_particle+1) :: Vp
  REAL(dp) :: dyn_visc  ! dynamic viscosity, kg/(m*s)
  REAL(dp) :: l_gas     ! Gas mean free path in air
  REAL(dp) :: a
  REAL(dp), INTENT(IN) :: dt_coag ! integration timestep for coagulation [s]
  REAL(dp), INTENT(INOUT) :: d_npar(:) ! relative change in particle number concentrations [1]

  REAL(dp) :: stt
  integer:: omp_rank
  ! dconc_coag = 0d0


  n_conc = get_conc()
  diameter = get_dp()
  volume = get_volume()
  mass = get_mass()
  ! print*, 'mass from volume', mass
  ! print*, 'mass from compos',sum(get_composition(), 2)
! print FMT_TIME, GTIME%hms
! print*, sum(n_conc)
! print*, sum(diameter)
! print*, sum(volume)
! print*, sum(mass)

  ! dp_max = diameter(n_bins_particle)*diameter(n_bins_particle)/diameter(n_bins_particle-1)
  ! Vp(1:n_bins_particle) = volume(1:n_bins_particle)
  ! Vp(n_bins_particle+1) = (pi*dp_max**3.)/6.

  ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603

  ! Dynamic viscosity of air
  dyn_visc = 1.8D-5*(GTEMPK/298.0d0)**0.85
  ! Gas mean free path in air (m)
  l_gas=2D0*dyn_visc/(GPRES*SQRT(8D0*Mair/(pi*R*GTEMPK)))
  ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)
  slip_correction = 1D0+(2D0*l_gas/(diameter)) * (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))
  ! Diffusivity for the different particle sizes m^2/s
  diffusivity = slip_correction*kb*GTEMPK/(3D0*pi*dyn_visc*diameter)

  ! Speed of particles (m/s)
  speed_p = SQRT(8D0*kb*GTEMPK/(pi*mass))
! print*, 'speed 55', SQRT(8D0*kb*GTEMPK/(pi*mass(55)))
  ! Particle mean free path (m)
  free_path_p = 8D0*diffusivity/(pi*speed_p)
! print*, 'fp 55', 8D0*diffusivity(55)/(pi*speed_p(55))
  ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)
  dist = (1D0/(3D0*diameter*free_path_p))*((diameter + free_path_p)**3D0 &
  - (diameter**2D0 + free_path_p**2D0)**(3D0/2D0)) - diameter

  DO i = 1,n_bins_particle
    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
     Beta_Fuchs(i,:) = 1D0/((diameter + diameter(i))/(diameter + diameter(i) +&
     2D0*(dist**2D0 + dist(i)**2D0)**0.5D0) + 8D0*(diffusivity + diffusivity(i))/&
     (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter + diameter(i))))

    ! coagulation rates between two particles of all size combinations  (m^3/s)
    coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs(i,:)*(diameter*diffusivity(i) &
                          + diameter*diffusivity + diameter(i)*diffusivity &
                          + diameter(i)*diffusivity(i))
  END DO
  ! print*, 'SB', coagulation_coef(55, :)
  !
  ! do i=1,n_bins_particle
  !     do j=i,n_bins_particle
  !            coagulation_coef(i,j)=calc_coag_coeff(i,j,diameter/2d0,mass,GTEMPK,GPRES)
  !     end do
  !     do j=1,i
  !         coagulation_coef(i,j)=coagulation_coef(j,i)
  !     end do
  ! end do
  ! print*, 'UH', coagulation_coef(55, :)
dconc_coag = 0d0
  if (USE_OPENMP) THEN
    call omp_set_num_threads(4)
    !$OMP PARALLEL shared(n_conc, coagulation_coef) private(omp_rank, stt,j,m)
    stt=omp_get_wtime()
    omp_rank = omp_get_thread_num()

    !$OMP PARALLEL DO
    do j =1, n_bins_particle
      do m = 1, j
        if (m==j) then
          a=0.5_dp
        else
          a= 1.0_dp
        END if
        ! print*, 'b',dconc_coag(j,m)
        ! dconc_coag(j,m) = 0d0 ! 1/m^3s
        dconc_coag(j,m) = a * coagulation_coef(m,j) * n_conc(j) * n_conc(m) ! 1/m^3s
        ! print*, 'a',dconc_coag(j,m)
      END DO
    END DO
    !$OMP END parallel DO
    !$omp END parallel


  ! If no OPENMP. NOTE the actual loops are identical
  ELSE
    do j =1, n_bins_particle
      do m = 1, j
        if (m==j) then
          a=0.5_dp
        else
          a= 1.0_dp
        END if

        ! if (GTIME%sec>3300) print*, 'bc',dconc_coag(j,m)
        dconc_coag(j,m) = a * coagulation_coef(m,j) * n_conc(j) * n_conc(m) ! 1/m^3s
        ! if (GTIME%sec>3300) print*, 'ac',dconc_coag(j,m)
      END DO
    END DO
  END IF
if (GTIME%sec>3300) THEN

end if
! Convert dconc_coag to units of timestep
dconc_coag = dconc_coag * dt_coag

!Check whether changes are within limits:
DO ii = 1, n_bins_particle
  IF (n_conc(ii) > 1.d0 .and. sum(dconc_coag(ii,ii:)) > 1.d-20) THEN
    d_npar(ii) = sum(dconc_coag(ii, ii:)) / n_conc(ii)
  ELSE
    d_npar(ii) = 0.d0
  END IF
END DO

END SUBROUTINE Coagulation_routine

!
!
! ! Uhma coag coeff
! function calc_coag_coeff(ii,jj,radius,mass,temp,pres)
!
! implicit none
!
! integer, intent(in) :: ii, jj
! real(dp),dimension(:),intent(in) :: radius,mass
! real(dp),intent(in) :: temp,pres
! real(dp) :: Pstand    = 1.01325e5
! real(dp) :: calc_coag_coeff
! real(dp) :: free_path, viscosity, vel_12, rad_12, dif_12
! real(dp) :: prov1, prov2, dist, continuum, free_molec
! real(dp), dimension(2) :: radii, m_part, knudsen, corr, diff, veloc, omega
! real(dp), parameter :: stick_prob=1. !sticking probability
!
! ! air mean free path (m) and viscosity (kg/ms) (where from?)
! free_path = (6.73d-8*temp*(1.+110.4/temp))/(296.*pres/pstand*1.373)
! viscosity = (1.832d-5*406.4*temp**1.5)/(5093*(temp+110.4))
!
! ! for both sections
! radii(1) = radius(ii)        ! radii of colliding particles (m)
! radii(2) = radius(jj)
! m_part(1) = mass(ii)        ! masses of colliding particles (kg)
! m_part(2) = mass(jj)
!
! knudsen = free_path/radii                                                                ! particle Knudsen number
! corr = 1. + knudsen*(1.142+0.558*exp(-0.999/knudsen))        ! Cunninghan correction factor (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
! diff = kb*temp*corr/(6*pi*viscosity*radii)                ! particle diffusion coefficient (m^2/s)
! veloc = sqrt((8.*kb*temp)/(pi*m_part))                        ! mean thermal velocity of a particle (m/s)
! omega = 8.*diff/(pi*veloc)             ! mean free path (m)
!
! vel_12 = sqrt(veloc(1)**2 + veloc(2)**2)        ! mean relative thermal velocity
! rad_12 = sum(radii)
! dif_12 = diff(1)+diff(2)               ! relative diffusion coefficient
! continuum = 4.*pi*rad_12*dif_12        ! flux in continuum regime
! free_molec = pi*vel_12*rad_12**2       ! flux in free molecular regime
!
! ! flux matching according to Fuchs (1964) (e.g. Seinfeld & Pandis p. 661)
! prov1 = (rad_12+omega(1))**3 - (rad_12**2 + omega(1)**2)**1.5
! prov1 = prov1/(3.*rad_12*omega(1)) - rad_12
! prov2 = (rad_12+omega(2))**3 - (rad_12**2 + omega(2)**2)**1.5
! prov2 = prov2/(3.*rad_12*omega(2)) - rad_12
! dist = sqrt(prov1**2 + prov2**2)        ! distance at which fluxes are matched
!
! ! coagulation coefficient between particles [m^3/s]
! calc_coag_coeff = stick_prob*continuum / (rad_12/(rad_12+dist) + continuum/free_molec)
!
! end function calc_coag_coeff
!
!

function collision_rate(jj,diameter, mass,vapour_prop)

    implicit none

    type(vapour_ambient), INTENT(IN) :: vapour_prop

    integer, INTENT(IN) :: jj                      ! the 'number' of condensing vapour
    REAL(dp),dimension(n_bins_particle) :: collision_rate
    REAL(dp), INTENT(IN) :: diameter(:), mass(:)
    ! integer :: ii
    REAL(dp) :: Diff_org

    REAL(dp) :: air_free_path,viscosity !dif_vap, r_vap, n1, , r_hyd, r_h2o,zero
    REAL(dp), DIMENSION(n_bins_particle) :: knudsen, Diff_par, slip_correction


    REAL(dp):: Diff_H2SO4_0    ! gas diffusivity for H2SO4 AT RH=0%

    REAL(dp) :: Diff_H2SO4,Keq1,Keq2,d_H2SO4,mH2SO4,speedH2SO4,&
                dorg,dens_air,speedorg

    REAL(dp),dimension(n_bins_particle) :: gasmeanfpH2SO4,speed_p,KnH2SO4,f_corH2SO4,DH2SO4eff,gasmeanfporg,&
                                        Knorg,f_cororg,Dorgeff

    ! viscosity of air, density oif air and mean free path in air
    viscosity     = 1.8D-5*(GTEMPK/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air      = Mair*GPRES/(R*GTEMPK)
    air_free_path = 2D0*viscosity/(GPRES*SQRT(8D0*Mair/(pi*R*GTEMPK))) ! gas mean free path in air
    ! knudsen number
    knudsen = 2D0 * air_free_path/diameter

    ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34 pg 407)
    slip_correction = 1D0 + knudsen * &
                      (1.257D0 + 0.4D0*exp(-1.1D0/(2D0*air_free_path/diameter)))

    ! particle diffusion coefficient (m^2 s^-1)
    Diff_par = (kb * GTEMPK * slip_correction) / (3D0 * pi * viscosity * diameter)

    speed_p = 0d0
    where (mass> 0)
      speed_p = SQRT(8D0*kb*GTEMPK/(pi*mass)) ! speed of particle unit(mass)=kg
    END where

   ! write(*,*) 'Carlton debug check molecvol line 276', molecvol(jj)
    IF (vapour_prop%cond_type(jj)==2) THEN
     ! For testing. In the main model we use RH dependence
     ! H2SO4 mass transfer rate (m3/s) ,RH dependent Diffusion, diameter and mass

     ! diffusivity of H2SO4 at RH=0%
      Diff_H2SO4_0 = 0.09D-4

      Keq1=0.13D0      ! Hanson and Eisele
      Keq2=0.016D0     ! Hanson and Eisele

    ! Diffusivity H2SO4 at ambient RH
       Diff_H2SO4 =  (Diff_H2SO4_0 + 0.85D0*Diff_H2SO4_0*Keq1*GRH + 0.76D0*Diff_H2SO4_0*Keq1*Keq2*GRH**2D0) &
        /(1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)

       ! RH dependent H2SO4 diameter
       d_H2SO4 = (((vapour_prop%molar_mass(jj)/vapour_prop%density(jj)/Na)*6D0/pi)**(1D0/3D0) + &
          Keq1*GRH*(((vapour_prop%molar_mass(jj) + 18D-3)/vapour_prop%density(jj)/Na)*6D0/pi)**(1D0/3D0) + &
          Keq1*Keq2*GRH**2D0*(((vapour_prop%molar_mass(jj) + 36D-3)/vapour_prop%density(jj)/Na)*6D0/pi)**(1D0/3D0)) / &
          (1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)
      !

      !  ! RH H2SO4 molecular mass
       mH2SO4 = (vapour_prop%molar_mass(jj)/Na + Keq1*GRH*(vapour_prop%molar_mass(jj) + 18D-3)/Na +  &
                Keq1*Keq2*GRH**2D0*(vapour_prop%molar_mass(jj) + 36D-3)/Na) / &
                (1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)

       ! speed of H2SO4 molecule
       speedH2SO4 = SQRT(8D0*kb*GTEMPK/(pi*mH2SO4))
       ! gasmeanfpH2SO4 = 0

        ! collision of H2SO4 molecule with particle. Only for bins where mass > 0
        gasmeanfpH2SO4 = 3D0*(Diff_H2SO4 + Diff_par)/SQRT(speedH2SO4**2D0 + speed_p**2D0)

       ! Knudsen number H2SO4
       KnH2SO4=2D0*gasmeanfpH2SO4/(diameter + d_H2SO4)

       ! Fuchs-Sutugin correction factor for transit
       f_corH2SO4 = (0.75*vapour_prop%alpha(jj)*(1D0 + KnH2SO4))/ &
                    (KnH2SO4**2D0 + KnH2SO4 + 0.283*KnH2SO4*vapour_prop%alpha(jj) + 0.75*vapour_prop%alpha(jj))

       DH2SO4eff = (Diff_H2SO4 + Diff_par)*f_corH2SO4    !m2/s

       collision_rate = 2D0*pi*(diameter + d_H2SO4)*DH2SO4eff  !mass transfer rate m3/s

       ! write(*,*) 'mH2SO4', mH2SO4

    ! Knudsen number for organic condensable compounds
  ELSE ! Organics mass transfer rate (m3/s)            print*, 'xxx1'


       ! diameter of organic compounds
       dorg= (6D0*vapour_prop%molec_volume(jj)/pi)**(1D0/3D0)               !estimated diameter (m)

       ! diffusivity organic compound
       if (Use_atoms) THEN
         Diff_org  = 1D-7 * GTempK**1.75D0 * SQRT( 1/(Mair*1D3) + 1/(1d3*VAPOUR_PROP%molar_mass(jj)) ) &
                   / (GPres/1.01325D5 * (Vol_org(jj)**(1D0/3D0) + 20.1**(1D0/3D0))**2D0)

       ELSE
         Diff_org=5D0/(16D0*Na*dorg**2D0*dens_air)*&
         sqrt(R*GTEMPK*Mair/(2D0*pi)*((vapour_prop%molar_mass(jj) + Mair)/vapour_prop%molar_mass(jj)))
       END IF

       speedorg=SQRT(8D0*Kb*GTEMPK/(pi*vapour_prop%molec_mass(jj))) !speed of organic molecules

       gasmeanfporg=3D0*(Diff_org + Diff_par)/SQRT(speedorg**2D0 + speed_p**2D0)

       ! Knudsen number organic comp
       Knorg=2D0*gasmeanfporg/(diameter + dorg)

       ! Fuchs-Sutugin correction factor for transit
       f_cororg=(0.75*vapour_prop%alpha(jj)*(1D0 + Knorg))/&
                 (Knorg**2D0 + Knorg+0.283*Knorg*vapour_prop%alpha(jj) + 0.75*vapour_prop%alpha(jj))

       Dorgeff = (Diff_org + Diff_par)*f_cororg                    ! m^2/s

       collision_rate = 2D0*pi*(diameter + dorg)*Dorgeff        ! mass transfer coefficient s^-1

     ENDIF

END function collision_rate



function collision_rate_uhma(jj,diameter,mass,vapour_prop)

    implicit none

    type(vapour_ambient),intent(in) :: vapour_prop
    integer, intent(in) :: jj                                                        ! the 'number' of condensing vapour
    real(dp),dimension(:),intent(in) :: diameter,mass
    real(dp), dimension(n_bins_particle) :: collision_rate_uhma, radius

    real(dp) :: air_free_path, r_vap, viscosity, r_h2o,temp,rh,pres
    real(dp), dimension(n_bins_particle) :: knudsen, corr, dif_part
    real(dp),dimension(n_cond_org+1) :: molecmass,molecvol,molarmass,molarmass1,alpha,dens
    !variables for new method(when input_flag=0)
    real(dp) :: DH2SO40,DH2SO4,Keq1,Keq2,d_H2SO4,mH2SO4,speedH2SO4,Dorg,dX,dens_air,speedorg
    real(dp),dimension(n_bins_particle) :: gasmeanfpH2SO4,speed_p,KnH2SO4,f_corH2SO4,DH2SO4eff,gasmeanfporg,&
                                        Knorg,f_cororg,Dorgeff
    radius = diameter/2d0
    temp=GTEMPK
    rh=grh ! RH unit %,values from 1-100
    pres=gpres
    molecvol=vapour_prop%molec_volume !m3
    molarmass=vapour_prop%molar_mass*1d3 !g/mol
    molarmass1=molarmass*1d-3
    molecmass=vapour_prop%molec_mass !kg

    alpha=vapour_prop%alpha
    dens=vapour_prop%density
   ! air mean free path (m) and viscosity (kg s^-1 m^-1) (where from?)
    air_free_path = (6.73d-8*temp*(1.+110.4/temp))/(296.*pres/101325d0*1.373)
    viscosity = (1.832d-5*406.4*temp**1.5)/(5093*(temp+110.4))
    r_h2o = (3*(18 * 1.d-3 / Na)/1000d0/(4*pi))**(1./3.)        ! radius of water molecule
    r_vap = (molecvol(jj)*3./(4.*pi))**(1./3.)        ! radius of condensable molecules (m)

    ! for particles in each size section
    knudsen = air_free_path/radius                                                        ! particle Knudsen number
    corr = 1. + knudsen*(1.142+0.558*exp(-0.999/knudsen))        ! Cunninghan correction factor (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
    dif_part = kb*temp*corr/(6*pi*viscosity*radius)        ! particle diffusion coefficient (m^2/s)

   IF (.TRUE.) THEN
     ! H2SO4 mass transfer rate (m3/s) ,RH dependent Diffusion, diameter and mass

     IF (jj .EQ. VAPOUR_PROP%vbs_bins )  THEN

       DH2SO40= 0.09D-4 ! gas diffusivity H2SO4 m2/s RH=0%
       Keq1=0.13D0
       Keq2=0.016D0
       DH2SO4 = (DH2SO40+0.85D0*DH2SO40*Keq1*rh+0.76D0*DH2SO40*Keq1*Keq2*rh**2D0)&
          /(1D0+Keq1*rh+Keq1*Keq2*rh**2D0) ! Diffusivity H2SO4 at ambient RH

       d_H2SO4=(((molarmass1(jj)/dens(jj)/Na)*6D0/pi)**(1D0/3D0)+&
          Keq1*rh*(((molarmass1(jj)+18D-3)/dens(jj)/Na)*6D0/pi)**(1D0/3D0)+&
          Keq1*Keq2*rh**2D0*(((molarmass1(jj)+36D-3)/dens(jj)/Na)*6D0/pi)**(1D0/3D0))/&
          (1D0+Keq1*rh+Keq1*Keq2*rh**2D0) ! RH dependent H2SO4 diameter

       mH2SO4=(molarmass1(jj)/Na+Keq1*rh*(molarmass1(jj)+18D-3)/Na+Keq1*Keq2*rh**2D0*(molarmass1(jj)+36D-3)/Na)/&
          (1D0+Keq1*rh+Keq1*Keq2*rh**2D0) ! RH H2SO4 molecular mass

       speedH2SO4=SQRT(8D0*kb*temp/(pi*mH2SO4)) ! speed of H2SO4 molecule
       gasmeanfpH2SO4=0
       where (mass>0)
        speed_p=SQRT(8D0*kb*temp/(pi*mass)) ! speed of particle unit(mass)=kg
        gasmeanfpH2SO4=3D0*(DH2SO4+dif_part)/SQRT(speedH2SO4**2D0+speed_p**2D0)
       endwhere



       KnH2SO4=2D0*gasmeanfpH2SO4/(2*radius+d_H2SO4) ! Knudsen number H2SO4

       f_corH2SO4=(0.75*alpha(jj)*(1D0+KnH2SO4))/(KnH2SO4**2D0+KnH2SO4+0.283*KnH2SO4*alpha(jj)+0.75*alpha(jj)) ! Fuchs-Sutugin correction factor for transit
       DH2SO4eff=(DH2SO4+dif_part)*f_corH2SO4    !m2/s
     !  where (mass>0)
         collision_rate_uhma=2D0*pi*(2.*radius+d_H2SO4)*DH2SO4eff  !mass transfer rate m3/s
     !  elsewhere
     !    collision_rate_uhma=0
     !  endwhere

     ELSE
      !............................................................................
      !Organics mass transfer rate (m3/s)

       dX= (6.*molecvol(jj)/pi)**(1./3.)    !estimated diameter ( m)
       dens_air=mair*1e-3* pres/(R*temp)     ! density of air  kg/m3
       Dorg=5./(16.*Na*dX**2.*dens_air)*&
          sqrt(R*temp*mair*1e-3/(2.*pi)*((molarmass1(jj)+mair*1e-3)/molarmass1(jj))) ! diffusivity organic compound

       speedorg=SQRT(8D0*kb*temp/(pi*molecmass(jj))) !spped of organic molecules
       gasmeanfporg=0
       where (mass>0)
        speed_p=SQRT(8D0*kb*temp/(pi*mass)) ! speed of particle unit(mass)=kg
        gasmeanfporg=3D0*(Dorg+dif_part)/SQRT(speedorg**2D0+speed_p**2D0)
       endwhere
       ! print*, 'UH', dif_part

       Knorg=2D0*gasmeanfporg/(2*radius+dX)       ! Knudsen number organic comp
       f_cororg=(0.75*alpha(jj)*(1D0+Knorg))/&
         (Knorg**2D0+Knorg+0.283*Knorg*alpha(jj)+0.75*alpha(jj))     ! Fuchs-Sutugin correction factor for transit
       Dorgeff=(Dorg+dif_part)*f_cororg                    ! m^2/s
      ! where (mass>0)
         collision_rate_uhma=2D0*pi*(2.*radius+dX)*Dorgeff        ! mass transfer coefficient s^-1
      ! elsewhere
      !   collision_rate_uhma=0
      ! endwhere
!         write(*,*) 'haha',collision_rate_uhma(1),collision_rate_uhma(30)
!         stop

     ENDIF



   ENDIF
end function collision_rate_uhma




END MODULE aerosol_dynamics
