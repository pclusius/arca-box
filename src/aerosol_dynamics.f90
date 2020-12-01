MODULE aerosol_dynamics

USE second_precision, ONLY: dp
USE constants
USE INPUT
USE PSD_scheme
USE auxillaries

IMPLICIT NONE

CONTAINS


! ======================================================================================================================
! Main condensation subroutine. We use the Analytical predictor of condensation scheme Jacobson 1997c, 2002 and
! fundamentals of atmospheric modelling. Dmass, which is the flux onto or removed from the particles is the outcome
! of this subroutine which is fed to PSD
! ======================================================================================================================
SUBROUTINE Condensation_apc(VAPOUR_PROP, conc_vap, dmass, dt_cond, d_dpar,d_vap)
    IMPLICIT NONE
    type(vapour_ambient), INTENT(IN)            :: VAPOUR_PROP  ! Properties of condensing vapours
    REAL(dp), INTENT(INOUT)                     :: conc_vap(:)  ! [#/m^3], condensing vapour concentrations, DIM(n_cond_tot)
    REAL(dp), INTENT(INOUT)                     :: dmass(:,:)   ! [kg/m^3] change of mass per particle in particle phase, DIM(n_bins_par, n_cond_tot)
    REAL(dp), INTENT(IN)                        :: dt_cond      ! integration timestep for condensation [s]
    REAL(dp), INTENT(INOUT)                     :: d_dpar(:)    ! relative change in particle diameter [1]
    REAL(dp), INTENT(INOUT)                     :: d_vap(:)     ! relative change in vapor concentrations [1]
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: CR           ! [/s] collisions per second, collision rate * concentration
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: Kelvin_Effect
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: kohler_effect
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: conc_pp_old
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: xorg         ! mole fraction of each organic compound in each bin
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: conc_pp_eq
    REAL(dp), DIMENSION(n_cond_tot)             :: conc_tot
    REAL(dp), DIMENSION(n_cond_tot)             :: conc_vap_old
    REAL(dp)                                    :: n_conc(n_bins_par)
    REAL(dp)                                    :: diameter(n_bins_par)
    REAL(dp)                                    :: volume(n_bins_par)
    REAL(dp)                                    :: mass(n_bins_par) ! [kg/m^3]
    REAL(dp)                                    :: conc_pp(n_bins_par,n_cond_tot) ! [#/m^3] Particle phase concentrations, DIM(n_bins_par,n_cond_tot)
    REAL(dp)                                    :: st           ! Surface tension for all compounds
    REAL(dp)                                    :: conc_guess
    REAL(dp)                                    :: sum_org      ! Sum of organic concentration in bin, transient variable

  INTEGER :: ii
  if (TEMP_DEP_SURFACE_TENSION) THEN
      st = (1D0/3D0) * ((76.1d0 - 0.155d0*(GTEMPK-273.15d0))*1d-3)
  ELSE
    st = VAPOUR_PROP%surf_tension(1)
  END IF
  conc_pp = 0d0
  conc_pp = get_composition()
  n_conc = get_conc()
  diameter = get_dp()
  volume = get_volume()
  mass = sum(conc_pp, 2)
  conc_pp_eq = 0d0
  ! Fill the number concentration composition matrix
  do ii=1, n_bins_par
    conc_pp(ii,:) = conc_pp(ii,:) * Na / VAPOUR_PROP%molar_mass * n_conc(ii)
  end do

  ! original gas phase concentration
  conc_vap_old = conc_vap
  ! original particle phase concentration
  conc_pp_old = conc_pp
  ! vapour phase + particle phase
  conc_tot = conc_vap + sum(conc_pp,1)

  xorg = 1d0
  if (use_raoult) THEN
    DO ii=1,n_bins_par
      if (volume(ii) >0.0 .and. diameter(ii) > 1D-9) then
        sum_org = sum(conc_pp(ii,:),DIM=1, Mask=(VAPOUR_PROP%cond_type==1))
        if (sum_org>0) THEN
          xorg(ii,:)= conc_pp(ii,:)/sum_org
        END if
      END if
    END DO
  END if

  ! Kelvin and Kohler factors.
  DO ii = 1, n_cond_tot
    ! Kelvin factor takes into account the curvature of particles. Unitless
    Kelvin_Effect(:,ii) = 1d0

    Kelvin_Effect(:,ii) = 1D0 + 2D0*st*VAPOUR_PROP%molar_mass(ii) &
                        / (Rg*GTEMPK*VAPOUR_PROP%density(ii)*diameter/2D0)

    ! Kohler factor the solute partial pressure effect. Unitless
    kohler_effect(:,ii) = Kelvin_Effect(:,ii)*xorg(:,ii)
  END DO
  ! Treat H2SO4 specially
  kohler_effect(:,n_cond_tot) = Kelvin_Effect(:,n_cond_tot)

  ! Approximate equilibrium concentration (#/m^3) of each compound in each size bin
  DO ii=1,n_bins_par
    conc_pp_eq(ii,1:n_cond_tot-1) = conc_vap_old(1:n_cond_tot-1)*SUM(conc_pp_old(ii,1:n_cond_tot-1)) &
                                    / (Kelvin_Effect(ii,1:n_cond_tot-1)*VAPOUR_PROP%c_sat(1:n_cond_tot-1))
  END DO

  DO ii = 1, n_cond_tot

      ! Total number of collisions/s
    CR(:,ii) = n_conc * collision_rate(ii,diameter,mass,VAPOUR_PROP)

    ! Update sulfuric acid condensation sink
    if (ii == VAPOUR_PROP%ind_H2SO4) GCS = sum(CR(:,ii))

    ! apc scheme here. NOTE that for sulfuric kohler_effect = Kelvin_Effect
    conc_guess = (conc_vap(ii) + dt_cond*sum(CR(:,ii)*kohler_effect(:,ii)*VAPOUR_PROP%c_sat(ii))) &
               / (1D0 + dt_cond*sum(CR(:,ii)))

    ! apc scheme here. NOTE that for sulfuric acid VAPOUR_PROP%c_sat(ii) = 0 so the last term will vanish
    conc_pp(:,ii) = conc_pp_old(:,ii) + dt_cond*CR(:,ii)*(MIN(conc_guess,conc_tot(ii)) &
                  - kohler_effect(:,ii)*VAPOUR_PROP%c_sat(ii))

    ! Prevent overestimation of evaporation by setting negative concentrations to zero
    WHERE ( conc_pp(:,ii)<0D0 ) conc_pp(:,ii) = 0D0

    ! Prevents particles to grow over the saturation limit. For organics, no acids or GENERIC included
    IF (ii < VAPOUR_PROP%ind_GENERIC) then
      WHERE (diameter>0d0 .and. conc_pp(:,ii)>conc_pp_eq(:,ii) &
                .AND. conc_vap_old(ii)<(Kelvin_Effect(:,ii) * VAPOUR_PROP%c_sat(ii))) &
                conc_pp(:,ii) = conc_pp_eq(:,ii)
    END IF

    ! XXX NOTE We change GENERIC back to what it was before
    conc_pp(:,VAPOUR_PROP%ind_GENERIC) = conc_pp_old(:,VAPOUR_PROP%ind_GENERIC)

    ! Update conc_vap: total conc - particle phase concentration
    conc_vap(ii) = conc_tot(ii) - sum(conc_pp(:,ii))
    ! if (ii==VAPOUR_PROP%ind_GENERIC) print*, 'vapor conc after', conc_vap(ii)


  END DO

  ! Calculate dmass for PSD
  DO ii = 1, n_cond_tot
    dmass(:,ii) = (conc_pp(:,ii) - conc_pp_old(:,ii))*VAPOUR_PROP%molar_mass(ii) / Na
  END DO
  ! dmass(:,n_cond_tot-1) = max(dmass(:,n_cond_tot-1), 0)
  do ii = 1, n_bins_par
    if (n_conc(ii)>1d-10) THEN
      dmass(ii,:) = dmass(ii,:) / n_conc(ii)
    else
      dmass(ii,:) = 0d0
    END if
  END DO

  ! derive diameter changes for integration time-step optimization
  DO ii = 1, n_bins_par
    IF (SUM(conc_pp_old(ii,:)) > 0.d0 .and. n_conc(ii) > 1.d-10) THEN
      d_dpar(ii) = (SUM(conc_pp_old(ii,:)*VAPOUR_PROP%molar_mass(:)) / Na /n_conc(ii) + SUM(dmass(ii,:))) &
                    / (SUM(conc_pp_old(ii,:)*VAPOUR_PROP%molar_mass(:)) / Na /n_conc(ii))
      if (d_dpar(ii)<0) THEN
          d_dpar(ii) = -1d0*((-1d0*d_dpar(ii)) ** (1.d0/3.d0)) - 1.d0
      ELSE
          d_dpar(ii) = (d_dpar(ii) ** (1.d0/3.d0)) - 1.d0
      END IF

      !PRINT*, 'i,mass, mass change, d_dp', ii, SUM(conc_pp_old(ii,:)*VAPOUR_PROP%molar_mass(:)) / Na /n_conc(ii), SUM(dmass(ii,:)), d_dpar(ii)
    ELSE
      d_dpar(ii) = 0.d0
    END IF
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
    IMPLICIT NONE
    REAL(dp), DIMENSION(n_bins_par,n_bins_par), intent(inout) :: dconc_coag    ! coagulation coefficients [m^3/s]
    REAL(dp), INTENT(IN)                        :: dt_coag ! integration timestep for coagulation [s]
    REAL(dp), INTENT(INOUT)                     :: d_npar(:) ! relative change in particle number concentrations [1]
    REAL(dp), DIMENSION(n_bins_par)             :: n_conc,diameter, mass
    REAL(dp), DIMENSION(n_bins_par)             :: volume
    integer                                     :: i,j,m,ii
    REAL(dp), DIMENSION(n_bins_par,n_bins_par)  :: coagulation_coef        ! coagulation coefficients [m^3/s]
    REAL(dp), DIMENSION(n_bins_par,n_bins_par)  :: sticking_prob           ! alpha in Fuchs' beta, S&P (2016) p. 550, eq. 13.56
    REAL(dp), DIMENSION(n_bins_par)             :: slip_correction
    REAL(dp), DIMENSION(n_bins_par)             :: diffusivity,dist
    REAL(dp), DIMENSION(n_bins_par)             :: speed_p
    REAL(dp), DIMENSION(n_bins_par)             :: free_path_p
    REAL(dp), DIMENSION(n_bins_par, n_bins_par) :: Beta_Fuchs
    REAL(dp)                                    :: dyn_visc  ! dynamic viscosity, kg/(m*s)
    REAL(dp)                                    :: l_gas     ! Gas mean free path in air
    REAL(dp)                                    :: a, comp

    sticking_prob = 1
    n_conc        = get_conc()
    diameter      = get_dp()
    volume        = get_volume()
    mass          = get_mass()

  ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603

  ! Dynamic viscosity of air
  dyn_visc = 1.8D-5*(GTEMPK/298.0d0)**0.85
  ! Gas mean free path in air (m)
  l_gas=2D0*dyn_visc/(GPRES*SQRT(8D0*Mair/(pi*Rg*GTEMPK)))
  ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)
  slip_correction = 1D0+(2D0*l_gas/diameter) * (1.257D0 + 0.4D0 * exp(-1.1D0*diameter/(2D0*l_gas)))
  ! Diffusivity for the different particle sizes m^2/s
  diffusivity = slip_correction*kb*GTEMPK/(3D0*pi*dyn_visc*diameter)

  ! Speed of particles (m/s)
  speed_p = SQRT(8D0*kb*GTEMPK/(pi*mass))
! print*, 'speed 55', SQRT(8D0*kb*GTEMPK/(pi*mass(55)))
  ! Particle mean free path (m)
  free_path_p = 8D0*diffusivity/(pi*speed_p)
! print*, 'fp 55', 8D0*diffusivity(55)/(pi*speed_p(55))
  ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)
  dist = (SQRT(2D0)/(3D0*diameter*free_path_p))*((diameter + free_path_p)**3D0 &
  - (diameter**2D0 + free_path_p**2D0)**(3D0/2D0)) - diameter

  DO i = 1,n_bins_par
    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
     Beta_Fuchs(i,:) = 1D0/((diameter + diameter(i))/(diameter + diameter(i) &
     + 2D0*(dist**2D0 + dist(i)**2D0)**0.5D0) + (8D0/sticking_prob(i,:))*(diffusivity + diffusivity(i)) &
     /(((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter + diameter(i))))

    ! coagulation rates between two particles of all size combinations  (m^3/s)
    coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs(i,:)*(diameter*diffusivity(i) &
                          + diameter*diffusivity + diameter(i)*diffusivity &
                          + diameter(i)*diffusivity(i))
  END DO

dconc_coag = 0d0


do j=1, n_bins_par
    do m = j, n_bins_par
        if (m==j) then
            a=0.5_dp
        else
            a= 1.0_dp
        END if
        IF (n_conc(j) > 1.d0 .and. n_conc(m) > 1.d0) dconc_coag(j,m) = a * coagulation_coef(j,m) * n_conc(j) * n_conc(m) * dt_coag ! 1/m^3 '/ dt
        if ((dconc_coag(j,m)>n_conc(j)) .or. (dconc_coag(j,m)>n_conc(m))) THEN
            print FMT_WARN0, 'It seems that the upper size range is too small and should be increased (a good start is by 200%)'
            print FMT_WARN0, i2chr(j)//', '//i2chr(m)//', '//f2chr(coagulation_coef(j,m))//', '//f2chr(n_conc(j))//', '//f2chr(n_conc(m))
        END IF
    END DO

    ! Check whether changes are within limits:
    IF (n_conc(j) > 1.d0 .and. sum(dconc_coag(j,j:)) > 1.d-20) THEN
        d_npar(j) = sum(dconc_coag(j, j:)) / n_conc(j)
    ELSE
        d_npar(j) = 0.d0
    END IF
END DO



END SUBROUTINE Coagulation_routine


! Loss rate calculation by P. Roldin
SUBROUTINE deposition_velocity(d_p,ustar,Av,Au,Ad,V_chamber,T,p,E_field,dn_par, dt)
    IMPLICIT NONE
    REAL(dp), DIMENSION(n_bins_par), INTENT(in) :: d_p
    REAL(dp), INTENT(in)                        :: ustar,Av,Au,Ad,V_chamber,T,p,dt, E_field
    REAL(dp), DIMENSION(:), INTENT(inout)       :: dn_par ! s^-1
    REAL(dp), DIMENSION(n_bins_par)             :: D,vs,Sc,r,aa,bb,Ii,vdv,vdu,vdd,Cc,dens !,ve
    REAL(dp), DIMENSION(size(dn_par))           :: k_dep ! s^-1
    REAL(dp)                                    :: dyn_visc,kin_visc,l_gas,dens_air
    ! REAL(dp), DIMENSION(3) :: n
    ! INTEGER :: j
    if (E_field > 0) continue

    dens = (get_mass())/(pi*(get_dp())**3/6)
    ! Dry deposition model from Lai and Nazaroff J. Aerosol Sci., 31, 463–476, 2000. !%%%%%%%%
    dens_air=Mair*p/(Rg*T)		! Air density
    dyn_visc=1.8D-5*(T/298.)**0.85 ! dynamic viscosity
    kin_visc=dyn_visc/dens_air ! Pa s kinematic viscosity of air at T=273.15 K
    l_gas=2D0*dyn_visc/(p*SQRT(8D0*Mair/(pi*Rg*T))) ! Gas mean free path in air (m)
    Cc = 1D0+(2D0*l_gas/(d_p))*(1.257+0.4*exp(-1.1/(2D0*l_gas/d_p))) ! Cunninghams correction factor (seinfeld and Pandis eq 9.34
    D = Cc*kb*T/(3D0*pi*dyn_visc*d_p)              ! Diffusivitys for the different particle sizes m^2/s
    vs=dens*d_p**2D0*g_0*Cc/(18.*dyn_visc)           ! gravitational setting velocity
    Sc=kin_visc*D**(-1D0)
    r=d_p*ustar*(2D0*kin_visc)**(-1D0)


    aa=0.5*LOG((10.92*Sc**(-1D0/3D0)+4.3)**3D0/(Sc**(-1D0)+0.0609))+SQRT(3D0)&
    *ATAN((8.6-10.92*Sc**(-1D0/3D0))/(SQRT(3D0)*10.92*Sc**(-1D0/3D0)))

    bb=0.5*LOG((10.92*Sc**(-1D0/3D0)+r)**3/(Sc**(-1D0)+7.669D-4*r**3D0))&
    +SQRT(3D0)*ATAN((2D0*r-10.92*Sc**(-1D0/3D0))/(SQRT(3D0)*10.92*Sc**(-1D0/3D0)))
    Ii=3.64*Sc**(2D0/3D0)*(aa-bb)+39.

    vdv=ustar/Ii ! Deposition velocity vertical surface
    vdu=vs/(1D0-EXP(-vs*Ii/ustar)) ! Deposition velocity upward horizontal surface

    WHERE ((vs*Ii/ustar)<10d0)
        vdd=vs/(EXP(vs*Ii/ustar)-1D0) ! Deposition velocity downward horizontal surface
    ELSEWHERE
        vdd=vs*(EXP(-vs*Ii/ustar)) ! Deposition velocity downward horizontal surface
    END WHERE
    k_dep=(vdv*Av+vdu*Au+vdd*Ad)/V_chamber ! First order loss coefficient (s^-1) for non-charged particles, rectangular cavity

    ! Enhanced deposition velocity of charged particles (McMurry and Rader, 1985):

    ! n=(/1D0, 2D0, 3D0/) ! number of elemental charges (1,2,3)
    ! e=1.602D-19 ! Columbs
    ! DO j=1,3
    ! ve=n(j)*e*Cc*E_field/(3D0*pi*dyn_visc*d_p) ! Characteristic average deposition velocity due to electrstatic forces
    ! k_dep(j+1,:)=ve*(Av+Au+Ad)/V_chamber+k_dep(1,:) ! First order loss rate coefficient (s^-1) due to particle charge (1,2 or 3), rectangular cavity
    ! END DO
    if (GTIME%printnow) print*, 'G-mean k_dep', EXP(SUM(LOG(k_dep))/n_bins_par)
    dn_par = (get_conc()) - (get_conc())*EXP(-k_dep*dt) ! calculates how the particle number concentration changes in each size bin due to dry deposition

END SUBROUTINE deposition_velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE UPDATE_MOLECULAR_DIFF_AND_CSPEED(VAPOUR_PROP)
    IMPLICIT NONE
    type(vapour_ambient), INTENT(INOUT) :: VAPOUR_PROP  ! VapourProp file that gets updated
    INTEGER                             :: om,sa        ! short for indices

    sa = VAPOUR_PROP%ind_H2SO4
    om = VAPOUR_PROP%n_condorg

    ! --- Molecule thermal speed (m/s) ALL MOLECULES -------------------------------------------------------------------
    VAPOUR_PROP%c_speed = SQRT(8D0*Kb*GTEMPK/(pi*VAPOUR_PROP%molec_mass))

    ! --- ORGANIC Molecular diffusivity, 2 methods, where Fullers is the default ---------------------------------------
    IF (Use_atoms) THEN
        VAPOUR_PROP%diff(1:om)  = 1D-7 * GTempK**1.75D0 * SQRT( 1/(Mair*1D3) + 1/(1d3*VAPOUR_PROP%molar_mass(1:om)) )   &
                                / (GPres/1.01325D5 * (VAPOUR_PROP%diff_vol(1:om)**(1D0/3D0) + 20.1**(1D0/3D0))**2D0)
    ELSE
        VAPOUR_PROP%diff(1:om)  = 5d0/(16d0*Na* VAPOUR_PROP%molec_dia(1:om)**2d0*(Mair * GPRES/(Rg*GTEMPK)))            &
                                * sqrt(Rg*GTEMPK*Mair/(2.*pi)*((VAPOUR_PROP%molar_mass(1:om)+Mair)/VAPOUR_PROP%molar_mass(1:om)))
    END IF


    ! --- Molecular diffusivity, sulfuric acid diffusivity is RH dependant ---------------------------------------------

    VAPOUR_PROP%diff(sa)  = 1D-7 * GTempK**1.75D0 * SQRT( 1/(Mair*1D3) + 1/(1d3*VAPOUR_PROP%molar_mass(sa)) )   &
                            / (GPres/1.01325D5 * (VAPOUR_PROP%diff_vol(sa)**(1D0/3D0) + 20.1**(1D0/3D0))**2D0)

    IF (USE_RH_CORRECTION) THEN
        ! Diffusivity H2SO4 at ambient RH, Diffusivity of H2SO4 at RH=0%: 0.09D-4, defined in Constants.f90
        ! VAPOUR_PROP%diff(sa) = (VAPOUR_PROP%diff(sa) + 0.85D0*VAPOUR_PROP%diff(sa)*Keq1*GRH + 0.76D0*VAPOUR_PROP%diff(sa)*Keq1*Keq2*GRH**2D0) &
        !                         /(1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)
        VAPOUR_PROP%diff(sa) = VAPOUR_PROP%diff(sa) * (1d0 + 0.85D0*Keq1*GRH + 0.76D0*Keq1*Keq2*GRH**2D0) &
                                /(1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)
        ! RH dependent H2SO4 diameter
        VAPOUR_PROP%diff_dia(sa) = (((VAPOUR_PROP%molar_mass(sa)/VAPOUR_PROP%density(sa)/Na)*6D0/pi)**(1D0/3D0) &
                                    + Keq1*GRH*(((VAPOUR_PROP%molar_mass(sa) + 18D-3)/VAPOUR_PROP%density(sa)/Na)*6D0/pi)**(1D0/3D0) &
                                    + Keq1*Keq2*GRH**2D0*(((VAPOUR_PROP%molar_mass(sa) + 36D-3) &
                                    / VAPOUR_PROP%density(sa)/Na)*6D0/pi)**(1D0/3D0)) / (1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)

        VAPOUR_PROP%wet_dia(sa)  = VAPOUR_PROP%diff_dia(sa)
        ! RH dependent H2SO4 molecular mass
        VAPOUR_PROP%wet_mass(sa) = (VAPOUR_PROP%molar_mass(sa)/Na + Keq1*GRH*(VAPOUR_PROP%molar_mass(sa) + 18D-3)/Na &
                                    + Keq1*Keq2*GRH**2D0*(VAPOUR_PROP%molar_mass(sa) + 36D-3)/Na) / (1D0 + Keq1*GRH + Keq1*Keq2*GRH**2D0)
        ! update sulfuric acid speed due to RH
        VAPOUR_PROP%c_speed(sa) = SQRT(8D0*Kb*GTEMPK/(pi*VAPOUR_PROP%wet_mass(sa)))
    END IF

END SUBROUTINE UPDATE_MOLECULAR_DIFF_AND_CSPEED


function collision_rate(jj,diameter, mass,VAPOUR_PROP)
    IMPLICIT NONE
    type(vapour_ambient), INTENT(IN) :: VAPOUR_PROP            ! The parameters used are Diffusivity, diffusion diameter, thermal speed, sticking factor (usually 1) and
    integer, INTENT(IN)             :: jj                      ! The index of condensing vapour
    REAL(dp), INTENT(IN)            :: diameter(:), mass(:)    ! Particle diameter and mass
    REAL(dp)                        :: air_free_path,viscosity ! Air properties
    REAL(dp), DIMENSION(n_bins_par) :: collision_rate          ! Collision rate of gas jj and particles
    REAL(dp), DIMENSION(n_bins_par) :: knudsen                 ! Knudsen number []
    REAL(dp), DIMENSION(n_bins_par) :: Diff_par                ! Particle diffusion coefficient [m^2 s^-1]
    REAL(dp), DIMENSION(n_bins_par) :: slip_correction         ! Cunningham slip correction []
    REAL(dp), DIMENSION(n_bins_par) :: speed_p
    REAL(dp), DIMENSION(n_bins_par) :: vapmeanfp,FS_corr,D_vap_eff

    ! viscosity of air and mean free path in air
    viscosity     = 1.8D-5*(GTEMPK/298D0)**0.85D0  ! dynamic viscosity of air
    air_free_path = 2D0*viscosity/(GPRES*SQRT(8D0*Mair/(pi*Rg*GTEMPK))) ! gas mean free path in air

    ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34 pg 407)
    slip_correction = 1D0 + 2D0*air_free_path/diameter * (1.257D0 + 0.4D0*exp(-1.1D0/(2D0*air_free_path/diameter)))

    Diff_par = (kb * GTEMPK * slip_correction) / (3D0 * pi * viscosity * diameter)

    where (mass> 0)
        speed_p = SQRT(8D0*kb*GTEMPK/(pi*mass)) ! Thermal speed of particle m/s
    ELSEWHERE
        speed_p = 0d0
    END where

    vapmeanfp = 3D0*(VAPOUR_PROP%diff(jj) + Diff_par)/SQRT(VAPOUR_PROP%c_speed(jj)**2D0 + speed_p**2D0)

    Knudsen = 2D0 * vapmeanfp/(diameter + VAPOUR_PROP%diff_dia(jj))

    ! Fuchs-Sutugin correction factor for transition regime
    FS_corr        = (0.75 * VAPOUR_PROP%alpha(jj) * (1D0 + Knudsen)) &
                     /(Knudsen**2D0 + Knudsen+0.283*Knudsen*VAPOUR_PROP%alpha(jj) + 0.75d0 * VAPOUR_PROP%alpha(jj))
    D_vap_eff      = (VAPOUR_PROP%diff(jj) + Diff_par) * FS_corr ! m^2/s
    collision_rate = 2D0*pi*(diameter + VAPOUR_PROP%diff_dia(jj)) * D_vap_eff       ! mass transfer coefficient s^-1

END function collision_rate
!



END MODULE aerosol_dynamics
