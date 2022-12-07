! ===================================================================================
! Atmospherically Relevant Chemistry and Aerosol box model
! Copyright (C) 2021  Multi-Scale Modelling group
! Institute for Atmospheric and Earth System Research (INAR), University of Helsinki
! Contact information arca@helsinki.fi
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ===================================================================================


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
SUBROUTINE Condensation_apc(VAPOUR_PROP, conc_vap, dmass, dt_cond, d_dpar,d_vap, kelvin_effect)
    IMPLICIT NONE
    type(vapour_ambient), INTENT(IN)            :: VAPOUR_PROP  ! Properties of condensing vapours
    REAL(dp), INTENT(INOUT)                     :: conc_vap(:)  ! [#/m^3], condensing vapour concentrations, DIM(n_cond_tot)
    REAL(dp), INTENT(INOUT)                     :: dmass(:,:)   ! [kg/m^3] change of mass per particle in particle phase, DIM(n_bins_par, n_cond_tot)
    REAL(dp), INTENT(IN)                        :: dt_cond      ! integration timestep for condensation [s]
    REAL(dp), INTENT(INOUT)                     :: d_dpar(:)    ! relative change in particle diameter [1]
    REAL(dp), INTENT(INOUT)                     :: d_vap(:)     ! relative change in vapor concentrations [1]
    REAL(dp), DIMENSION(:,:), INTENT(IN)        :: Kelvin_effect
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: kohler_effect
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: conc_pp_old
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: conc_pp_eq
    REAL(dp), DIMENSION(n_bins_par,n_cond_tot)  :: conc_pp      ! [#/m^3] Particle phase concentrations, DIM(n_bins_par,n_cond_tot)
    REAL(dp), DIMENSION(n_bins_par)             :: CR           ! [/s] collisions per second, collision rate * concentration
    REAL(dp), DIMENSION(n_cond_tot)             :: conc_tot
    REAL(dp), DIMENSION(n_cond_tot)             :: conc_vap_old
    REAL(dp), DIMENSION(n_bins_par)             :: n_conc
    ! REAL(dp), DIMENSION(n_bins_par)             :: conc_pp_eq_bin
    REAL(dp), DIMENSION(n_bins_par)             :: diameter
    REAL(dp), DIMENSION(n_bins_par)             :: volume
    REAL(dp), DIMENSION(n_bins_par)             :: mass         ! [kg/m^3]
    ! REAL(dp), DIMENSION(n_bins_par)             :: term2
    ! REAL(dp)                                    :: st,term1,term3         ! Surface tension for all compounds
    REAL(dp)                                    :: conc_guess
    REAL(dp)                                    :: sum_org      ! Sum of organic concentration in bin, transient variable
    INTEGER                                     :: ic           ! index for condesables
    INTEGER                                     :: ip           ! index for particles


    ! st = VAPOUR_PROP%surf_tension(1)

    conc_pp = 0d0
    conc_pp = get_composition()

    n_conc = get_conc()
    diameter = get_dp()
    volume = get_volume()
    mass = sum(conc_pp, 2)
    conc_pp_eq = 0d0

    ! Fill the number concentration composition matrix
    kohler_effect = 1d0
    DO ip=1,n_bins_par
        conc_pp(ip,:) = conc_pp(ip,:) * Na / VAPOUR_PROP%molar_mass * n_conc(ip)
        if (use_raoult) THEN
            if (volume(ip)>0.0) then
                sum_org = sum(conc_pp(ip,:),DIM=1, Mask=(VAPOUR_PROP%cond_type==1))
                if (sum_org>0) THEN
                    ! mole fraction of each organic compound in each bin
                    kohler_effect(ip,:)= conc_pp(ip,:)/sum_org
                END if
            END if
        END if
    END DO

    ! original gas phase concentration
    conc_vap_old = conc_vap
    ! original particle phase concentration
    conc_pp_old = conc_pp
    ! vapour phase + particle phase
    conc_tot = conc_vap + sum(conc_pp,1)



    ! Kelvin and Kohler factors.
    kohler_effect = kohler_effect * Kelvin_Effect

    ! Make sure kohler effect for H2SO4 is 1
    kohler_effect(:,VAPOUR_PROP%ind_H2SO4) = Kelvin_Effect(:,VAPOUR_PROP%ind_H2SO4)

    ! Approximate equilibrium concentration (#/m^3) of each compound in each size bin
    DO ip=1,n_bins_par
        conc_pp_eq(ip,1:n_cond_org) = conc_vap_old(1:n_cond_org)*SUM(conc_pp_old(ip,1:n_cond_org)) &
                                        / (Kelvin_Effect(ip,1:n_cond_org)*VAPOUR_PROP%c_sat(1:n_cond_org))
    END DO

    DO ic = 1, n_cond_tot

        ! NOTE GENERIC does not exist in gas phase
        if (ic == VAPOUR_PROP%ind_GENERIC) cycle

        ! Collision rate of particles and gas molecules * concentration -> total number of collisions/s
        CR = n_conc * collision_rate(ic,diameter,mass,VAPOUR_PROP)

        ! Update sulfuric acid condensation sink
        if (ic == VAPOUR_PROP%ind_H2SO4) GCS = sum(CR)

        ! apc scheme here. NOTE that for sulfuric kohler_effect = Kelvin_Effect
        conc_guess = (conc_vap(ic) + dt_cond*sum(CR*kohler_effect(:,ic)*VAPOUR_PROP%c_sat(ic))) / (1D0 + dt_cond*sum(CR))


        ! conc_pp_eq_bin(:) = conc_guess*SUM(conc_pp_old(ip,1:n_cond_org)) &
        !                                     / (Kelvin_Effect(ip,1:n_cond_org)*VAPOUR_PROP%c_sat(1:n_cond_org))


        ! apc scheme cont. NOTE that for sulfuric acid VAPOUR_PROP%c_sat(ic) ~ 0 so the last term will vanish
        conc_pp(:,ic) = conc_pp_old(:,ic) + dt_cond*CR*(MIN(conc_guess,conc_tot(ic)) &
                        - kohler_effect(:,ic)*VAPOUR_PROP%c_sat(ic))

        ! Prevent overestimation of evaporation by setting negative concentrations to zero
        WHERE ( conc_pp(:,ic)<0D0 ) conc_pp(:,ic) = 0D0

        ! Prevents particles to grow over the saturation limit. For organics, no acids or GENERIC included
        IF (ic < VAPOUR_PROP%ind_GENERIC) then
            WHERE (   ( diameter>0d0 )  .and.  ( conc_pp(:,ic)>conc_pp_eq(:,ic) )                 &
                    .AND. ( conc_vap_old(ic)  <  (Kelvin_Effect(:,ic) * VAPOUR_PROP%c_sat(ic)) )  &
                    )
                    conc_pp(:,ic) = conc_pp_eq(:,ic)
            END WHERE
        END IF


        ! Update conc_vap: total conc - particle phase concentration. Sulfuric acid is handled in chemistry if it has H2SO4
        if (ic /= VAPOUR_PROP%ind_H2SO4 .or. H2SO4_ind_in_chemistry<1) conc_vap(ic) = conc_tot(ic) - sum(conc_pp(:,ic))

        ! derive the relative changes in the vapor phase. MIN_CONCTOT_CC_FOR_DVAP is in NML_CUSTOM
        IF ( conc_tot(ic) > MIN_CONCTOT_CC_FOR_DVAP * 1d6 ) d_vap(ic) = ( conc_vap(ic) - conc_vap_old(ic) ) / conc_tot(ic)

        ! Calculate dmass for PSD
        dmass(:,ic) = (conc_pp(:,ic) - conc_pp_old(:,ic))*VAPOUR_PROP%molar_mass(ic) / Na

    END DO


    ! Only apply condensation if there are particles
    do ip = 1, n_bins_par
        if (n_conc(ip)>1d-10) THEN
            dmass(ip,:) = dmass(ip,:) / n_conc(ip)
        else
            dmass(ip,:) = 0d0
        END if

        ! Derive diameter changes for integration time-step optimization
        if ( ip < n_bins_par ) THEN
            ! only if there are particles and if they had a composition in the last timestep
            IF (SUM(conc_pp_old(ip,:)) > 0.d0 .and. n_conc(ip) > 9.9d-7) THEN
                ! print*, dmass(ip,:)
                ! d_dpar(ip) = (SUM(conc_pp_old(ip,:)*VAPOUR_PROP%molar_mass(:)) / Na /n_conc(ip) + SUM(dmass(ip,:))) &
                !             / (SUM(conc_pp_old(ip,:)*VAPOUR_PROP%molar_mass(:)) / Na /n_conc(ip))
                d_dpar(ip) = 1d0 + SUM(dmass(ip,:)) / (SUM(conc_pp_old(ip,:)*VAPOUR_PROP%molar_mass(:)) / Na / n_conc(ip))

                if (d_dpar(ip)<0d0.and.OPTIMIZE_DT) THEN
                    call SET_ERROR(PRC%cch, 'Particles shrink beyond their size')
                    return
                ELSE
                    d_dpar(ip) = ( d_dpar(ip)**(1.d0/3.d0) ) - 1.d0
                END IF

            ELSE
                d_dpar(ip) = 0.d0
            END IF
        END IF
    END DO

END SUBROUTINE Condensation_apc


SUBROUTINE AGING(composition, halflife, dt)
    real(dp), intent(INOUT) :: composition(:,:)
    real(dp), intent(IN   ) :: halflife
    real(dp), intent(IN   ) :: dt
    real(dp)                :: decay, k, agedmass(n_bins_par)
    INTEGER                 :: ic,ip

    agedmass = 0d0
    do ic=1,VAPOUR_PROP%n_cond_org
        k = log(2d0)/halflife
        decay = exp(-k*dt)
        composition(:,ic) = composition(:,ic) * decay
        agedmass = agedmass + composition(:,ic)*(1-decay)
    END DO
    composition(:,VAPOUR_PROP%ind_GENERIC) = composition(:,VAPOUR_PROP%ind_GENERIC) + agedmass

END SUBROUTINE AGING


! -------------------------------------------------------------------------------------------------------------
! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
SUBROUTINE UPDATE_COAG_COEF()
    IMPLICIT NONE
    integer                                     :: i
    REAL(dp), DIMENSION(n_bins_par)             :: diameter, mass
    REAL(dp), DIMENSION(n_bins_par,n_bins_par)  :: sticking_prob ! alpha in Fuchs' beta, S&P (2016) p. 550, eq. 13.56
    REAL(dp), DIMENSION(n_bins_par)             :: slip_correction
    REAL(dp), DIMENSION(n_bins_par)             :: diffusivity,dist
    REAL(dp), DIMENSION(n_bins_par)             :: speed_p
    REAL(dp), DIMENSION(n_bins_par)             :: free_path_p
    REAL(dp), DIMENSION(n_bins_par)             :: Beta_Fuchs
    REAL(dp)                                    :: dyn_visc     ! dynamic viscosity, kg/(m*s)
    REAL(dp)                                    :: l_gas        ! Gas mean free path in air

    sticking_prob = alpha_coa
    diameter      = get_dp()
    mass          = get_mass()

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

    ! Particle mean free path (m)
    free_path_p = 8D0*diffusivity/(pi*speed_p)

    ! Mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)
    dist = (SQRT(2D0)/(3D0*diameter*free_path_p))*((diameter + free_path_p)**3D0 &
           - sqrt((diameter**2D0 + free_path_p**2D0)**(3D0))) - diameter

    DO i = 1,n_bins_par
        ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
        Beta_Fuchs = 1D0/((diameter + diameter(i))/(diameter + diameter(i) &
                        + 2D0*sqrt(dist**2D0 + dist(i)**2D0)) + (8D0/sticking_prob(i,:))*(diffusivity + diffusivity(i)) &
                        /((sqrt(speed_p**2D0+speed_p(i)**2D0))*(diameter + diameter(i))))

        ! coagulation rates between two particles of all size combinations  (m^3/s)
        Gcoag_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i) &
                        + diameter*diffusivity + diameter(i)*diffusivity &
                        + diameter(i)*diffusivity(i))
    END DO

END SUBROUTINE UPDATE_COAG_COEF


SUBROUTINE Coagulation_routine(dconc_coag, dt_coag, d_npar, update_k) ! Add more variables if you need it
    IMPLICIT NONE
    REAL(dp), DIMENSION(n_bins_par,n_bins_par), intent(inout) :: dconc_coag    ! coagulation coefficients [m^3/s]
    REAL(dp), INTENT(IN)                        :: dt_coag ! integration timestep for coagulation [s]
    REAL(dp), INTENT(INOUT)                     :: d_npar(:) ! relative change in particle number concentrations [1]
    REAL(dp), DIMENSION(n_bins_par)             :: n_conc,diameter, mass
    REAL(dp), DIMENSION(n_bins_par)             :: volume
    integer                                     :: i,j,m
    REAL(dp)                                    :: aa
    LOGICAL, intent(INOUT)                      :: update_k

    n_conc = get_conc()

    ! Left here for documentation purposes only. For ARCA paper, to reproduce
    ! Seinfelf & Pandis fig 13.6, which used analytical solution and K_ij = K
    ! if (gtime%sec==0) Gcoag_coef = 1d-15
    ! update_k = .false.

    if (update_k) THEN
        ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
        call UPDATE_COAG_COEF()
        update_k = .false.
    END IF
dconc_coag = 0d0

if (GTIME%savenow) G_COAG_SINK = 0d0

do j=1, n_bins_par
    do m = j, n_bins_par
        if (m==j) then
            aa=0.5_dp
        else
            aa= 1.0_dp
        END if
        if (GTIME%savenow) G_COAG_SINK(j) = G_COAG_SINK(j) + aa * Gcoag_coef(j,m) * n_conc(m)
        IF (n_conc(j) > 1.d0 .and. n_conc(m) > 1.d0) dconc_coag(j,m) = aa * Gcoag_coef(j,m) * n_conc(j) * n_conc(m) * dt_coag ! 1/m^3 '/ dt
    END DO

    ! Check whether changes are worth to be applied:
    IF (n_conc(j) > 1.d0 .and. sum(dconc_coag(j,j:)) > 1.d-20) THEN
        d_npar(j) = sum(dconc_coag(j, j:)) / n_conc(j)
    ELSE
        d_npar(j) = 0.d0
        dconc_coag(j,j:) = 0.d0
    END IF
END DO

END SUBROUTINE Coagulation_routine


! Loss rate calculation by P. Roldin
SUBROUTINE deposition_velocity(d_p,ustar,Av,Au,Ad,V_chamber,T,p,E_field,dn_par, dt,k_dep)
    IMPLICIT NONE
    REAL(dp), DIMENSION(n_bins_par), INTENT(in) :: d_p
    REAL(dp), INTENT(in)                        :: ustar,Av,Au,Ad,V_chamber,T,p,dt, E_field
    REAL(dp), DIMENSION(:), INTENT(inout)       :: dn_par ! s^-1
    REAL(dp), DIMENSION(n_bins_par)             :: D,vs,Sc,r,aa,bb,Ii,vdv,vdu,vdd,Cc,dens !,ve
    REAL(dp), INTENT(INOUT)                     :: k_dep(:) ! s^-1
    REAL(dp)                                    :: dyn_visc,kin_visc,l_gas,dens_air

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
    if (GTIME%printnow) print FMT_SUB, 'Particle loss rate (geom.mean):'//TRIM(f2chr( EXP(SUM(LOG(k_dep))/n_bins_par)))&
                    //'/s, '//TRIM(f2chr( EXP(SUM(LOG(k_dep))/n_bins_par)*3600d0))//'/h'
    dn_par = get_conc() * ( 1 - EXP(-k_dep*dt)) ! calculates how the particle number concentration changes in each size bin due to dry deposition

END SUBROUTINE deposition_velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE UPDATE_MOLECULAR_DIFF_AND_CSPEED(VAPOUR_PROP)
    IMPLICIT NONE
    type(vapour_ambient), INTENT(INOUT) :: VAPOUR_PROP  ! VapourProp file that gets updated
    INTEGER                             :: om,sa        ! short for indices

    sa = VAPOUR_PROP%ind_H2SO4
    om = VAPOUR_PROP%n_cond_org

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

subroutine avg3(arr, refdvap, irefdvap)
    implicit none
    REAL(dp) :: refdvap
    integer  :: irefdvap
    real(dp),intent(in)::arr(:)
    REAL(dp) :: ret(2)
    integer :: i
    ret(1) = avg3neg(pack( arr,arr<0d0 ))
    ret(2) = avg3pos(arr)
    refdvap = ret(maxloc(abs(ret),1))
    if (maxloc(abs(ret),1)==1) THEN
        irefdvap = minloc(arr,1)
    ELSE
        irefdvap = maxloc(arr,1)
    END IF
end subroutine avg3

real(dp) function avg3pos(arr)
    ! returns average of three largest (absolute magnitude) values of arr
    implicit none
    real(dp),intent(in)::arr(:)
    real(dp) :: Z(size(arr,1)-1)
    integer :: ii(2), l,i
    l = size(arr,1)
    avg3pos = 0d0
    if (l==0) return
    if (l<4) THEN
        avg3pos = sum(arr)/l
    ELSE
        avg3pos =  (maxval(arr,1))
        ii(1) = (maxloc(arr,1))
        Z = ( arr( pack([(i,i=1,l)], [(i/=ii(1),i=1,l)]) )     )
        avg3pos = avg3pos + maxval( Z )
        ii(2) =  maxloc( Z ,1)
        avg3pos =  avg3pos + maxval(  Z( pack([(i,i=1,l-1)], [(i/=ii(2),i=1,l-1)]) )     )
        avg3pos = avg3pos/3_dp
    END IF
end function avg3pos

real(dp) function avg3neg(arr)
    ! returns average of three largest (absolute magnitude) values of arr
    implicit none
    real(dp),intent(in)::arr(:)
    real(dp) :: Z(size(arr,1)-1)
    integer :: ii(2), l,i
    l = size(arr,1)
    avg3neg = 0d0
    if (l==0) return
    if (l<4) THEN
        avg3neg = sum(arr)/l
    ELSE
        avg3neg =  (minval(arr,1))
        ii(1) = (minloc(arr,1))
        Z = ( arr( pack([(i,i=1,l)], [(i/=ii(1),i=1,l)]) )     )
        avg3neg = avg3neg + minval( Z )
        ii(2) =  minloc( Z ,1)
        avg3neg =  avg3neg + minval(  Z( pack([(i,i=1,l-1)], [(i/=ii(2),i=1,l-1)]) )     )
        avg3neg = avg3neg/3_dp
    END IF
end function avg3neg


END MODULE aerosol_dynamics
