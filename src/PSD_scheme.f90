MODULE PSD_scheme
  USE second_precision, ONLY: dp
  USE constants
  USE INPUT

  IMPLICIT NONE

  ! START: variables that will be defined outside
  INTEGER :: n_cond_tot  ! number of species that can go to the particle phase
  INTEGER :: n_cond_org  ! number of organic species that can go to the particle phase
  ! END variables that will be defined outside

  REAL(dp), ALLOCATABLE :: dconc_coag(:,:)  ! coagulation: collision number matrix: nr_bins * nr_bins
  REAL(dp), ALLOCATABLE :: pre_Kelvin(:,:)      ! Precalculated Kelvin terms
  REAL(dp), ALLOCATABLE :: Kelvin_eff(:,:)      ! Precalculated Kelvin terms
  REAL(dp), ALLOCATABLE :: dmass(:,:)       ! change in particle mass (due to condensation or mixing) (nr_bins,n_cond_tot)
  REAL(dp), ALLOCATABLE :: dconc_dep_mix(:) ! change in particle concentration (e.g. reduced coagulatio or mixing) (nr_bins)
  REAL(dp), ALLOCATABLE :: G_COAG_SINK(:) ! change in particle concentration (e.g. reduced coagulatio or mixing) (nr_bins)
  REAL(dp)              :: mix_ratio        ! gives the rate ratio: added volume over present volume per timestep
  REAL(dp)              :: bin_ratio        ! gives relative bin width dp(2)/dp(1)
  CHARACTER(len=15)     :: process          ! defines the process that passes information to subroutine Mass_Number_Change (coagulation, condensation, mixing)
  type(PSD)             :: current_PSD      ! Main PSD container. This variable stores the current timestep concentrations
  type(PSD) :: new_PSD, mix_PSD, interm_PSD ! Variables that store PSD values during the calculations
  type(PSD) :: old_PSD                      ! This is used to save the current state and to restore it in case of an error related to timestep handling

CONTAINS

! ==================================================================================
! This subroutine calls other subroutines from module PSD to:
!  1) extract some infromation from input
!  2) allocate the array sizes for PSD
!  3) generate the PSD related arrays
SUBROUTINE initialize_PSD
    IMPLICIT NONE

    CALL PSD_get_input()      ! get some variables from input module needed for PSD representation
    CALL PSD_Allocate()       ! allocate the variables representing the PSD
    CALL GeneratePSDarrays()  ! Generate the PSD arrays for modelling
    new_PSD = current_PSD
    mix_PSD = current_PSD
    interm_PSD = current_PSD

END SUBROUTINE initialize_PSD


! ==================================================================================
! This subroutine picks inputs needed from input module to initialize the PSD module
! ==================================================================================
SUBROUTINE PSD_get_input()
    IMPLICIT NONE
    INTEGER :: nr_noncond = 0      ! number of species on the particle that are non-volatile

    ! SET mode of particle size distribution representation:
    !   PSD_style: options for particle size distribution representation
    !   1 ... fully stationary (working)
    !   2 ... reserved for fully stationary(coagulation)/full moving(growth)
    !   3 ... reserved for fully moving
    !   4 ... reserved for hybrid
    current_PSD%PSD_style = PSD_MODE

    ! Set the number of bins for model array:
    current_PSD%nr_bins = n_bins_par
    ! Set the range the particle array should feature [m]
    current_PSD%dp_range = (/min_particle_diam,max_particle_diam/)

    !set the number of species considered in the particle phase:
    !H2SO4 is added to the particle phase => +1
    !non condensables initially on the particle phase
    !vapour_number is the stuff that has a quantifiable vapour pressure > 0
    if (ALLOCATED(XTRAS)) nr_noncond = size(XTRAS)
    n_cond_tot = VAPOUR_PROP%n_condtot
    n_cond_org = VAPOUR_PROP%n_condorg


  END SUBROUTINE PSD_get_input


! ===========================================================================
! This subroutine allocates the particle size distribution relevant variables
! based on choice of representation type
! Further it might be used to minimizes the dimensions of varibales used in
! other PSD representation types to save computation time
! ===========================================================================
SUBROUTINE PSD_Allocate()
    IMPLICIT NONE

    ALLOCATE(pre_Kelvin(n_bins_par,n_cond_tot))
    ALLOCATE(Kelvin_eff(n_bins_par,n_cond_tot))

    ! ALLOCATE Change vectors for FS and MA
    IF (current_PSD%PSD_style == 1 .or. current_PSD%PSD_style == 2 ) THEN
        ALLOCATE(dmass(current_PSD%nr_bins,n_cond_tot))
        ALLOCATE(dconc_dep_mix(current_PSD%nr_bins))
        ALLOCATE(G_COAG_SINK(current_PSD%nr_bins))
        ALLOCATE(dconc_coag(current_PSD%nr_bins,current_PSD%nr_bins))

    END IF

    ! FULLY STATIONARY representation !
    IF (current_PSD%PSD_style == 1) then
        ALLOCATE(current_PSD%diameter_fs(current_PSD%nr_bins))
        ALLOCATE(current_PSD%volume_fs(current_PSD%nr_bins))
        ALLOCATE(current_PSD%conc_fs(current_PSD%nr_bins))
        ALLOCATE(current_PSD%composition_fs(current_PSD%nr_bins,n_cond_tot))
        ALLOCATE(current_PSD%density_fs(n_cond_tot))
        ALLOCATE(current_PSD%particle_density_fs(current_PSD%nr_bins))
        ALLOCATE(current_PSD%particle_mass_fs(current_PSD%nr_bins))
        ALLOCATE(current_PSD%dp_dry_fs(current_PSD%nr_bins))

    ! MOVING AVERAGE, FIXED GRID representation !
    ELSE IF (current_PSD%PSD_style == 2) THEN
        ALLOCATE(current_PSD%diameter_ma(current_PSD%nr_bins))
        ALLOCATE(current_PSD%volume_ma(current_PSD%nr_bins))
        ALLOCATE(current_PSD%conc_ma(current_PSD%nr_bins))
        ALLOCATE(current_PSD%composition_ma(current_PSD%nr_bins,n_cond_tot))
        ALLOCATE(current_PSD%density_ma(n_cond_tot))
        ALLOCATE(current_PSD%grid_ma(current_PSD%nr_bins+1))
        ALLOCATE(current_PSD%particle_density_ma(current_PSD%nr_bins))
        ALLOCATE(current_PSD%particle_mass_ma(current_PSD%nr_bins))
        ALLOCATE(current_PSD%dp_dry_ma(current_PSD%nr_bins))

    ELSE
        PRINT*,'  This PSD representation is not yet defined:', current_PSD%PSD_style
        STOP
    END IF
END SUBROUTINE PSD_Allocate


! ===========================================================================
! Generates the particle related array content based on the choice of representation type
! Concentration is set to zero here -> it is set to input values in subroutine "GeneratePSDfromInput"
SUBROUTINE GeneratePSDarrays()
    IMPLICIT NONE
    INTEGER :: i  ! integer for incrementation

    ! --------------------------
    ! FULLY STATIONARY METHOD!
    IF (current_PSD%PSD_style == 1) THEN
        ! Define diameter ratio from lower and upper limits
        bin_ratio = (log(current_PSD%dp_range(2)/current_PSD%dp_range(1)))/(current_PSD%nr_bins-1)
        ! Define diameter array as exponent of a geometric series
        current_PSD%diameter_fs =  exp([(i*bin_ratio + log(current_PSD%dp_range(1)), i=0,current_PSD%nr_bins-1)])
        current_PSD%diameter_fs(1) = current_PSD%dp_range(1)
        ! Define volume array:
        current_PSD%volume_fs = 1D0/6D0 * pi * (current_PSD%diameter_fs)**3.d0   ! Single particle volume (m^3)
        ! particle_density needed for diffusion
        current_PSD%particle_density_fs = 1.4D3
        ! particle_density needed for diffusion
        current_PSD%particle_mass_fs = current_PSD%volume_fs * current_PSD%particle_density_fs
        ! partice phase density for condensables
        current_PSD%density_fs=VAPOUR_PROP%density
        ! get bin ratio:
        bin_ratio = current_PSD%diameter_fs(2)/current_PSD%diameter_fs(1)
        ! Initialize concentration
        current_PSD%conc_fs = 0.d0

    ! --------------------------
    ! MOVING AVERAGE, FIXED GRID
    ELSE IF (current_PSD%PSD_style == 2) THEN

        ! Define diameter array:
        current_PSD%diameter_ma(1) = current_PSD%dp_range(1)
        do i = 2,current_PSD%nr_bins
            current_PSD%diameter_ma(i) = current_PSD%diameter_ma(i-1) &
                                        * (current_PSD%dp_range(2) / current_PSD%dp_range(1)) &
                                        ** (1.D0/(current_PSD%nr_bins-1))
        end do
        ! Define VOLUME array:
        current_PSD%volume_ma = 1D0/6D0 * pi * current_PSD%diameter_ma**3.d0
        ! Set up grid array (i.e. bin borders):
        current_PSD%grid_ma(1:current_PSD%nr_bins) = current_PSD%diameter_ma &
          * (1.d0 - 0.5d0 * (current_PSD%diameter_ma(2) / current_PSD%diameter_ma(1) - 1.d0))
        current_PSD%grid_ma(current_PSD%nr_bins+1) = current_PSD%diameter_ma(current_PSD%nr_bins) &
          * (1.d0 + 0.5d0 * (current_PSD%diameter_ma(2) / current_PSD%diameter_ma(1) - 1.d0))
        ! Particle_density needed for diffusion
        current_PSD%particle_density_ma = 1.4D3
        ! Particle_density needed for diffusion
        current_PSD%particle_mass_ma = current_PSD%volume_ma * current_PSD%particle_density_ma
        ! Partice phase density for condensables
        current_PSD%density_ma=VAPOUR_PROP%density
        ! Get bin ratio:
        bin_ratio = current_PSD%grid_ma(2)/current_PSD%grid_ma(1)
        ! Initialize concentration
        current_PSD%conc_ma = 0.d0

    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF

END SUBROUTINE GeneratePSDarrays


! ===========================================================================
! Generates y distribution as a function of diameter from input (dp_from,
! conc_from by fitting the property y (e.g. particle concentration) Note:
! monodisperse peaks are captured Important: fits tend to underestimate the
! concentration if not monodisperse: if dNdlogdp(i) > 0. and dNdlogdp(i +1 or -1) = 0.
! Then the fit in between is 0!
! Note this routine is not tested for cases where input size resolution is higher than in model
SUBROUTINE GeneratePSDfromInput(dp_from,conc_from,conc_out)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: dp_from(:), conc_from(:)  !arrays: 1) dp  2)other property as input for fitting
    REAL(dp), INTENT(INOUT) :: conc_out(:) ! Final output, fitted concentration array

    REAL(dp) :: dp_sim(current_PSD%nr_bins)
    REAL(dp) :: ddNdlogdp   !difference in dNdlogdp between two channels
    REAL(dp) :: dp_diff(size(dp_from))  !parameter to capture difference between model diameter array and input diameter

    LOGICAL :: mono         ! if mono is true -> capture monodisperse peak, else: ignore (to avoid to capture it several times)
    INTEGER :: n_channels   ! Dimensions of the input PSD
    INTEGER :: channel      ! bin number which is close and smaller than to input diameter
    INTEGER :: j,k          ! integers for loops

    conc_out = 0d0

    dp_sim = get_dp()
    n_channels = size(dp_from)

    ! -------------------------------------------
    ! Fully stationary and moving average  method
    IF (current_PSD%PSD_style == 1 .or. current_PSD%PSD_style == 2) THEN
        mono = .true.
        DO j = 1, current_PSD%nr_bins
            DO k = 1, n_channels
                dp_diff(k) = dp_sim(j) - dp_from(k)
                IF (dp_diff(k) < 0.d0) dp_diff(k) = 1.d10
            END DO

            ! in case the diameter is smaller or larger than the smallest or largest channel diameter => bin conc = 0.d0
            ! else: assign lower channel for fitting
            IF (minval(dp_diff) > 1.d9) THEN
                conc_out(j) = 0.d0
                channel = 0
                mono = .true.
            ELSE
                channel = minloc((dp_diff), DIM=1)
            END IF

            IF (channel > 0 .and. channel < n_channels) THEN
                IF (conc_from(channel) > 0.d0 .and. conc_from(channel+1) > 0.d0) THEN
                    ! Input peak ranges over at least 2 input channels -> fit required
                      ddNdlogdp = conc_from(channel+1) - conc_from(channel)
                      conc_out(j) = conc_from(channel) + (ddNdlogdp) * LOG10(dp_sim(j) &
                                      & / dp_from(channel)) / LOG10(dp_from(channel+1) / dp_from(channel))

                      mono = .true.
                END IF
                IF (conc_from(channel) <= 0.d0) mono = .true.
            END IF

            IF (channel > 0) THEN
                IF (mono) THEN ! Check for monodisperse peaks
                    IF (channel == 1) THEN
                        IF (conc_from(channel) > 0.d0 .and. conc_from(channel+1) <= 0.d0) THEN
                            ! In case there is a monodisperse peak in channel 1:
                            conc_out(j) = conc_from(channel)
                            mono = .false.
                        END IF
                    END IF
                    IF (channel > 1 .and. channel < n_channels ) THEN
                        IF (conc_from(channel-1) <= 0.d0 .and. conc_from(channel) > 0.d0 .and. conc_from(channel+1) <= 0.d0) THEN
                            ! In case there is a monodisperse peak somewhere in the middle:
                            mono = .false.
                            conc_out(j) = conc_from(channel)
                        END IF
                    END IF
                    IF (channel == n_channels) THEN
                        IF (conc_from(channel-1) <= 0.d0 .and. conc_from(channel) > 0.d0) THEN
                            ! In case there is a monodisperse peak at the end:
                            mono = .false.
                            conc_out(j) = conc_from(channel)

                        END IF
                    END IF
                END IF
            END IF

        END DO

        ! Convert from dNdlogdp to concentration:
        conc_out = conc_out * LOG10(dp_sim(2)/dp_sim(1))

    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF

END SUBROUTINE GeneratePSDfromInput



! ===================================================================================================================
! Applies changes to the particle size distribution amd composition according to input-change-vector:
! 1) dconc_coag matrix for coagulation;
! 2) dconc_dep_mix and dmass for mixing
! 3) dmass for condensation
! Based on current and absolute change in single particle mass/number [units similar to composition and concentration]
! the new distribution is determined
! Subroutine input is change-array dmass (n_cond_tot,nr_bins) for mass
! Subroutine input is change-array dconc_dep_mix(nr_bins) or dconc_coag(nr_bins,nr_bins) for number
! ===================================================================================================================
SUBROUTINE Mass_Number_Change(process)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)  :: process

    ! FULLY STATIONARY METHOD & MOVING AVERAGE FIXED GRID
    IF (current_PSD%PSD_style == 1 .or. current_PSD%PSD_style == 2) THEN

        IF (process == 'condensation') THEN ! Condensation
            CALL PSD_Change_Condensation()

        ELSE IF (process == 'coagulation') THEN !Coagulation
            CALL PSD_Change_coagulation()

        ELSE IF (process == "mixing") THEN !Mix it
            CALL PSD_Change_mixing()

        ELSE IF (process == 'deposition') THEN !Deposition
            CALL PSD_Change_deposition()
        END IF

    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF
END SUBROUTINE Mass_Number_Change


! ====================================================================================================================
! Determines the changes in the PSD due to condensation based on dmass. Condensation is continuous process: all
! particles grow at the same rate. ATTENTION: The composition in the largest bin is wrong if growth happens via a==0!!
! This might cause errors -> find some solution at some point
! ....................................................................................................................
SUBROUTINE PSD_Change_condensation()
    IMPLICIT NONE
    INTEGER ::   ii     !some integer for looping
    !------------------------
    !FULLY STATIONARY METHOD!
    IF (current_PSD%PSD_style == 1) THEN

        ! Initiate new_PSD
        new_PSD%composition_fs = current_PSD%composition_fs  ! Set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_fs = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        !Find new volume in case of condensation
        DO ii = 1, current_PSD%nr_bins
            mix_PSD%volume_fs(ii) = current_PSD%volume_fs(ii) + &
                                       SUM(dmass(ii,:) / current_PSD%density_fs(:))

           where(dmass(ii,:) + current_PSD%composition_fs(ii,:)<0d0)
               mix_PSD%composition_fs(ii,:) = 0d0
           ELSEWHERE
               mix_PSD%composition_fs(ii,:) = dmass(ii,:) + current_PSD%composition_fs(ii,:)
           end where

            mix_PSD%conc_fs(ii) = current_PSD%conc_fs(ii)
            IF ( mix_PSD%conc_fs(ii) > 1.d-100 ) CALL bin_redistribute_fs(ii)
            mix_PSD%conc_fs(ii) = 0.d0

        END DO

    !---------------------------
    ! Moving average, fixed grid
    ELSE IF  (current_PSD%PSD_style == 2) THEN
        ! Initiate new_PSD
        new_PSD%composition_ma = current_PSD%composition_ma  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_ma = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        new_PSD%diameter_ma = current_PSD%diameter_ma
        ! Find new volume in case of condensation
        DO ii = 1, current_PSD%nr_bins

            mix_PSD%volume_ma(ii) = SUM((current_PSD%composition_ma(ii,:) + dmass(ii,:)) / current_PSD%density_ma(:))! + SUM(dmass(ii,:) / current_PSD%density_ma(:))

            mix_PSD%composition_ma(ii,:) = dmass(ii,:) + current_PSD%composition_ma(ii,:)

            mix_PSD%conc_ma(ii) = current_PSD%conc_ma(ii)
            where(mix_PSD%composition_ma(ii,:)<0d0) mix_PSD%composition_ma(ii,:) = 0d0
            IF ( mix_PSD%conc_ma(ii) > 0d0 ) CALL bin_redistribute_ma(ii)
            mix_PSD%conc_ma(ii) = 0.d0
        END DO
    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF
END SUBROUTINE PSD_Change_condensation


! ====================================================================================================================
! This subroutine uses the dconc_coag array (nr_bins * nr_bins) which is determined by coagulation
! In the first row there is 1 entry at (1,1): concentration loss of smallest particles due to collision with each other
! In line n there are n entries(n,1:n): concentration loss of smallest particles due to collision with particle of size
! n, second smallest,... last entry: collisions of n with n
! Result is the new particle size distribution and composition: new_PSD%conc_fs and new_PSD%composition_fs
SUBROUTINE PSD_Change_coagulation
    IMPLICIT NONE
    INTEGER :: ii,jj  !some integer for incrementation

    !-----------------------
    ! FULLY STATIONARY METHOD
    IF (current_PSD%PSD_style == 1) THEN
        ! Initiate new_PSD and mix_PSD
        new_PSD%composition_fs = current_PSD%composition_fs  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_fs = current_PSD%conc_fs  ! set the initial value of new particle concentration to the current
        mix_PSD%conc_fs = 0.d0  !initial value is zero -> content is determined below
        ! Apply changes for all combinations of i and j
        DO ii = 1, current_PSD%nr_bins
            DO jj = ii, current_PSD%nr_bins
                IF (dconc_coag(ii,jj) > 1.d-100) THEN
                    IF (.not. use_speed) THEN
                        IF (dconc_coag(ii,jj) > new_PSD%conc_fs(ii)) dconc_coag(ii,jj) = new_PSD%conc_fs(ii) * 0.9
                        IF (dconc_coag(ii,jj) > new_PSD%conc_fs(jj)) dconc_coag(ii,jj) = new_PSD%conc_fs(jj) * 0.9
                    END IF
                    ! Reduce the new particle concentration by the number of particles that are lost by coagulation in i and j -> they will be added later to the new bin
                    new_PSD%conc_fs(ii) = new_PSD%conc_fs(ii) - dconc_coag(ii,jj)  !reduce number in ii
                    new_PSD%conc_fs(jj) = new_PSD%conc_fs(jj) - dconc_coag(ii,jj)  !reduce number in jj (if ii=ij we have to reduce twice (which is done here) as 1 collision removes 2 particles)

                    ! Determine new mass compositions: (total mass of collision products (i+j)
                    mix_PSD%composition_fs(ii,:) = (current_PSD%composition_fs(ii,:) + current_PSD%composition_fs(jj,:))  !composition of of collision result bins: i + j
                    ! Update concentration in the mix_PSD
                    mix_PSD%conc_fs(ii) = dconc_coag(ii,jj)
                    ! Determine volume of the mixing aerosol
                    mix_PSD%volume_fs(ii) = SUM(mix_PSD%composition_fs(ii,:) / current_PSD%density_fs(:))   !Determine the volume of particles in the mix distribution

                    CALL bin_redistribute_fs(ii)
                    ! if (PRC%err) return
                    mix_PSD%conc_fs(ii) = 0.d0  !reset value to zero -> changes have been applied
              END IF
            END DO
        END DO
        if (minval(new_PSD%composition_fs)<0d0.and.Use_speed) THEN
            CALL SET_ERROR(PRC%coa, 'Getting negative composition in Mass_Number_Change')
        END IF

    !-----------------------
    ! Moving average, fixed grid  !!
    ELSE IF  (current_PSD%PSD_style == 2) THEN
        ! Initiate new_PSD and mix_PSD
        new_PSD%composition_ma = current_PSD%composition_ma ! Set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_ma = current_PSD%conc_ma  ! Set the initial value of new particle concentration to the current
        mix_PSD%conc_ma = 0.d0  ! Initial value is zero -> content is determined below
        new_PSD%diameter_ma = current_PSD%diameter_ma
        ! Apply changes for all combinations of i and j
        DO ii = 1, current_PSD%nr_bins
            DO jj = ii, current_PSD%nr_bins
                IF (dconc_coag(ii,jj) > 0d0) THEN
                    IF (.not. use_speed) THEN
                        IF (dconc_coag(ii,jj) > new_PSD%conc_ma(ii)) dconc_coag(ii,jj) = new_PSD%conc_ma(ii) * 0.9
                        IF (dconc_coag(ii,jj) > new_PSD%conc_ma(jj)) dconc_coag(ii,jj) = new_PSD%conc_ma(jj) * 0.9
                    END IF

                    ! Reduce the new particle concentration by the number of particles that are lost by coagulation in i and j -> they will be added later to the new bin
                    new_PSD%conc_ma(ii) = new_PSD%conc_ma(ii) - dconc_coag(ii,jj)  !reduce number in i
                    new_PSD%conc_ma(jj) = new_PSD%conc_ma(jj) - dconc_coag(ii,jj)  !reduce number in j (if i=j we have to reduce twice (which is done here) as 1 collision removes 2 particles)

                    ! Determine new mass compositions: (total mass of collision products (i+j) + total mass already in mix bin i) / devided by sum of concentration (collisions + mix)
                    mix_PSD%composition_ma(ii,:) = (current_PSD%composition_ma(ii,:) + current_PSD%composition_ma(jj,:))  !composition of of collision result bins: i + j

                    ! Update concentration in the mix_PSD
                    mix_PSD%conc_ma(ii) = dconc_coag(ii,jj)

                    ! Determine volume of the mixing aerosol
                    mix_PSD%volume_ma(ii) = SUM(mix_PSD%composition_ma(ii,:) / current_PSD%density_ma(:))   !Determine the volume of particles in the mix distribution

                    CALL bin_redistribute_ma(ii)
                    if (PRC%err) return
                    mix_PSD%conc_ma(ii) = 0.d0  !reset value to zero -> changes have been applied
                END IF
            END DO
        END DO
    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF
END SUBROUTINE PSD_Change_coagulation


! ====================================================================================================================
! This subroutine uses the dconc_dep_mix and dmass vectors provided by mixing subroutine
! dconc_dep_mix give the size dependent particle number that is added: (nr_bins)
! dmass gives the composition of those particles: (nr_bins,n_cond_tot)
! Result is the new particle size distribution and composition: new_PSD%conc_fs and new_PSD%composition_fs
! mix_ratio dtermines the ratio between present and added volume
SUBROUTINE PSD_Change_mixing

    IMPLICIT NONE
    INTEGER ::  ii  !some integer for incrementation

    !------------------------
    !FULLY STATIONARY METHOD!
    IF (current_PSD%PSD_style == 1) THEN

        !Initialize new_PSD and mix_PSD:

        new_PSD%composition_fs = current_PSD%composition_fs  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        mix_PSD%conc_fs = 0.d0  !initial value is zero -> content is determined below
        new_PSD%conc_fs = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        !apply changes for all bins
        DO ii = 1, current_PSD%nr_bins
            IF (dconc_dep_mix(ii) > 1.d-100) THEN
                !Set mass of mixing aerosol:
                mix_PSD%composition_fs(ii,:) = dmass(ii,:)  !composition of i
                !set concentration in the mix_PSD
                mix_PSD%conc_fs(ii) = dconc_dep_mix(ii) * ABS(mix_ratio)
                !Determine new composition and particle concentrations
                new_PSD%conc_fs(ii) = (current_PSD%conc_fs(ii) + mix_PSD%conc_fs(ii) * ABS(mix_ratio)) / &
                                     (1.d0 + MAX(0d0,mix_ratio))
                new_PSD%composition_fs(ii,:) = (current_PSD%composition_fs(ii,:) * current_PSD%conc_fs(ii) + &
                                              mix_PSD%composition_fs(ii,:) * ABS(mix_ratio) * mix_PSD%conc_fs(ii)) / &
                                              (current_PSD%conc_fs(ii) + ABS(mix_ratio) * mix_PSD%conc_fs(ii))
                mix_PSD%conc_fs(ii) = 0.d0  !reset value to zero -> changes have been applied
            ELSE
                new_PSD%conc_fs(ii) = current_PSD%conc_fs(ii) !nothing to mix with => new concentration is current concentration; composition is already new (see subroutine Mass_Number_Change)
            END IF
        END DO

    !---------------------------
    ! Moving average, fixed grid
    ELSE IF  (current_PSD%PSD_style == 2) THEN
        ! Initialize new_PSD and mix_PSD:
        new_PSD%composition_ma = current_PSD%composition_ma  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no particles )
        mix_PSD%conc_ma = 0.d0  !initial value is zero -> content is determined below
        new_PSD%conc_ma = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        new_PSD%diameter_ma = current_PSD%diameter_ma
        mix_PSD%diameter_ma(1:current_PSD%nr_bins-1) = (current_PSD%grid_ma(1:current_PSD%nr_bins-1) + &
                current_PSD%grid_ma(2:current_PSD%nr_bins)) / 2.d0  !the mixing aerosol diameter is in between the grid-limits
                ! XXX This will give 1 bin less than total and reallocate the variable accordingly, probably not what is wanted

        ! Apply changes for all bins
        DO ii = 1, current_PSD%nr_bins
            IF (dconc_dep_mix(ii) > 1.d-100) THEN
                ! Set mass of mixing aerosol:
                mix_PSD%composition_ma(ii,:) = dmass(ii,:)  !composition of i
                ! Set concentration in the mix_PSD
                mix_PSD%conc_ma(ii) = dconc_dep_mix(ii) * ABS(mix_ratio)

                ! print*, 'in PSD', mix_PSD%conc_ma(ii)
                ! Determine new composition and particle concentrations
                new_PSD%conc_ma(ii) = (current_PSD%conc_ma(ii) + mix_PSD%conc_ma(ii) * ABS(mix_ratio)) / &
                             (1.d0 + max(0d0,mix_ratio))

                 ! print*, 'in PSD', new_PSD%conc_ma(ii), mix_PSD%conc_ma(ii) * ABS(mix_ratio),  1.d0 + max(0d0,mix_ratio)
                ! New composition based on total mass composition of both aerosols of bin ii
                new_PSD%composition_ma(ii,:) = (current_PSD%composition_ma(ii,:) * current_PSD%conc_ma(ii) + &
                                      mix_PSD%composition_ma(ii,:) * ABS(mix_ratio) * mix_PSD%conc_ma(ii)) / &
                                      (current_PSD%conc_ma(ii) + ABS(mix_ratio) * mix_PSD%conc_ma(ii))
                ! New diameter based on total volume of both aerosols of bin ii
                new_PSD%diameter_ma(ii) = ((current_PSD%diameter_ma(ii) ** 3.d0 * current_PSD%conc_ma(ii) + &
                                  mix_PSD%diameter_ma(ii) ** 3.d0 * ABS(mix_ratio) * mix_PSD%conc_ma(ii)) / &
                                  (current_PSD%conc_ma(ii) + ABS(mix_ratio) * mix_PSD%conc_ma(ii))) ** (1.d0/3.d0)
                mix_PSD%conc_ma(ii) = 0.d0  ! Reset value to zero -> changes have been applied
            ELSE
                new_PSD%conc_ma(ii) = current_PSD%conc_ma(ii) ! Nothing to mix with => new concentration is current concentration; composition is already new (see subroutine Mass_Number_Change)
            END IF
        END DO
    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF
END SUBROUTINE PSD_Change_mixing


! ==================================================================
! This subroutine uses the dconc_dep_mix to deposit particles
! Result is the new particle size distribution: new_PSD%conc_fs
SUBROUTINE PSD_Change_deposition
    IMPLICIT NONE
    INTEGER :: ii  ! some integer for incrementation

    IF (current_PSD%PSD_style == 1) THEN
        !---------------------------
        ! FULLY STATIONARY METHOD
        ! Initialize new_PSD and mix_PSD:
        new_PSD%composition_fs = current_PSD%composition_fs  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_fs = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        !apply changes for all bins
        DO ii = 1, current_PSD%nr_bins
            IF (dconc_dep_mix(ii) > 1.d-100) THEN
                new_PSD%conc_fs(ii) = current_PSD%conc_fs(ii) - dconc_dep_mix(ii) !reduce number of particles due to deposition; composition remains the same
            ELSE IF (dconc_dep_mix(ii) < 0.d0) THEN
                ! Error: negative deposition
                Print*, 'Error in deposition: negative deposition'
                Print*, 'Bin/value:', ii, dconc_dep_mix(ii)
            ELSE
                new_PSD%conc_fs(ii) = current_PSD%conc_fs(ii) !nothing to mix with => new concentration is current concentration; composition is already new (see subroutine Mass_Number_Change)
            END IF
        END DO
    ELSE IF  (current_PSD%PSD_style == 2) THEN
        !---------------------------
        ! Moving average, fixed grid
        ! Initialize new_PSD and mix_PSD:
        new_PSD%composition_ma = current_PSD%composition_ma  ! Set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
        new_PSD%conc_ma = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
        new_PSD%diameter_ma = current_PSD%diameter_ma
        ! Apply changes for all bins
        DO ii = 1, current_PSD%nr_bins
            IF (dconc_dep_mix(ii) > 1.d-100) THEN
                new_PSD%conc_ma(ii) = current_PSD%conc_ma(ii) - dconc_dep_mix(ii) !reduce number of particles due to deposition; composition remains the same
            ELSE IF (dconc_dep_mix(ii) < 0.d0) THEN
                ! Error: negative deposition
                Print*, 'Error in deposition: negative deposition'
                Print*, 'Bin/value:', ii, dconc_dep_mix(ii)
            ELSE
                new_PSD%conc_ma(ii) = current_PSD%conc_ma(ii) !nothing to mix with => new concentration is current concentration; composition is already new (see subroutine Mass_Number_Change)
            END IF
        END DO
    ELSE
        print*,'choose other form of representation:'
        print*,current_PSD%PSD_style,'is not defined yet'
    END IF
END SUBROUTINE PSD_Change_deposition


! =================================================================
! This subroutine a) redistributes the particles that are newly formed by coagulation into the bins or b) redistributes
! particles that changed their seize due to condensation/evaporation
! It also determines the new composition. Results are saved in new_PSD%conc_fs and new_PSD%composition_fs
! FULLY STATIONARY METHOD!
! =================================================================
SUBROUTINE bin_redistribute_fs(ind)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ind   ! The bin which is to be redistributed -> it generally has different diameter than in diameter_fs!
    INTEGER             :: n_b   ! Short for number of bins
    INTEGER             :: aa    ! aa and (aa+-1): bin number where the new concentration moves to
    REAL(dp)            :: r1,r2 ! Contraction fractions that move to bin a and a-1, respectively

    !Find the bin numbers (aa-1, a) where the content goes to
    aa = MINLOC(current_PSD%volume_fs-mix_PSD%volume_fs(ind),1, &
    mask = (current_PSD%volume_fs-mix_PSD%volume_fs(ind)) >= 0.0_dp)


    ! Move to largest bin or beyond
    IF (aa < 1) THEN
        n_b = new_PSD%nr_bins  ! to save some space
        ! New composition in new_PSD%nr_bins
        !new_PSD%composition_fs(n_b,:) = (new_PSD%composition_fs(n_b,:) * new_PSD%conc_fs(n_b) &
        !                        + mix_PSD%composition_fs(ind,:) * mix_PSD%conc_fs(ind)) &
        !                        / (new_PSD%conc_fs(n_b) + mix_PSD%conc_fs(ind))

        ! Particle concentration should be updated
        IF ( current_PSD%volume_fs(n_b) - mix_PSD%volume_fs(ind) < 0.d0 ) THEN  !the growth is into last bin and beyond => leave particles in last bin
          new_PSD%conc_fs(n_b) = new_PSD%conc_fs(n_b) + mix_PSD%conc_fs(ind)
        ELSE !the growth is into second last and last bin => redistribute particles and update composition in second-last bin
          ! Fraction of particles in size bin aa-1
          r1 = (current_PSD%volume_fs(n_b) - mix_PSD%volume_fs(ind)) &
          / (current_PSD%volume_fs(n_b) - current_PSD%volume_fs(n_b-1))
          ! Fraction of particles in size bin (aa)
          r2 = 1.0_dp - r1
          ! Determine new particle number concentrations:
          new_PSD%conc_fs(n_b-1) = new_PSD%conc_fs(n_b-1) + r1 * mix_PSD%conc_fs(ind)
          new_PSD%conc_fs(n_b) = new_PSD%conc_fs(n_b) + r2 * mix_PSD%conc_fs(ind)
          ! Don't update the composition any more in aa-1; would be possible, but what for?
        END IF

    ELSE IF (aa > 1 .and. aa <= current_PSD%nr_bins) THEN

        ! Fraction of particles in size bin aa-1
        r1 = (current_PSD%volume_fs(aa) - mix_PSD%volume_fs(ind)) &
        / (current_PSD%volume_fs(aa) - current_PSD%volume_fs(aa-1))
        ! Fraction of particles in size bin (aa)
        r2 = 1.0_dp - r1
        interm_PSD%conc_fs(aa) = r2 * mix_PSD%conc_fs(ind)
        interm_PSD%conc_fs(aa-1) = r1 * mix_PSD%conc_fs(ind)

        ! Determine new mass compositions: composition of fraction that goes to a-1
        interm_PSD%composition_fs(aa-1,:) = mix_PSD%composition_fs(ind,:) * current_PSD%volume_fs(aa-1) / mix_PSD%volume_fs(ind)

        ! new composition in aa-1
        if (abs(new_PSD%conc_fs(aa-1) + interm_PSD%conc_fs(aa-1)) > 0d0) THEN
        new_PSD%composition_fs(aa-1,:) = (new_PSD%composition_fs(aa-1,:) * new_PSD%conc_fs(aa-1) &
                                        + interm_PSD%composition_fs(aa-1,:) * interm_PSD%conc_fs(aa-1)) &
                                        / (new_PSD%conc_fs(aa-1) + interm_PSD%conc_fs(aa-1))
        ELSE
            new_PSD%composition_fs(aa-1,:) = 0d0
        END IF

        ! composition of fraction that goes to aa
        interm_PSD%composition_fs(aa,:) = mix_PSD%composition_fs(ind,:) * current_PSD%volume_fs(aa) / mix_PSD%volume_fs(ind)

        ! new composition in aa
        if (abs(new_PSD%conc_fs(aa) + interm_PSD%conc_fs(aa)) > 0d0) THEN
            new_PSD%composition_fs(aa,:) = (new_PSD%composition_fs(aa,:) * new_PSD%conc_fs(aa) &
                                        +  interm_PSD%composition_fs(aa,:) * interm_PSD%conc_fs(aa)) &
                                        /  (new_PSD%conc_fs(aa) + interm_PSD%conc_fs(aa))
        ELSE
            new_PSD%composition_fs(aa,:) = 0d0
        END IF

        ! Determine new particle number concentrations:
        new_PSD%conc_fs(aa-1) = new_PSD%conc_fs(aa-1) + interm_PSD%conc_fs(aa-1)
        new_PSD%conc_fs(aa) = new_PSD%conc_fs(aa) + interm_PSD%conc_fs(aa)

    ! The particles stay in the same bin without any changes (should not happen)
    ELSE
        new_PSD%conc_fs(ind) = new_PSD%conc_fs(ind) + mix_PSD%conc_fs(ind)
        ! print*, 'Trouble ahead: ',ind, aa, GTIME%sec, 'mix conc:', mix_PSD%conc_fs(ind)
    END IF
END SUBROUTINE bin_redistribute_fs


!=================================================================
! This subroutine a) re1ribtes the particles that are newly formed by coagulation into the bins or b) redistributes
! particles that changed their seize due to condensation/evaporation
! It also determines the new composition. Results are saved in new_PSD%conc_ma and new_PSD%composition_ma, new_PSD%diameter_ma
! MOVING Average, FIXED GRID METHOD!
! ==============================================================
SUBROUTINE bin_redistribute_ma(ind)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ind   ! The bin which is to be redistributed -> it generally has different diameter than in diameter_ma!
    INTEGER             :: n_b   ! Short for number of bins
    INTEGER             :: aa    ! aa and (aa+-1): bin number where the new concentration moves to
    REAL(dp)            :: dp_ind  !diameter of the mix_PSD that is to be redistributed [m]
    REAL(dp)            :: old_new_dp !the diameter of new_PSD before it considers changes [m]

    ! XXX added max to keep volumes values out
    dp_ind = (6.d0 * max(0d0,mix_PSD%volume_ma(ind)) / pi) ** (1.d0/3.d0)

    ! Find the bin numbers (aa-1, a) where the content goes to
    aa = MINLOC(current_PSD%grid_ma - dp_ind,1, &
         mask = (current_PSD%grid_ma- dp_ind) >= 0.0_dp) - 1    !find the bin the content of mix_PSD is in

    ! Move beyond largest bin : don't change diameter here
    IF (aa == -1) THEN
        n_b = new_PSD%nr_bins  ! to save some space
        ! New composition in new_PSD%nr_bins
        new_PSD%composition_ma(n_b,:) = (new_PSD%composition_ma(n_b,:) * new_PSD%conc_ma(n_b) &
                                    + mix_PSD%composition_ma(ind,:) * mix_PSD%conc_ma(ind)) &
                                    / (new_PSD%conc_ma(n_b) + mix_PSD%conc_ma(ind))
        ! Determine new particle concentration in largest bin
        new_PSD%conc_ma(new_PSD%nr_bins) = new_PSD%conc_ma(new_PSD%nr_bins) &
                                                     + mix_PSD%conc_ma(ind)
    ! Move stuff to bins 1 to nr_bins
    ELSE IF (aa >=  1 .and. aa <= current_PSD%nr_bins) THEN

        ! Diameter before changes:
        old_new_dp = (6.d0 * SUM(new_PSD%composition_ma(aa,:) / current_PSD%density_ma(:)) / pi)

        ! New composition in aa
        new_PSD%composition_ma(aa,:) = (new_PSD%composition_ma(aa,:) * new_PSD%conc_ma(aa) &
                                    + mix_PSD%composition_ma(ind,:) * mix_PSD%conc_ma(ind)) &
                                    / (new_PSD%conc_ma(aa) + mix_PSD%conc_ma(ind))

        ! CHECK FOR NEGATIVES
        if (minval(new_PSD%composition_ma(aa,:))<0d0.and.Use_speed) THEN
            call SET_ERROR (PRC%coa,'Getting negative composition in MA Mass_Number_Change')   !d_dpar(maxloc(abs(d_dpar))
            return
        ELSE
            old_new_dp = old_new_dp ** (1.d0/3.d0)
        END IF
            ! Determine new diameter in bin aa
            ! new diameter based on total volume of both aerosols of bin ii
            new_PSD%diameter_ma(aa) = ((  old_new_dp ** 3.d0 * new_PSD%conc_ma(aa) + &
                                    dp_ind ** 3.d0 * mix_PSD%conc_ma(ind)) / &
                                    (new_PSD%conc_ma(aa) + mix_PSD%conc_ma(ind))) ** (1.d0/3.d0)
            ! Determine new particle concentration in bin aa
            new_PSD%conc_ma(aa) = new_PSD%conc_ma(aa) + mix_PSD%conc_ma(ind)

    ! The particles shrink below the lower size limit -> change concentration but not composition or diameter (should not happen)
    ELSE
        new_PSD%conc_ma(1) = new_PSD%conc_ma(1) + mix_PSD%conc_ma(ind)
        new_PSD%diameter_ma(1) = new_PSD%grid_ma(1)
    END IF
END SUBROUTINE bin_redistribute_ma

! =====================================================================================================================
! Updates the particle number concentration with all PSD_styles
! Updates only those bins that are within the given size range [lower_limit,upper_limit]
SUBROUTINE send_conc(lower_limit, upper_limit,conc_fit)
    IMPLICIT NONE
    INTEGER :: ii
    REAL(dp) :: lower_limit ! Lower limit of the size range that is considered for replacing current concentration [m]
    REAL(dp) :: upper_limit         !upper limit of the size range that is considered for replacing current concentration [m]
    real(dp), INTENT(IN) :: conc_fit(:)
    IF (current_PSD%PSD_style == 1) THEN
        DO ii = 1, current_PSD%nr_bins
            IF (current_PSD%diameter_fs(ii) >= lower_limit .and. current_PSD%diameter_fs(ii) <= upper_limit) THEN
                current_PSD%conc_fs(ii) = conc_fit(ii)
            END IF
        END DO
    ELSE IF (current_PSD%PSD_style == 2) THEN
        DO ii = 1, current_PSD%nr_bins
            IF (current_PSD%diameter_ma(ii) >= lower_limit .and. current_PSD%diameter_ma(ii) <= upper_limit) THEN
                current_PSD%conc_ma(ii) = conc_fit(ii)
            END IF
        END DO
    END IF
END SUBROUTINE send_conc


! =====================================================================================================================
! Updates the particle composition with all PSD_styles
SUBROUTINE set_composition(ip, gdp, useold)
    IMPLICIT NONE
    INTEGER :: ip
    real(dp):: gdp, sum_org
    real(dp):: comp_pp(VAPOUR_PROP%n_condtot)
    logical :: useold

    IF (current_PSD%PSD_style == 1) THEN
        if (useold) THEN
            comp_pp = old_PSD%composition_fs(ip,:) * Na / VAPOUR_PROP%molar_mass * old_PSD%conc_fs(ip)
              if (current_PSD%volume_fs(ip) >0.0 .and. current_PSD%diameter_fs(ip) > 1D-9) then
                sum_org = sum(comp_pp,DIM=1)
                if (sum_org>0) THEN
                  VAPOUR_PROP%mfractions= comp_pp/sum_org
                END if
              END if
        end if

        current_PSD%composition_fs(ip,:) = VAPOUR_PROP%mfractions * current_PSD%volume_fs(ip) * VAPOUR_PROP%density

    ELSE IF (current_PSD%PSD_style == 2) THEN
        if (useold) THEN
            comp_pp = old_PSD%composition_ma(ip,:) * Na / VAPOUR_PROP%molar_mass * old_PSD%conc_fs(ip)

              if (current_PSD%volume_ma(ip) >0.0 .and. current_PSD%diameter_ma(ip) > 1D-9) then
                sum_org = sum(comp_pp,DIM=1)
                if (sum_org>0) THEN
                  VAPOUR_PROP%mfractions= comp_pp/sum_org
                END if
              END if
        end if

        current_PSD%composition_ma(ip,:) = VAPOUR_PROP%mfractions * (pi*gdp**3d0)/6d0 * VAPOUR_PROP%density
    END IF
END SUBROUTINE set_composition


! =====================================================================================================================
! Returns a composition array of size(nr_bins) that is independant of PSD_style
PURE FUNCTION get_composition(parts)
    IMPLICIT none
    REAL(dp) :: get_composition(current_PSD%nr_bins,n_cond_tot)
    type(PSD),INTENT(IN), optional :: parts
    type(PSD) :: d_PSD
    if (PRESENT(parts)) THEN
        d_PSD = parts
    ELSE
        d_PSD = current_PSD
    END IF

    IF (d_PSD%PSD_style == 1) THEN
      get_composition = d_PSD%composition_fs
  ELSE IF (d_PSD%PSD_style == 2) THEN
      get_composition = d_PSD%composition_ma
    END IF
END FUNCTION get_composition


! =====================================================================================================================
! Returns a diameter array of size(nr_bins) that is independant of PSD_style
PURE FUNCTION get_dp(parts)
    IMPLICIT none
    REAL(dp) :: get_dp(current_PSD%nr_bins)
    type(PSD),INTENT(IN), optional :: parts
    type(PSD) :: d_PSD
    if (PRESENT(parts)) THEN
        d_PSD = parts
    ELSE
        d_PSD = current_PSD
    END IF

    IF (d_PSD%PSD_style == 1) THEN
        get_dp = d_PSD%diameter_fs
    ELSE IF (d_PSD%PSD_style == 2) THEN
        get_dp = exp(log(d_PSD%grid_ma(:(size(d_PSD%grid_ma)-1))) + 0.5*log(d_PSD%grid_ma(2)/d_PSD%grid_ma(1)))
    END IF
END FUNCTION get_dp


! =====================================================================================================================
! Returns a volume array of size(nr_bins) that is independant of PSD_style
PURE FUNCTION get_volume(parts)
    IMPLICIT none
    REAL(dp) :: get_volume(current_PSD%nr_bins)
    type(PSD),INTENT(IN), optional :: parts
    type(PSD) :: d_PSD
    if (PRESENT(parts)) THEN
        d_PSD = parts
    ELSE
        d_PSD = current_PSD
    END IF

    IF (d_PSD%PSD_style == 1) THEN
        get_volume = d_PSD%volume_fs
    ELSE IF (d_PSD%PSD_style == 2) THEN
        get_volume = d_PSD%volume_ma
    END IF
END FUNCTION get_volume


! =====================================================================================================================
! Returns a mass array of size(nr_bins) that is independant of PSD_style
PURE FUNCTION get_mass(parts) result(out)
    IMPLICIT none
    REAL(dp) :: out(current_PSD%nr_bins)
    type(PSD),INTENT(IN), optional :: parts
    type(PSD) :: d_PSD
    if (PRESENT(parts)) THEN
        d_PSD = parts
    ELSE
        d_PSD = current_PSD
    END IF

    IF (d_PSD%PSD_style == 1) THEN
        out = d_PSD%volume_fs * d_PSD%particle_density_fs
    ELSE IF (d_PSD%PSD_style == 2) THEN
        out = d_PSD%volume_ma * d_PSD%particle_density_ma
    END IF
END FUNCTION get_mass


! =====================================================================================================================
! Returns a concentration array of size(nr_bins) that is independant of PSD_style
PURE FUNCTION get_conc(parts) result(out)
    IMPLICIT none
    REAL(dp) :: out(current_PSD%nr_bins)
    type(PSD),INTENT(IN), optional :: parts
    type(PSD) :: d_PSD
    if (PRESENT(parts)) THEN
        d_PSD = parts
    ELSE
        d_PSD = current_PSD
    END IF
    IF (d_PSD%PSD_style == 1) THEN
        out = d_PSD%conc_fs
    ELSE IF (d_PSD%PSD_style == 2) THEN
          out = d_PSD%conc_ma
    END IF
END FUNCTION get_conc


END MODULE PSD_scheme
