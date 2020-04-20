MODULE ParticleSizeDistribution
! Probably I can remove the do loop for nr_times: it will be outside this module
! When used as a module, read in is:
! nr_bins,PSD_style,nr_species,dp_range,nr_channels,process
! change vectors: dconc_coag,dmass,dconc_dep_mix,mix_ratio

  ! Question: should we have the density temperature dependent?
  USE second_precision, ONLY: dp
  USE constants
  USE INPUT

  IMPLICIT NONE

  ! START: variables that will be defined outside
  INTEGER :: nr_species_P  ! number of species that can go to the particle phase
  INTEGER :: nr_cond
  ! END variables that will be defined outside

  REAL(dp), ALLOCATABLE :: dconc_coag(:,:)  ! coagulation: collision number matrix: nr_bins * nr_bins
  REAL(dp), ALLOCATABLE :: dmass(:,:)       ! change in particle mass (due to condensation or mixing) (nr_bins,nr_species_P)
  REAL(dp), ALLOCATABLE :: dconc_dep_mix(:) ! change in particle concentration (e.g. reduced coagulatio or mixing) (nr_bins)
  REAL(dp), ALLOCATABLE :: mix_ratio(:)     ! gives the rate ratio: added volume over present volume per timestep
  REAL(dp), ALLOCATABLE :: conc_pp(:,:)     ! concentration of all species in the particle phase [#/m⁻³] (nr_bins,nr_species_P)
  CHARACTER(len=15) :: process              ! defines the process that passes information to subroutine Mass_Number_Change (coagulation, condensation, mixing)
  type(PSD) :: current_PSD                  ! Main PSD container. This variable stores the current timestep concentrations
  type(PSD) :: new_PSD, mix_PSD, interm_PSD ! Variables that store PSD values during the calculations



  CONTAINS

  SUBROUTINE initialize_PSD
    !This subroutine calls other subroutines from module PSD to:
    !  1) extract some infromation from input
    !  2) allocate the array sizes for PSD
    !  3) generate the PSD related arrays
    IMPLICIT NONE

    CALL PSD_get_input()      ! get some variables from input module needed for PSD representation
    CALL PSD_Allocate()       ! allocate the variables representing the PSD
    CALL GeneratePSDarrays()  ! Generate the PSD arrays for modelling

  END SUBROUTINE initialize_PSD


  ! ==================================================================================
  ! This subroutine picks inputs needed from input MODULE to initialize the PSD module
  ! ==================================================================================
  SUBROUTINE PSD_get_input()

    IMPLICIT NONE

    INTEGER :: nr_noncond = 0      !number of species on the particle that are non-volatile

    !SET mode of particle size distribution representation:
    !   PSD_style: options for particle size distribution representation
    !   1 ... fully stationary (working)
    !   2 ... reserved for fully stationary(coagulation)/full moving(growth)
    !   3 ... reserved for fully moving
    !   4 ... reserved for hybrid
    current_PSD%PSD_style = PSD_MODE

    !Set the number of bins for model array:
    current_PSD%nr_bins = n_bins_particle
    new_PSD%nr_bins = n_bins_particle
    !Set the range the particle array should feature [m]
    current_PSD%dp_range = (/min_particle_diam,max_particle_diam/)

    !set the number of species considered in the particle phase:
    !H2SO4 is added to the particle phase => +1
    !non condensables initially on the particle phase
    !vapour_number is the stuff that has a quatifiable vapour pressure > 0
    if (ALLOCATED(XTRAS)) nr_noncond = size(XTRAS)
    nr_species_P = VAPOUR_PROP%vapour_number + 1 + nr_noncond
    nr_cond      = VAPOUR_PROP%vapour_number

  END SUBROUTINE PSD_get_input


! ===========================================================================
! This subroutine allocates the particle size distribution relevant variables
! based on choice of representation type
! Further it might be used to minimizes the dimensions of varibales used in
! other PSD representation types to save computation time
! ===========================================================================
SUBROUTINE PSD_Allocate()
  IMPLICIT NONE

  ALLOCATE(conc_pp(current_PSD%nr_bins,nr_species_P))
  IF (current_PSD%PSD_style == 1) THEN
      ! FULLY STATIONARY representation !

    ALLOCATE(current_PSD%diameter_fs(current_PSD%nr_bins))
    ALLOCATE(current_PSD%volume_fs(current_PSD%nr_bins))
    ALLOCATE(current_PSD%conc_fs(current_PSD%nr_bins))
    ALLOCATE(current_PSD%composition_fs(current_PSD%nr_bins,nr_species_P))
    ALLOCATE(current_PSD%density_fs(nr_species_P))
    ALLOCATE(current_PSD%particle_density_fs(current_PSD%nr_bins))
    ALLOCATE(current_PSD%particle_mass_fs(current_PSD%nr_bins))

    ALLOCATE(dmass(current_PSD%nr_bins,nr_species_P))
    ALLOCATE(dconc_dep_mix(current_PSD%nr_bins))
    ALLOCATE(dconc_coag(current_PSD%nr_bins,current_PSD%nr_bins))
    ALLOCATE(mix_ratio(current_PSD%nr_bins))

    new_PSD = current_PSD
    mix_PSD = current_PSD
    interm_PSD = current_PSD

  ELSE
    PRINT*,'  This PSD representation is not yet defined:', current_PSD%PSD_style
    STOP
  END IF
  END SUBROUTINE PSD_Allocate


! Generates the particle related array content based on the choice of representation type
! Concentration is set to zero here -> it is set to input values in subroutine "GeneratePSDfromInput"
SUBROUTINE GeneratePSDarrays()
  IMPLICIT NONE
  INTEGER :: i  ! integer for incrementation
  REAL(dp) :: binw_ratio, binwl(2), factor_voldiam,factor_areadiam

  current_PSD%conc_fs = 0.d0
  conc_pp = 0.d0

  IF (current_PSD%PSD_style == 1) THEN
    !=====================================================================================================
    ! FULLY STATIONARY METHOD!

    ! Define diameter ratio from lower and upper limits
    binw_ratio = (log(current_PSD%dp_range(2))-log(current_PSD%dp_range(1)))/(current_PSD%nr_bins-1)
    ! Define diameter array as exponent of a geometric series
    current_PSD%diameter_fs =  exp([(i*binw_ratio + log(current_PSD%dp_range(1)), i=0,current_PSD%nr_bins-1)])
    current_PSD%diameter_fs(1) = current_PSD%dp_range(1)
    ! calculate first bin width
    binwl =  exp([((i+0.5)*binw_ratio + log(current_PSD%dp_range(1)), i=-1,0)])
    ! calculate the ratio of the representative area diameter to nominal diameter
    factor_areadiam = ( 1.d0/(3*(binwl(2)-binwl(1))) * (binwl(2)**3 - binwl(1)**3)  )**(1/2.d0) / current_PSD%diameter_fs(1)
    ! calculate the ratio of the representative volume diameter to nominal diameter
    factor_voldiam = ( 1.d0/(4*(binwl(2)-binwl(1))) * (binwl(2)**4 - binwl(1)**4)  )**(1/3.d0) / current_PSD%diameter_fs(1)

    if (factor_voldiam > 1.01 .or. factor_areadiam > 1.01) THEN
      print FMT_WARN1, 'You should increase particle bins: volume diameter to nominal fraction: ', factor_voldiam
      print FMT_WARN1, 'You should increase particle bins: area diameter to nominal fraction:   ', factor_areadiam
    end if
    ! Define volume array:
    current_PSD%volume_fs = 1D0/6D0 * pi * (current_PSD%diameter_fs*factor_voldiam)**3.d0   ! Single particle volume (m^3)
    ! particle_density needed for diffusion
    current_PSD%particle_density_fs = 1.4D3
    ! particle_density needed for diffusion
    current_PSD%particle_mass_fs = current_PSD%volume_fs * current_PSD%particle_density_fs
    ! partice phase density for condensables
    current_PSD%density_fs=VAPOUR_PROP%density

  ELSE
    !=====================================================================================================
    ! SOME OTHER METHOD

    print*,'choose other form of representation:'
    print*,current_PSD%PSD_style,'is not defined yet'
  END IF

END SUBROUTINE GeneratePSDarrays


! Generates y distribution as a function of diameter from input (dp_from, conc_from by fitting the property y (e.g. particle concentration)
! Note: monodisperse peaks are captured
! Important: fits tend to underestimate the concentration if not monodisperse:
! if dNdlogdp(i) > 0. and dNdlogdp(i +1 or -1) = 0. then the fit in between is 0!
! Note this routine is not tested for cases where input size resolution is higher than in model

SUBROUTINE GeneratePSDfromInput(dp_from,conc_from,conc_out)
  IMPLICIT NONE
  REAL(dp), INTENT(IN) :: dp_from(:), conc_from(:)  !arrays: 1) dp  2)other property as input for fitting
  REAL(dp), INTENT(INOUT) :: conc_out(:) ! Final output, fitted concentration array

  REAL(dp) :: ddNdlogdp   !difference in dNdlogdp between two channels
  REAL(dp) :: dp_diff(size(dp_from))  !parameter to capture difference between model diameter array and input diameter

  LOGICAL :: mono         ! if mono is true -> capture monodisperse peak, else: ignore (to avoid to capture it several times)
  INTEGER :: n_channels   ! Dimensions of the input PSD
  INTEGER :: channel      ! bin number which is close and smaller than to input diameter
  INTEGER :: j,k          ! integers for loops

  n_channels = size(dp_from)
  IF (current_PSD%PSD_style == 1) THEN
    ! ==========================================================================
    ! FULLY STATIONARY METHOD

    mono = .true.
    DO j = 1, current_PSD%nr_bins
      DO k = 1, n_channels
        dp_diff(k) = current_PSD%diameter_fs(j) - dp_from(k)
        IF (dp_diff(k) < 0.d0) dp_diff(k) = 1.d10
      END DO

      ! in case the the diameter is smaller or larger than the smallest or largest channel diameter => bin conc = 0.d0
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
        !Input peak ranges over at least 2 input channels -> fit required
          ddNdlogdp = conc_from(channel+1) - conc_from(channel)
          conc_out(j) = conc_from(channel) + (ddNdlogdp) * LOG10(current_PSD%diameter_fs(j) &
                          & / dp_from(channel)) / LOG10(dp_from(channel+1) / dp_from(channel))

          mono = .true.
        END IF
        IF (conc_from(channel) <= 0.d0) mono = .true.
      END IF

      IF (channel > 0) THEN
        IF (mono) THEN ! check for monodisperse peaks
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
              !In case there is a monodisperse peak at the end:
              mono = .false.
              conc_out(j) = conc_from(channel)

            END IF
          END IF
        END IF
      END IF

    END DO
    ! Convert from dNdlogdo to concentration:
    conc_out = conc_out * LOG10(current_PSD%diameter_fs(2)/current_PSD%diameter_fs(1))

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
! Subroutine input is change-array dmass (nr_species_P,nr_bins) for mass
! Subroutine input is change-array dconc_dep_mix(nr_bins) or dconc_coag(nr_bins,nr_bins) for number
! ===================================================================================================================
SUBROUTINE Mass_Number_Change(process)
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(IN)  :: process

  IF (current_PSD%PSD_style == 1) THEN
    ! FULLY STATIONARY METHOD =============================================

    new_PSD%composition_fs = current_PSD%composition_fs  !set the new particle composition to the current -> needed in some cases (e.g.: if there is no growth)
    IF (process == 'condensation') THEN  !Condensation
      new_PSD%conc_fs = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
      ! new_PSD%composition_fs = 0.d0
      CALL PSD_Change_Condensation()
    ELSE IF (process == 'coagulation') THEN !Coagulation
      new_PSD%conc_fs = current_PSD%conc_fs  ! set the initial value of new particle concentration to the current
      CALL PSD_Change_coagulation()
    ELSE IF (process == "mixing") THEN
      new_PSD%conc_fs = 0.d0  !Set some initial values for the number concentration (also needed for mass content)
      CALL PSD_Change_mixing()
    ELSE IF (process == 'deposition') THEN !Deposition simply reduces the particle number by the number that is passed via dconc_dep_mix - that's it
      new_PSD%conc_fs = current_PSD%conc_fs - dconc_dep_mix !reduce number of particles due to deposition; composition remains the same
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
SUBROUTINE PSD_Change_Condensation()
  IMPLICIT NONE
  INTEGER   :: i ! integer for looping

  IF (current_PSD%PSD_style == 1) THEN
    ! FULLY STATIONARY METHOD! ============================================
    DO i = 1, current_PSD%nr_bins
      if (current_PSD%conc_fs(i)>0d0) THEN
        mix_PSD%volume_fs(i) =  current_PSD%volume_fs(i) &
          + SUM(dmass(i,:) / current_PSD%density_fs(:) / current_PSD%conc_fs(i)) ! in m^3
        ELSE
          mix_PSD%volume_fs(i) = 0d0
        END IF


      mix_PSD%composition_fs(i,:) = dmass(i,:) + current_PSD%composition_fs(i,:)
      mix_PSD%conc_fs(i) = current_PSD%conc_fs(i)

      IF ( mix_PSD%conc_fs(i) > 0.d0 ) CALL bin_redistribute_fs(i)

      mix_PSD%conc_fs(i) = 0.d0
    END DO
  ELSE
    print*,'choose other form of representation:'
    print*,current_PSD%PSD_style,'is not defined yet'
  END IF
END SUBROUTINE PSD_Change_Condensation


! ====================================================================================================================
! This subroutine uses the dconc_coag array (nr_bins * nr_bins) which is determined by coagulation
! In the first row there is 1 entry at (1,1): concentration loss of smallest particles due to collision with each other
! In line n there are n entries(n,1:n): concentration loss of smallest particles due to collision with particle of size
! n, second smallest,... last entry: collisions of n with n
! Result is the new particle size distribution and composition: new_PSD%conc_fs and new_PSD%composition_fs
SUBROUTINE PSD_Change_coagulation
  IMPLICIT NONE
  INTEGER :: i,j  ! integer for incrementation

  IF (current_PSD%PSD_style == 1) THEN
    ! =========================================================
    ! FULLY STATIONARY METHOD!
    mix_PSD%conc_fs = 0.d0  !initial value is zero -> content is determined below

    ! Apply changes for all combinations of i and j
    DO i = 1, current_PSD%nr_bins-1 ! changed by Carlton. going to nr bins gives unrealistic number concentration in nr_bins
      DO j = 1, i
        ! Reduce the new particle concentration by the number of particles that are lost by coagulation in i and j
        ! -> they will be added later to the new bin
        IF (dconc_coag(i,j) > 0.d0) THEN

          ! reduce number in i
          new_PSD%conc_fs(i) = new_PSD%conc_fs(i) - dconc_coag(i,j)*GTIME%dt_aer

          ! reduce number in j (if i=j we have to reduce twice (which is done here) as 1 collision removes 2 particles)
          new_PSD%conc_fs(j) = new_PSD%conc_fs(j) - dconc_coag(i,j)*GTIME%dt_aer

          ! Determine new mass compositions: (total mass of collision products (i+j) + total mass already in mix bin i)
          ! / devided by sum of concentration (collisions + mix)
          ! composition of of collision result bins: i + j
          mix_PSD%composition_fs(i,:) = (current_PSD%composition_fs(i,:) + current_PSD%composition_fs(j,:))

          ! Update concentration in the mix_PSD
          mix_PSD%conc_fs(i) = dconc_coag(i,j)*GTIME%dt_aer

          ! Determine volume of the mixing aerosol
          mix_PSD%volume_fs(i) = current_PSD%volume_fs(i) + current_PSD%volume_fs(j)

          CALL bin_redistribute_fs(i)

        END IF
      END DO
    END DO

  ELSE
    print*,'choose other form of representation:'
    print*,current_PSD%PSD_style,'is not defined yet'
  END IF
END SUBROUTINE PSD_Change_coagulation






  SUBROUTINE PSD_Change_mixing
    !This subroutine uses the dconc_dep_mix and dmass vectors provided by mixing subroutine
    !dconc_dep_mix give the size dependent particle number that is added: (nr_bins)
    !dmass gives the composition of those particles: (nr_bins,nr_species_P)
    !Result is the new particle size distribution and composition: new_PSD%conc_fs and new_PSD%composition_fs
    !mix_ratio dtermines the ratio between present and added volume

    IMPLICIT NONE
    INTEGER :: i ! integer for incrementation

    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      mix_PSD%conc_fs = 0.d0  !initial value is zero -> content is determined below

      !apply changes for all combinations of i and j

      DO i = 1, current_PSD%nr_bins
        IF (dconc_dep_mix(i) > 0.d0) THEN
          !Set mass of mixing aerosol:
          mix_PSD%composition_fs(i,:) = dmass(i,:)  !composition of i
          !set concentration in the mix_PSD
          mix_PSD%conc_fs(i) = dconc_dep_mix(i) * mix_ratio(i)
          !Determine new composition and particle concentrations
          new_PSD%conc_fs(i) = (current_PSD%conc_fs(i) + mix_PSD%conc_fs(i) * mix_ratio(i)) / &
                               (1.d0 + mix_ratio(i))
          new_PSD%composition_fs(i,:) = (current_PSD%composition_fs(i,:) + mix_PSD%composition_fs(i,:) * mix_ratio) / &
                                        (1.d0 + mix_ratio)
          mix_PSD%conc_fs(i) = 0.d0  !reset value to zero -> changes have been applied
        ELSE
          new_PSD%conc_fs(i) = current_PSD%conc_fs(i) !nothing to mix with => new concentration is current concentration; composition is already new (see subroutine Mass_Number_Change)
        END IF
      END DO
    ELSE
      print*,'choose other form of representation:'
      print*,current_PSD%PSD_style,'is not defined yet'
    END IF
  END SUBROUTINE PSD_Change_mixing


  !=================================================================
  ! This subroutine a) redistribtes the particles that are newly formed by coagulation into the bins or b) redistributes
  ! particles that changed their seize due to condensation/evaporation
  ! It also determines the new composition. Results are saved in new_PSD%conc_fs and new_PSD%composition_fs
  ! FULLY STATIONARY METHOD!
  ! ==============================================================
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
    IF (aa == 0) THEN
      n_b = new_PSD%nr_bins  ! to save some space
      ! New composition in new_PSD%nr_bins
      new_PSD%composition_fs(n_b,:) = (new_PSD%composition_fs(n_b,:) * new_PSD%conc_fs(n_b) &
                                    + mix_PSD%composition_fs(ind,:) * mix_PSD%conc_fs(ind)) &
                                    / (new_PSD%conc_fs(n_b) + mix_PSD%conc_fs(ind))
      ! Determine new particle concentration in largest bin
      new_PSD%conc_fs(new_PSD%nr_bins) = new_PSD%conc_fs(new_PSD%nr_bins) &
                                                     + current_PSD%conc_fs(ind)
    ! Coagulation: (partly) leads to change in bin; condensation: growth or shrinkage
    ELSE IF (aa > 1 .and. aa <= current_PSD%nr_bins) THEN
      ! Fraction of particles in size bin aa-1
      r1 = (current_PSD%volume_fs(aa) - mix_PSD%volume_fs(ind)) / &
           (current_PSD%volume_fs(aa) - current_PSD%volume_fs(aa-1))
      ! Fraction of particles in size bin (aa)
      r2 = 1.0_dp - r1
      interm_PSD%conc_fs(aa) = r2 * mix_PSD%conc_fs(ind)
      interm_PSD%conc_fs(aa-1) = r1 * mix_PSD%conc_fs(ind)

      ! Determine new mass compositions: composition of fraction that goes to a-1
      interm_PSD%composition_fs(aa-1,:) = mix_PSD%composition_fs(ind,:) * current_PSD%volume_fs(aa-1) / mix_PSD%volume_fs(ind)
      ! new composition in aa-1
      new_PSD%composition_fs(aa-1,:) = (new_PSD%composition_fs(aa-1,:) * new_PSD%conc_fs(aa-1) &
                                     + interm_PSD%composition_fs(aa-1,:) * interm_PSD%conc_fs(aa-1)) &
                                     / (new_PSD%conc_fs(aa-1) + interm_PSD%conc_fs(aa-1))
      ! composition of fraction that goes to a
      interm_PSD%composition_fs(aa,:) = mix_PSD%composition_fs(ind,:) * current_PSD%volume_fs(aa) / mix_PSD%volume_fs(ind)
      ! new composition in a
      new_PSD%composition_fs(aa,:) = (new_PSD%composition_fs(aa,:) * new_PSD%conc_fs(aa) + &
                                    interm_PSD%composition_fs(aa,:) * interm_PSD%conc_fs(aa)) / &
                                   (new_PSD%conc_fs(aa) + interm_PSD%conc_fs(aa))

      ! Determine new particle number concentrations:
      new_PSD%conc_fs(aa-1) = new_PSD%conc_fs(aa-1) + interm_PSD%conc_fs(aa-1)
      new_PSD%conc_fs(aa) = new_PSD%conc_fs(aa) + interm_PSD%conc_fs(aa)

    ! The particles stay in the same bin without any changes (should not happen)
    ELSE
      new_PSD%conc_fs(ind) = new_PSD%conc_fs(ind) + mix_PSD%conc_fs(ind)
    END IF

    current_PSD%conc_fs(ind) = new_PSD%conc_fs(ind)
    current_PSD%composition_fs(ind,:)  = new_psd%composition_fs(ind,:)

    ! Currently we can't handle zero conc properly
    if (current_PSD%conc_fs(ind) < 1d-200) then
      print '("| ",t3, a,i0,a,es10.3,t100, "|")', 'zero or less particles in bin ', ind,': ', current_PSD%conc_fs(ind)
      current_PSD%conc_fs(ind)=1D-18
    end if

  END SUBROUTINE bin_redistribute_fs



END MODULE ParticleSizeDistribution
