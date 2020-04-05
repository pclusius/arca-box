MODULE ParticleSizeDistribution
! Probably I can remove the do loop for nr_times: it will be outside this module
! When used as a module, read in is:
! nr_bins,PSD_style,nr_species,dp_range,nr_channels,process
! change vectors: dconc_coag,dmass,dconc_dep_mix,mix_ratio
! PSD_in
! USE SECOND_PRECISION,  ONLY : dp
! Question: should we have the density temperature dependent?

  USE second_precision, ONLY: dp
  USE constants
  USE INPUT


  IMPLICIT NONE

! Private
!
! Public ::
! START: variables that will be defined outside
! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)
! REAL(dp) :: pi = 2D0*ASIN(1D0)
  INTEGER :: nr_time_in  !the line of the input array (i.e.: specifies the time) which is used to generate the model PSD
  INTEGER :: nr_times      !number of points in time for PSD_in
  INTEGER :: nr_channels   !number of diameter channels for PSD_in
  INTEGER :: nr_species_P  !number of species that can go to the particle phase
  INTEGER :: nr_cond
  INTEGER :: nr_bins

  REAL(dp), ALLOCATABLE :: PSD_in(:,:)
  !END variables that will be defined outside

  REAL(dp), ALLOCATABLE :: dconc_coag(:,:)   !coagulation: collision number matrix: nr_bins * nr_bins
  REAL(dp), ALLOCATABLE :: dmass(:,:)        !change in particle mass (due to condensation or mixing) (nr_bins,nr_species_P)
  REAL(dp), ALLOCATABLE :: dconc_dep_mix(:)  !change in particle concentration (e.g. reduced coagulatio or mixing) (nr_bins)
  REAL(dp), ALLOCATABLE :: mix_ratio(:)      !gives the rate ratio: added volume over present volume per timestep
  REAL(dp), ALLOCATABLE :: dummy_property(:) !has dimension of diameter array and can be used for various things (e.g.: fitting)

  CHARACTER(len=15) :: process !defines the process that passes information to subroutine Mass_Number_Change (coagulation, condensation, mixing)

  real(dp), allocatable :: conc_pp(:,:)        !concentration of all species in the particle phase [#/m⁻³] (nr_bins,nr_species_P)
  real(dp), allocatable :: par_conc(:)

  type(PSD) :: current_PSD
  type(PSD) :: new_PSD
  type(PSD) :: old_PSD
  type(PSD) :: mix_PSD
  type(PSD) :: interm_PSD
  type(generic_PSD) :: par_dynamics
  real(dp), allocatable :: dia_particle(:), vol_particle(:)



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

    INTEGER :: nr_noncond      !number of species on the particle that are non-volatile

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
    nr_noncond = size(XTRAS)  !this will have to be derived from input: XTRAS -> don't know how at the moment
    nr_species_P = vapours%vapour_number + 1 !+ nr_noncond
    nr_cond      = vapours%vapour_number
    !input particle size distribution:
    nr_times = size(par_data(2:,1))
    nr_channels = size(par_data(1,2:))
    nr_bins = current_PSD%nr_bins
    !write(*,*) current_PSD%PSD_style, current_PSD%nr_bins, current_PSD%dp_range
    !write(*,*) nr_species_P, nr_noncond, nr_times, nr_channels
    !PAUSE

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
    ALLOCATE(dia_particle(current_PSD%nr_bins + 2))
    ALLOCATE(vol_particle(current_PSD%nr_bins + 2))

    par_dynamics%nr_bins   = current_PSD%nr_bins
    par_dynamics%dp_range  = current_PSD%dp_range

    ALLOCATE(par_dynamics%diameter(par_dynamics%nr_bins))
    ALLOCATE(par_dynamics%volume(par_dynamics%nr_bins))
    ALLOCATE(par_dynamics%conc(par_dynamics%nr_bins))
    ALLOCATE(par_dynamics%composition(par_dynamics%nr_bins,nr_species_P))
    ALLOCATE(par_dynamics%density(nr_species_P))
    ALLOCATE(par_dynamics%particle_density(par_dynamics%nr_bins))
    ALLOCATE(par_dynamics%particle_mass(par_dynamics%nr_bins))

    IF (current_PSD%PSD_style == 1) THEN
      ! FULLY STATIONARY representation !
      ALLOCATE(current_PSD%diameter_fs(current_PSD%nr_bins))
      ALLOCATE(current_PSD%volume_fs(current_PSD%nr_bins))
      ALLOCATE(current_PSD%conc_fs(current_PSD%nr_bins))
      ALLOCATE(par_conc(current_PSD%nr_bins))
      ALLOCATE(current_PSD%composition_fs(current_PSD%nr_bins,nr_species_P))

      ALLOCATE(current_PSD%density_fs(nr_species_P))
      ALLOCATE(current_PSD%particle_density_fs(current_PSD%nr_bins))
      ALLOCATE(current_PSD%particle_mass_fs(current_PSD%nr_bins))
      ALLOCATE(dmass(current_PSD%nr_bins,nr_species_P))
      ALLOCATE(dconc_dep_mix(current_PSD%nr_bins))
      ALLOCATE(dconc_coag(current_PSD%nr_bins,current_PSD%nr_bins))
      ALLOCATE(mix_ratio(current_PSD%nr_bins))
      ALLOCATE(dummy_property(current_PSD%nr_bins))

      ALLOCATE(new_PSD%diameter_fs(current_PSD%nr_bins))
      ALLOCATE(new_PSD%volume_fs(current_PSD%nr_bins))
      ALLOCATE(new_PSD%conc_fs(current_PSD%nr_bins))
      ALLOCATE(new_PSD%composition_fs(current_PSD%nr_bins,nr_species_P))
      ALLOCATE(new_PSD%density_fs(nr_species_P))
      ALLOCATE(new_PSD%particle_density_fs(current_PSD%nr_bins))
      ALLOCATE(new_PSD%particle_mass_fs(current_PSD%nr_bins))

      ALLOCATE(mix_PSD%diameter_fs(current_PSD%nr_bins))
      ALLOCATE(mix_PSD%volume_fs(current_PSD%nr_bins))
      ALLOCATE(mix_PSD%conc_fs(current_PSD%nr_bins))
      ALLOCATE(mix_PSD%composition_fs(current_PSD%nr_bins,nr_species_P))
      ALLOCATE(mix_PSD%density_fs(nr_species_P))
      ALLOCATE(mix_PSD%particle_density_fs(current_PSD%nr_bins))
      ALLOCATE(mix_PSD%particle_mass_fs(current_PSD%nr_bins))

      ALLOCATE(interm_PSD%diameter_fs(current_PSD%nr_bins))
      ALLOCATE(interm_PSD%volume_fs(current_PSD%nr_bins))
      ALLOCATE(interm_PSD%conc_fs(current_PSD%nr_bins))
      ALLOCATE(interm_PSD%composition_fs(current_PSD%nr_bins,nr_species_P))
      ALLOCATE(interm_PSD%density_fs(nr_species_P))
      ALLOCATE(interm_PSD%particle_density_fs(current_PSD%nr_bins))
      ALLOCATE(interm_PSD%particle_mass_fs(current_PSD%nr_bins))

  ELSE
      PRINT*,'  This PSD representation is not yet defined:', current_PSD%PSD_style
      STOP
    END IF
  END SUBROUTINE PSD_Allocate






  SUBROUTINE GeneratePSDarrays()
    !Generates the particle related array content based on the choice of representation type
    !Concentration is set to zero here -> it is set to input values in subroutine "GeneratePSDfromInput"
    IMPLICIT NONE
    INTEGER :: &
               i  !some integer for incrementation

    INTEGER, PARAMETER  ::  modes=5
    REAL(dp), DIMENSION(modes) :: s, N_modes, dm

    REAL(dp), DIMENSION(current_PSD%nr_bins) :: d_p, d_g, dlogDp
    current_PSD%conc_fs = 0.d0
    conc_pp = 0.d0

    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! for testing

     if (distribution1) then
        ! Particle diameters between 2D-9 and 2.5D-6 m:
       current_PSD%diameter_fs(1)=2D-9
       DO i=2,current_PSD%nr_bins
        current_PSD%diameter_fs(i)=current_PSD%diameter_fs(i-1)*(2.5D-6/current_PSD%diameter_fs(1))**(1D0/(current_PSD%nr_bins-1))
       END DO
       current_PSD%conc_fs = 1D0 ! Assume an initial particle number concentration of 1 [# m-3]
       ! current_PSD%conc_fs(current_PSD%nr_bins) = 2d8
       where((abs(current_PSD%diameter_fs-2D-7)-MINVAL(abs(current_PSD%diameter_fs-2D-7)))<1D-20)  current_PSD%conc_fs=2D8  ! add 200 [# cm-3] to 200 nm sized accumulation mode particles
       current_PSD%volume_fs = 1D0/6D0 * pi * current_PSD%diameter_fs**3.d0   ! Single particle volume (m^3)

     elseif(distribution2) then ! use some arbitary particle size distribution
       dia_particle(1)=2D-9 !1.95D-9
       vol_particle(1)=(pi*dia_particle(1)**3)/6
       DO i=2,current_PSD%nr_bins+1
          vol_particle(i)=vol_particle(i-1)*1.2!048!9!1.199
       END DO
       dia_particle = (vol_particle*6./pi)**(1D0/3D0)

       d_p = (dia_particle(1:current_PSD%nr_bins) + dia_particle(2:current_PSD%nr_bins+1))/2D0    ! Arithmetic mean diameter in each
       d_g = (dia_particle(1:current_PSD%nr_bins)*dia_particle(2:current_PSD%nr_bins+1))**0.5;    ! Geometric mean diameter in each size bin
       dlogDp=log10(dia_particle(2:current_PSD%nr_bins+1))-log10(dia_particle(1:current_PSD%nr_bins))

       dm = (/ 1.9317115D-8, 3.2812671D-8, 4.8989795D-8, 1.2D-7, 3.2719494D-7 /) ! Mode diameter size distribution (m)
       s = (/ 1.1833333, 1.2666667, 1.4333333, 1.5, 1.2666667 /)  ! Standard deviation for each mode
       N_modes = (/ 1D3, 1D3, 1D3, 1D9, 1D1 /) ! Number concentration modes (#/m^3)
       ! dm = (/ 3.2812671D-8, 4.8989795D-8, 3.2719494D-7 /) ! Mode diameter size distribution (m)
       ! s = (/ 1.2666667, 1.4333333, 1.2666667 /)  ! Standard deviation for each mode
       ! N_modes = (/ 1D3, 1D9, 1D3 /) ! Number concentration modes (#/m^3)

       current_PSD%diameter_fs = d_p !arithmetic mean diameter
       CALL dNdlogDp(d_g,d_p,dlogDp,modes,s,N_modes,dm,current_PSD%conc_fs) ! Initial size dist
       ! where (current_PSD%diameter_fs<1D0) current_PSD%diameter_fs=1D0
       ! write(*,*) 'l 285 sum(current_PSD%conc_fs)', sum(current_PSD%conc_fs)
       current_PSD%volume_fs = 1D0/6D0 * pi * current_PSD%diameter_fs**3.d0
      else
       current_PSD%diameter_fs(1) = current_PSD%dp_range(1)
       do i = 2,current_PSD%nr_bins
        current_PSD%diameter_fs(i) = current_PSD%diameter_fs(i-1) * &
        (current_PSD%dp_range(2) / current_PSD%dp_range(1)) ** &
        (1.D0/(current_PSD%nr_bins-1))
       end do
       current_PSD%volume_fs = 1D0/6D0 * pi * current_PSD%diameter_fs**3.d0   ! Single particle volume (m^3)
     end if
      !Define diameter array:

      !! particle_density .. neede for diffusion
      current_PSD%particle_density_fs = 1.4D3
      current_PSD%particle_mass_fs = 1D0/6D0 * pi * current_PSD%diameter_fs**3 * current_PSD%particle_density_fs
      par_conc=current_PSD%conc_fs
      current_PSD%density_fs=vapours%density

    ELSE
      print*,'choose other form of representation:'
      print*,current_PSD%PSD_style,'is not defined yet'
    END IF
    !PRINT*,current_PSD%diameter_fs
    !PRINT*,current_PSD%volume_fs
  END SUBROUTINE GeneratePSDarrays





  SUBROUTINE GeneratePSDfromInput(dp_fit,y_fit)
    !Generates y distribution as a function of diameter from input (dp_fit, y_fit -> input_dp_y) by fitting the property y (e.g. particle concentration)
    !Note: monodisperse peaks are captured
    !Important: fits tend to underestimate the concentration if not monodisperse:
    !if dNdlogdp(i) > 0. and dNdlogdp(i +1 or -1) = 0. then the fit in between is 0!
    !Note this routine is not tested for cases where input size resolution is higher than in model
    IMPLICIT NONE
    LOGICAL :: mono         ! if mono is true -> capture monodisperse peak, else: ignore (to avoid to capture it several times)
    LOGICAL :: warn_range   ! if false, checks whether the range for simulation and input fit -> otherwise we might loose data from input distribution

    ! INTEGER :: i            ! position in PSD_in -> 2 ... initial; 3 ... after 1 time interval (1 is the diamter array)
    INTEGER :: nr_channels  ! Dimensions of the input PSD
    INTEGER :: j,k          ! some integers for loops
    INTEGER :: channel      ! bin number which is close and smaller than to input diameter
    ! INTEGER :: nr_fitted    ! number of bins for fitted array (at the moment -> always nr_bins)

    REAL(dp) :: ddNdlogdp   !difference in dNdlogdp between two channels
    REAL(dp),ALLOCATABLE :: &
                            dp_diff(:), &  !parameter to capture difference between model diameter array and input diameter
                            input_dp_y(:,:)    !input array consisting of diameter in first row and property y (e.g.: particle concentration) for semilog fitting
    REAL(dp), ALLOCATABLE, INTENT(IN) :: dp_fit(:), y_fit(:)  !arrays: 1) dp  2)other property as input for fitting
!  REAL(dp), ALLOCATABLE :: y_fitted(:)  !arrays: 1) dp  2)other property as input for fitting

  !Allocate arrays dealing with inputs
    nr_channels = size(dp_fit)
    ALLOCATE(input_dp_y(2,nr_channels))
    ALLOCATE(dp_diff(nr_channels))
    input_dp_y(1,:) = dp_fit
    input_dp_y(2,:) = y_fit
    ! PRINT*, 'dp in: ',input_dp_y(1,:)
    ! PRINT*, 'conc in: ',input_dp_y(2,:)
    ! PRINT*,nr_channels
  !Allocate arrays containing fitted information:
  !  ALLOCATE(y_fitted(nr_fitted))
    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      mono = .true.
      warn_range = .true.
      DO j = 1, current_PSD%nr_bins
        DO k = 1, nr_channels
          dp_diff(k) = current_PSD%diameter_fs(j) - input_dp_y(1,k)
          IF (dp_diff(k) < 0.d0) dp_diff(k) = 1.d10!!
          !PRINT*,dp_diff(k-1)
        END DO
        IF ( warn_range .eqv. .false.) THEN
          IF (current_PSD%diameter_fs(current_PSD%nr_bins) < input_dp_y(1,nr_channels)) THEN
            PRINT*, 'Attention: the chosen diameter range for simulation does not cover the whole input range!'
            PRINT*, 'dp simulation max:', current_PSD%diameter_fs(current_PSD%nr_bins)
            PRINT*, 'max input diameter:', input_dp_y(1,nr_channels)
            PRINT*,current_PSD%diameter_fs(current_PSD%nr_bins) < input_dp_y(1,nr_channels),warn_range
            warn_range = .true.
            ! PAUSE ! Lukas, this is a deprecated feature
          END IF
        END IF
        !PRINT*, current_PSD%diameter_fs(j),input_dp_y(1,1:),dp_diff

        !in case the the diameter is smaller or larger than the smalles channel diameter => bin conc = 0.d0
        !else: assign lower channel for fitting
        IF (minval(dp_diff) > 1.d9) THEN
          !PRINT*,'not'
          dummy_property(j) = 0.d0
          channel = 0
          mono = .true.
        ELSE
          channel = minloc((dp_diff), DIM=1)
          !PRINT*,'  channel   :',channel
        END IF

        ddNdlogdp = 0.d0
        IF (channel > 0 .and. channel < nr_channels) THEN
          IF (channel < nr_channels .and. input_dp_y(2,channel) > 0.d0 .and. input_dp_y(2,channel+1) > 0.d0) THEN
          !Input peak ranges over at least 2 input channels -> fit required
            ddNdlogdp = input_dp_y(2,channel+1) - input_dp_y(2,channel)
            dummy_property(j) = input_dp_y(2,channel) + (ddNdlogdp) * &
            !linear fit (on linear diameter scale):
            !  (current_PSD%diameter_fs(j) / input_dp_y(1,channel)) / &
            !  (input_dp_y(1,channel+1) / input_dp_y(1,channel))
            !linear fit on log diameter scale:
            LOG10(current_PSD%diameter_fs(j) / input_dp_y(1,channel)) / &
            LOG10(input_dp_y(1,channel+1) / input_dp_y(1,channel))
            mono = .true.
          END IF
          IF (input_dp_y(2,channel) <= 0.d0) mono = .true.
        END IF
        IF (channel > 0) THEN
          IF (mono) THEN !check for monodisperse peaks
            IF (channel == 1) THEN
              IF (input_dp_y(2,channel) > 0.d0 .and. input_dp_y(2,channel+1) <= 0.d0) THEN
                !In case there is a monodisperse peak in channel 1:
                dummy_property(j) = input_dp_y(2,channel)
                mono = .false.
              END IF
            END IF
            IF (channel > 1 .and. channel < nr_channels ) THEN
              IF (input_dp_y(2,channel-1) <= 0.d0 .and. input_dp_y(2,channel) > 0.d0 .and. input_dp_y(2,channel+1) <= 0.d0) THEN
                !In case there is a monodisperse peak somewhere in the middle:
                mono = .false.
                dummy_property(j) = input_dp_y(2,channel)
              END IF
            END IF
            IF (channel == nr_channels) THEN
              IF (input_dp_y(2,channel-1) <= 0.d0 .and. input_dp_y(2,channel) > 0.d0) THEN
                !In case there is a monodisperse peak at the end:
                mono = .false.
                dummy_property(j) = input_dp_y(2,channel)
              END IF
            END IF
          END IF
        END IF

        !Convert from dNdlogdo to concentration:
        !dummy_property(j) = dummy_property(j) * &
        !LOG10(current_PSD%diameter_fs(2)/current_PSD%diameter_fs(1))
      END DO
    ELSE
      print*,'choose other form of representation:'
      print*,current_PSD%PSD_style,'is not defined yet'
    END IF
    DEALLOCATE(input_dp_y)
    DEALLOCATE(dp_diff)
  END SUBROUTINE GeneratePSDfromInput






  SUBROUTINE Mass_Number_Change(process)
    !Applies changes to the particle size distribution amd composition according to input-change-vector:
    ! 1) dconc_coag matrix for coagulation;
    ! 2) dconc_dep_mix and dmass for mixing
    ! 3) dmass for condensation
    !Based on current and absolute change in single particle mass/number [units similar to composition and concentraion] the new distribution is determined
    !Subroutine input is change-array dmass (nr_species_P,nr_bins) for mass
    !Subroutine input is change-array dconc_dep_mix(nr_bins) or dconc_coag(nr_bins,nr_bins) for number
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)  :: process
    ! Real(dp), dimension(par_dynamics%nr_bins, nr_species_P), intent(in) :: dmass
    ! INTEGER :: i     !some integer for looping

    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! PRINT*,process
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






  SUBROUTINE PSD_Change_condensation()
    !Determines the changes in the PSD due to condensation based on dmass
    !Condensation is continuous process: all particles grow at the same rate
    !ATTENTION: The composition in the largest bin is wrong if growth happens via a==0!!
    !           This might cause errors -> find some solution at some point
    IMPLICIT NONE
    INTEGER   :: i     !some integer for looping
    ! REAL(dp)  :: r1,r2    !contration fractions that move to bin a and a-1, respectively
    ! real(dp) :: test(current_PSD%nr_bins)

    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! write(*,*) 'in PSD_Change_Condensation dmass', sum(dmass)
      !Find new volume in case of condensation
      DO i = 1, current_PSD%nr_bins
        mix_PSD%volume_fs(i) =  current_PSD%volume_fs(i) + &
                                 SUM(dmass(i,:) / current_PSD%density_fs(:) / current_PSD%conc_fs(i)) ! in m^3


        mix_PSD%composition_fs(i,:) = dmass(i,:) + current_PSD%composition_fs(i,:)
        mix_PSD%conc_fs(i) = current_PSD%conc_fs(i)
        IF ( mix_PSD%conc_fs(i) > 0.d0 ) CALL bin_redistribute_fs(i)
        mix_PSD%conc_fs(i) = 0.d0
       END DO
    ELSE
      print*,'choose other form of representation:'
      print*,current_PSD%PSD_style,'is not defined yet'
    END IF
  END SUBROUTINE PSD_Change_condensation






  SUBROUTINE PSD_Change_coagulation
    !This subroutine uses the dconc_coag array (nr_bins * nr_bins) which is determined by coagulation
    !In the first row there is 1 entry at (1,1): concentration loss of smallest particles due to collision with each other
    !In line n there are n entries(n,1:n): concentration loss of smallest particles due to collision with particle of size n, second smallest,... last entry: collisions of n with n
    !Result is the new particle size distribution and composition: new_PSD%conc_fs and new_PSD%composition_fs

    IMPLICIT NONE
    INTEGER :: i,j  !some integer for incrementation

    IF (current_PSD%PSD_style == 1) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !FULLY STATIONARY METHOD!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      mix_PSD%conc_fs = 0.d0  !initial value is zero -> content is determined below
      ! new_PSD%conc_fs = current_PSD%conc_fs  !the initial concentration is like the old -> changes are applied below
      ! mix_PSD%composition_fs =0.0d0
      ! new_PSD%composition_fs = 0.0d0

      !apply changes for all combinations of i and j
      DO i = 1, current_PSD%nr_bins-1 ! changed by Carlton. going to nr bins gives unrealistic number concentration in nr_bins
        DO j = 1, i
          IF (dconc_coag(i,j) > 0.d0) THEN
            !PRINT*,'coag',i,j
            ! PRINT*, mix_PSD%composition_fs(i,:), mix_PSD%conc_fs(i)
            !Reduce the new particle concentration by the number of particles that are lost by coagulation in i and j -> they will be added later to the new bin
            new_PSD%conc_fs(i) = new_PSD%conc_fs(i) - dconc_coag(i,j)*MODELTIME%dt_aero  !reduce number in i
            new_PSD%conc_fs(j) = new_PSD%conc_fs(j) - dconc_coag(i,j)*MODELTIME%dt_aero  !reduce number in j (if i=j we have to reduce twice (which is done here) as 1 collision removes 2 particles)
            !Determine new mass compositions: (total mass of collision products (i+j) + total mass already in mix bin i) / devided by sum of concentration (collisions + mix)
            mix_PSD%composition_fs(i,:) = (current_PSD%composition_fs(i,:) + current_PSD%composition_fs(j,:))  !composition of of collision result bins: i + j
            !update concentration in the mix_PSD
            mix_PSD%conc_fs(i) = dconc_coag(i,j)*MODELTIME%dt_aero

            !Determine volume of the mixing aerosol

            mix_PSD%volume_fs(i) = current_PSD%volume_fs(i) +current_PSD%volume_fs(j)!SUM(mix_PSD%composition_fs(i,:) / current_PSD%density_fs(:) / current_PSD%conc_fs(i))   !Determine the volume of particles in the mix distribution
            ! original one
            ! mix_PSD%volume_fs(i) = SUM(mix_PSD%composition_fs(i,:) / current_PSD%density_fs(:) / current_PSD%conc_fs(i))
            !PRINT*,current_PSD%composition_fs(i,:) + current_PSD%composition_fs(j,:), dconc_coag(i,j)
            !PRINT*, mix_PSD%composition_fs(i,:), mix_PSD%conc_fs(i)
            !PAUSE
            ! print*, mix_PSD%volume_fs(i), '',SUM(mix_PSD%composition_fs(i,:) / current_PSD%density_fs(:) / current_PSD%conc_fs(i))
            CALL bin_redistribute_fs(i)
            mix_PSD%conc_fs(i) = 0.d0  !reset value to zero -> changes have been applied
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
    INTEGER :: i ! some integer for incrementation

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




  SUBROUTINE bin_redistribute_fs(i)
    !!!!!!!!!!!!!!!!!!!!!!!!!
    !FULLY STATIONARY METHOD!
    !!!!!!!!!!!!!!!!!!!!!!!!!
    !This subroutine a) redistribtes the particles that are newly formed by coagulatio into the bins or b) redistributes
    !particles that changed their seize due to condensation/evaporation
    !It also determines the new composition
    !Results are saved in new_PSD%conc_fs and new_PSD%composition_fs
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::   i  !the bin which is to be redistributed -> it generally has different diameter than in diameter_fs!
    INTEGER :: &
               j,k,  &!some integer for incrementation
               a     !a and (a+-1): bin number where the new concentration moves to
    REAL(dp) :: &
               r1,r2   !contration fractions that move to bin a and a-1, respectively


    !Find the bin numbers (a-1, a) where the content goes to
    a = MINLOC(current_PSD%volume_fs-mix_PSD%volume_fs(i),1, &
        mask = (current_PSD%volume_fs-mix_PSD%volume_fs(i)) >= 0.0_dp)
    IF (a == 0) THEN !move to largest bin or beyond
      j = new_PSD%nr_bins  !to save some space
      new_PSD%composition_fs(j,:) = (new_PSD%composition_fs(j,:) * new_PSD%conc_fs(j) + &    !new composition in new_PSD%nr_bins
                                    mix_PSD%composition_fs(i,:) * mix_PSD%conc_fs(i)) / &
                                    (new_PSD%conc_fs(j) + mix_PSD%conc_fs(i))
      !Determine new particle concentration in largest bin
      new_PSD%conc_fs(new_PSD%nr_bins) = new_PSD%conc_fs(new_PSD%nr_bins) &
                                                     + current_PSD%conc_fs(i)
    ELSE IF (a > 1 .and. a <= current_PSD%nr_bins) THEN !coagulation: (partly) leads to change in bin; condensation: growth or shrinkage
      r1 = (current_PSD%volume_fs(a) - mix_PSD%volume_fs(i)) / &
           (current_PSD%volume_fs(a) - current_PSD%volume_fs(a-1))  ! Fraction of particles in size bin a-1
      r2 = 1.0_dp - r1  ! Fraction of particles in size bin (a)
      interm_PSD%conc_fs(a) = r2 * mix_PSD%conc_fs(i)
      interm_PSD%conc_fs(a-1) = r1 * mix_PSD%conc_fs(i)
      !Determine new mass compositions:
      interm_PSD%composition_fs(a-1,:) = mix_PSD%composition_fs(i,:) * current_PSD%volume_fs(a-1) / mix_PSD%volume_fs(i)   !composition of fraction that goes to a-1
      new_PSD%composition_fs(a-1,:) = (new_PSD%composition_fs(a-1,:) * new_PSD%conc_fs(a-1) + &    !new composition in a-1
                                    interm_PSD%composition_fs(a-1,:) * interm_PSD%conc_fs(a-1)) / &
                                   (new_PSD%conc_fs(a-1) + interm_PSD%conc_fs(a-1))
      interm_PSD%composition_fs(a,:) = mix_PSD%composition_fs(i,:) * current_PSD%volume_fs(a) / mix_PSD%volume_fs(i) !composition of fraction that goes to a
      new_PSD%composition_fs(a,:) = (new_PSD%composition_fs(a,:) * new_PSD%conc_fs(a) + &   !new composition in a
                                    interm_PSD%composition_fs(a,:) * interm_PSD%conc_fs(a)) / &
                                   (new_PSD%conc_fs(a) + interm_PSD%conc_fs(a))
      !Determine new particle number concentrations:
      new_PSD%conc_fs(a-1) = new_PSD%conc_fs(a-1) + interm_PSD%conc_fs(a-1)
      new_PSD%conc_fs(a) = new_PSD%conc_fs(a) + interm_PSD%conc_fs(a)
    ELSE !The particles stay in the same bin without any changes (should not happen)
      new_PSD%conc_fs(i) = new_PSD%conc_fs(i) + mix_PSD%conc_fs(i)
    END IF

    current_PSD%conc_fs(i) = new_PSD%conc_fs(i)
    current_PSD%composition_fs(i,:)  = new_psd%composition_fs(i,:)


  ! if particle concentations  < 0d0
    if (current_PSD%conc_fs(i) <= 0D0) then
      current_PSD%conc_fs(i)=1D0
    end if

  END SUBROUTINE bin_redistribute_fs

  SUBROUTINE dNdlogDp(d_g,d_p,dlogDp,modes,s,N_modes,dm,N_bins)
     INTEGER ::  ii
     INTEGER, INTENT(IN) :: modes
     REAL(dp), DIMENSION(modes), INTENT(IN) :: s,N_modes,dm
     REAL(dp), DIMENSION(modes,current_PSD%nr_bins)   :: dNdlogDp_modes
     REAL(dp), DIMENSION(current_PSD%nr_bins), INTENT(IN) :: d_g,d_p,dlogDp
     REAL(dp), DIMENSION(current_PSD%nr_bins)  :: Vp
     REAL(dp), DIMENSION(current_PSD%nr_bins), INTENT(OUT) :: N_bins

     DO ii = 1,modes
           dNdlogDp_modes(ii,:) = ( N_modes(ii) / 2.506628 / LOG10(s(ii)))*&
             EXP(-((LOG10(d_g) - log10(dm(ii)))**2 ) / (2*LOG10(s(ii))**2 ))
             ! write(*,*) 'L704 psd', sum(EXP(-((LOG10(d_g) - log10(dm(ii)))**2 ) / (2*LOG10(s(ii))**2 )) ),sum((LOG10(d_g))), log10(dm(ii)), &
             ! sum(((LOG10(d_g) - log10(dm(ii)))**2))
     END DO

     N_bins=SUM(dNdlogDp_modes, DIM=1)*dlogDp+1D-3


     ! Vp=d_p**3.0*pi/6.0
     ! V_bins=N_bins*Vp


  END SUBROUTINE dNdlogDp

END MODULE ParticleSizeDistribution
