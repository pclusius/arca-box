PROGRAM Supermodel

USE second_MAIN                         ! Main second file
USE second_PARAMETERS                   ! CH_NSPEC (originally NSPEC) and chemical indices, ind_xxxx, come from here
USE second_Precision,  ONLY : dp        ! KPP Numerical type
USE second_Monitor,    ONLY : SPC_NAMES ! Names of chemicals from KPP
USE Chemistry
USE constants
USE AUXILLARIES
USE INPUT
use OUTPUT
USE ACDC_NH3
USE ACDC_DMA
USE SOLVEBASES
USE aerosol_auxillaries
USE ParticleSizeDistribution
USE aerosol_dynamics

IMPLICIT NONE


! ==================================================================================================================
! VARIABLE DECLARATION. MOST OF THE GLOBAL VARIABLES ARE DEFINED IN INPUT.F90 and CONSTANTS.F90
! ==================================================================================================================

REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)    ! Array to hold input variables for the current timestep
REAL(DP), ALLOCATABLE :: CH_GAS(:)        ! Array to hold all chemistry compounds
REAL(dp), ALLOCATABLE :: conc_fit(:)      ! An array giving particle conc independant of PSD_style [m⁻³]

INTEGER               :: dmps_ln = 0      ! line number from where background particles are read from
INTEGER               :: dmps_sp_min = 0, dmps_sp_max = 0 ! Indices for dmps_special

! Variables related to chemistry module
REAL(dp) :: CH_RO2      ! RO2 concentration in [molecules / cm^3]
REAL(dp) :: CH_H2O      ! H20 concentration in [molecules / cm^3]
REAL(dp) :: CH_Beta     ! solar zenit angle
REAL(dp) :: EW          ! Water content in Pa
REAL(dp) :: ES          ! Saturation vapour pressure
! Variables related to aerosol module
INTEGER, DIMENSION(:), allocatable :: index_cond
REAL(dp), dimension(:), allocatable:: conc_vapour

! Transient variables
CHARACTER(:), allocatable:: RUN_OUTPUT_DIR ! Saves the output directory relative to this executable
CHARACTER(1000) :: inibuf       ! Buffer to save backp from the INITfile that was called
INTEGER         :: I,J          ! Loop indices
INTEGER         :: ioi          ! iostat variable
REAL(dp)        :: cpu1, cpu2   ! CPU time in seconds

! Temporary variables -> will be replaced
REAL(dp), ALLOCATABLE :: nominal_dp(:)    ! array with nominal diameters. Stays constant independent of PSD_style

! speed_up: factor for increasing integration time step for individual prosesses
! (1): Photolysis, (2): chemistry; (3):Nucleation; (4):Condensation; (5): Coagulation; (6): Deposition
INTEGER          :: speed_up(10) = 1
TYPE(error_type) :: error

! ==================================================================================================================

! Welcoming message
print'(a,t35,a)', achar(10),  '--~:| HLS BOX v.0.3 |:~--'//achar(10)

! ==================================================================================================================
! Declare most variables and read user input and options in input.f90
CALL READ_INPUT_DATA
! ==================================================================================================================

! ==================================================================================================================
IF (Chemistry_flag) THEN
    ! Check that the input exists in chemistry, or if not, print warning
    CALL CHECK_INPUT_AGAINST_KPP
    ! This only called once for KPP in the beginning
    CALL KPP_SetUp
ENDIF
! ==================================================================================================================

! ==================================================================================================================
! If particles are considered -> initialize a particle representation, set initial PSD and determine composition
IF (Aerosol_flag) THEN
    ! OMP is parallel processing and currently has only marginal effect
    if (USE_OPENMP) write(*,FMT_SUB) 'Using openmp'

    ! Initialize the Particle representation
    CALL INITIALIZE_PSD

    ! Initialize the nominal diameter vector
    ALLOCATE(nominal_dp(n_bins_particle))
    nominal_dp = get_dp()

    ! Initialize the dummy for dmps fitting
    ALLOCATE(conc_fit(n_bins_particle))

    if (use_dmps_special .and. use_dmps) THEN
        write(*, FMT_MSG) 'Using dmps_special:'

        ! both limits have to be far enough for the smallest and largest diameters
        if (dmps_lowband_upper_limit > nominal_dp(2)) THEN
            dmps_sp_min = minloc(abs(get_dp()-dmps_lowband_upper_limit),1)
            print FMT_SUB, 'Lower bins replaced from 1 to '//i2chr(dmps_sp_min)
        ELSE IF (dmps_lowband_upper_limit > 0d0) THEN
            print FMT_FAT0, 'dmps_lowband_upper_limit is is too small and will not be used '
            STOP
        END IF

        if (dmps_highband_lower_limit < nominal_dp(n_bins_particle-1) &
         .and. dmps_highband_lower_limit > nominal_dp(2)) THEN
            dmps_sp_max = minloc(abs(get_dp()-dmps_highband_lower_limit),1)
            print FMT_SUB, 'Upper bins replaced from '//i2chr(dmps_sp_max)//' to '//i2chr(n_bins_particle)
        ELSE IF (dmps_highband_lower_limit > 0d0) THEN
            print FMT_FAT0, 'dmps_highband_lower_limit is is too large and will not be used'
            STOP
        END IF

    END IF

    ! Here "Hash table" for condensing vapours is created
    ALLOCATE(index_cond(n_cond_org))
    ALLOCATE(conc_vapour(n_cond_tot))

    index_cond=0
    ! check how many species we have in common
    DO i = 1,size(VAPOUR_PROP%vapour_names)
        j = IndexFromName(VAPOUR_PROP%vapour_names(i),SPC_NAMES)
        IF (j>0) index_cond(i) = j
    END DO

    print FMT_LEND,
    END IF ! IF (Aerosol_flag)

! Allocate and initiate the timestep vector for interpolated input
ALLOCATE(TSTEP_CONC(N_VARS))
TSTEP_CONC = 0

! Allocate and initiate the timestep vector for chemical compounds
ALLOCATE(CH_GAS(size(SPC_NAMES)))
CH_GAS = 0


! ==================================================================================================================
! Prepare output files etc.
! ==================================================================================================================
write(*,*)
print FMT_HDR, 'INITIALIZING OUTPUT '

! All run output goes to RUN_OUTPUT_DIR directory. RUN_OUTPUT_DIR is allocated to correct length so we dont need to TRIM every time
ALLOCATE(CHARACTER(len=LEN(TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(INDEX)//'/'//TRIM(RUN_NAME))) :: RUN_OUTPUT_DIR)
RUN_OUTPUT_DIR = TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(INDEX)//'/'//TRIM(RUN_NAME)

if (Aerosol_flag) THEN

    ! Open text files to save particle size ditribution in easily accessible format
    ! Sumfile that uses same format as SMEAR sumfiles, (with exception of time resolution; here model save interval is used)
    ! Format is:
    ! 000------------------000------------bin_1 diameter---bin_2 diameter---bin_3 diameter . bin_n diameter
    ! Time 0 (s)--Total particle number---dN/log10(dDp)---dN/log10(dDp)-----dN/log10(dDp) . . dN/log10(dDp)
    ! Time 1 (s)--Total particle number---dN/log10(dDp)---dN/log10(dDp)-----dN/log10(dDp) . . dN/log10(dDp)
    OPEN(101,file=RUN_OUTPUT_DIR//"/particle_conc.sum",status='replace',action='write')
    WRITE(101,*) 0d0,0d0,get_dp()

    ! Datfile where PSD is in linear form. Format is:
    ! 000----------bin_1 diameter---bin_2 diameter---bin_3 diameter . bin_n diameter
    ! Time 0 (s)---------dN--------------dN----------------dN . . . . . . . dN
    ! Time 1 (s)---------dN--------------dN----------------dN . . . . . . . dN
    OPEN(104,file=RUN_OUTPUT_DIR//"/particle_conc.dat",status='replace',action='write')
    WRITE(104,*) 0d0,get_dp()

END IF

! Save backup from the initfile
OPEN(UNIT=504, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=ioi)
open(unit=334, file=RUN_OUTPUT_DIR//'/InitBackup.txt')
DO while (ioi == 0)
    read(504,'(a)', iostat = ioi) inibuf
    write(334,'(a)') TRIM(inibuf)
end do
close(504)
close(334)

open(unit=333, file=RUN_OUTPUT_DIR//'/error_output.txt')

! Open netCDF files
CALL OPEN_FILES( RUN_OUTPUT_DIR, Description, MODS, CH_GAS, VAPOUR_PROP)

! If wait_for was defined in user options, wait for a sec
CALL PAUSE_FOR_WHILE(wait_for)

if (TRIM(INITIALIZE_WITH) /= '') THEN
    if (current_PSD%psd_style==1) CALL INITIALIZE_WITH_LAST(CURRENT_PSD%composition_fs, CURRENT_PSD%conc_fs,CH_GAS)
    if (current_PSD%psd_style==2) CALL INITIALIZE_WITH_LAST(CURRENT_PSD%composition_ma, CURRENT_PSD%conc_ma,CH_GAS)
END IF

write(*,*) ''
print FMT_HDR, 'Beginning simulation'
if (start_time_s>0) then
    write(*, '(a)', advance='no') 'Starting simulation at: '
    print FMT_TIME, GTIME%hms
end if
call cpu_time(cpu1) ! For efficiency calculation
! =================================================================================================
DO WHILE (GTIME%SIM_TIME_S - GTIME%sec > -1d-12) ! MAIN LOOP STARTS HERE
! =================================================================================================

    ! =================================================================================================
    ! =================================================================================================
    ! PREPARE AND PRINT OUT SOME FUNDAMENTAL STUFF
    ! =================================================================================================
    ! =================================================================================================

    ! Print time header in a nice way, start with empty row
    if (GTIME%printnow) THEN
        print *, ! Print time
        print FMT_TIME, GTIME%hms
        call cpu_time(cpu2) ! To compare real time vs simulatoin time, timer is stopped in the beginning
    ! if --gui flag was used, print a dot and EOL so the STDOUT-reading in Python GUI will be smoother.
    ELSE
      if (ingui) print'(a)', '.'
    END IF

    ! GC_AIR_NOW, TEMPK, PRES and RH are calculated as global variables and are available everywhere.
    ! USE GC_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3
    GTEMPK = interp(timevec, CONC_MAT(:,inm_tempK)) .mod. MODS(inm_tempK)
    GPRES = interp(timevec, CONC_MAT(:,inm_PRES)) .mod. MODS(inm_PRES)
    GRH = interp(timevec, CONC_MAT(:,inm_RH)) .mod. MODS(inm_RH)
    GC_AIR_NOW = C_AIR_cc(GTEMPK, GPRES)

    ! Calculate Water vapour pressure and concentration
    CALL WATER(ES,EW,CH_H2O)

    ! Assign values to input variables. N_VARS will cycle through all variables that user can provide
    ! or tamper, and leave zero if no input or mod was provided
    DO I = 1, N_VARS
        TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I)) .mod. MODS(I)
    END DO
    ! =================================================================================================


    ! =================================================================================================
    ! =================================================================================================
    ! Chemistry
    ! =================================================================================================
    ! =================================================================================================
    IF (Chemistry_flag) THEN
        DO I = 1, N_VARS ! <-- N_VARS will cycle through all input variables
            IF (INDRELAY_CH(I)>0) THEN ! <-- this will pick those that were paired in CHECK_INPUT_AGAINST_KPP
                IF (I == inm_H2SO4 .and. .not. model_H2SO4) THEN
                    CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                ELSE IF (I /= inm_H2SO4) THEN
                    CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                END IF
            END IF
        END DO
        ! Solar Zenith angle. For this to properly work, lat, lon and Date need to be defined in INIT_FILE
        call BETA(CH_Beta)
        Call CHEMCALC(CH_GAS, GTIME%sec, (GTIME%sec + GTIME%dt_chm), GTEMPK, TSTEP_CONC(inm_swr), CH_Beta,  &
                    CH_H2O, GC_AIR_NOW, TSTEP_CONC(inm_CS), TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2)

        if (model_H2SO4) TSTEP_CONC(inm_H2SO4) = CH_GAS(ind_H2SO4)

    END IF ! IF (Chemistry_flag)
    ! =================================================================================================


    ! =================================================================================================
    ! =================================================================================================
    ! NUCLEATION
    ! =================================================================================================
    ! =================================================================================================
    IF (.not. error%error_state) THEN
        if (ACDC) THEN
            CALL ACDC_J(TSTEP_CONC, acdc_iterations)
            J_TOTAL = J_ACDC_DMA + J_ACDC_NH3 ! [particles/s/m^3]
            IF (.not. RESOLVE_BASE) J_TOTAL = J_TOTAL + TSTEP_CONC(inm_JIN) * 1d6 ! [particles/s/m^3]
        else
          ! Only use input format rate
          J_TOTAL = TSTEP_CONC(inm_JIN) * 1d6
        END if
    END if
    if (GTIME%savenow .and. RESOLVE_BASE) CALL Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)
    ! =================================================================================================


    ! =================================================================================================
    ! =================================================================================================
    ! AEROSOL PART STARTS HERE
    ! =================================================================================================
    ! =================================================================================================
    IF (Aerosol_flag.and.(.not. error%error_state)) THEN

        ! Read in background particles
        if (use_dmps .and. GTIME%min >= (dmps_ln*dmps_tres_min) .and. GTIME%min/60d0 < DMPS_read_in_time) THEN
            ! Print user info
            if (gtime%printnow) print FMT_SUB, 'Reading in background particles from line '//i2chr(dmps_ln)
            ! Fit measured particle size distribution to the current model PSD
            CALL GeneratePSDfromInput( par_data(1,2:),  par_data(min(dmps_ln+2, size(par_data, 1)),2:), conc_fit )
            ! Sets the concentration and composition to the added particles PSD based on the mfractions.
            ! Sumfile is typically in particles /cm^3, so make sure dmps_multi is correct (in NML_CUSTOM)
            CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit*dmps_multi)
            do i = 1, n_bins_particle
                CALL set_composition(i,nominal_dp(i))
            end do
            ! Next time next line is read, so increase index with 1
            dmps_ln = dmps_ln + 1
        END IF

        if (use_dmps .and. use_dmps_special .and. (GTIME%min >= (dmps_ln*dmps_tres_min)) .and. (GTIME%min/60d0 .ge. DMPS_read_in_time)) THEN
            ! Print user info
            if (gtime%printnow) print FMT_MSG, 'Reading in background particles partially'
            ! Fit the measured particle size distribution to the current model PSD
            CALL GeneratePSDfromInput( par_data(1,2:),  par_data(min(dmps_ln+2, size(par_data, 1)),2:), conc_fit )
            ! Sumfile is typically in particles /cm^3, so make sure dmps_multi is correct
            conc_fit = conc_fit * dmps_multi

            ! if small particles are replaced with input
            if (dmps_sp_min>0) THEN
                ! Print user info
                if (gtime%printnow) print FMT_SUB, 'Lowband upper limit: '//i2chr(dmps_sp_min)
                ! Sets the concentration and composition to the partially added particles PSD based on the mfractions.
                CALL send_conc(current_PSD%dp_range(1),nominal_dp(dmps_sp_min),conc_fit)
                do i = 1, dmps_sp_min
                    CALL set_composition(i,nominal_dp(i))
                end do
            END IF

            ! If large particles are replaced with input
            if (dmps_sp_max>0) THEN
                ! Print user info
                if (gtime%printnow) print FMT_SUB, 'Highband lower limit: '//i2chr(dmps_sp_max)
                ! Sets the concentration and composition to the partially added particles PSD based on the mfractions.
                CALL send_conc(nominal_dp(dmps_sp_max),current_PSD%dp_range(2),conc_fit)
                do i = dmps_sp_max, n_bins_particle
                    CALL set_composition(i,nominal_dp(i))
                end do
            END IF

            ! Next time next line from measurement etc. is read
            dmps_ln = dmps_ln + 1

            END IF


        ! ADD NUCLEATED PARTICLES TO PSD, IN THE 1st BIN.  This using the mixing method, where two particle distros
        ! are mixed together. NOTE that this step is always done if Aerosol_flag is on, even if ACDC_flag is off. In
        ! case no NPF is wanted, turn off ACDC and make sure NUC_RATE_IN is 0
        dmass = 0d0
        dmass(1,VAPOUR_PROP%ind_HOA) = nominal_dp(1)**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_HOA)
        dconc_dep_mix = 0d0
        dconc_dep_mix(1) = J_TOTAL*GTIME%dt_aer
        ! Negative mixing ratio makes this nucleation
        mix_ratio = -1d0
        CALL Mass_Number_Change('mixing')
        ! Update current_psd
        current_PSD = new_PSD


        ! ..........................................................................................................
        ! CONDENSATION
        if (Condensation) THEN
            ! Pick the condensables from chemistry and change units from #/cm^3 to #/m^3
            conc_vapour = 0d0
            do i = 1, n_cond_org
              if (index_cond(i) /= 0) conc_vapour(i) =  CH_GAS(index_cond(i))*1D6 ! mol/m3
            end do
            ! Poor sulfuric acid always needs special treatment
            conc_vapour(n_cond_tot) = CH_GAS(ind_H2SO4)*1d6

            ! Update vapour pressures
            CALL Calculate_SaturationVapourConcentration(VAPOUR_PROP, GTEMPK, VP_MULTI)
            ! Solve mass flux
            dmass = 0d0
            CALL Condensation_apc(VAPOUR_PROP,conc_vapour,dmass)
            ! Distribute mass
            CALL Mass_Number_Change('condensation')
            ! Update current_psd
            current_PSD = new_PSD

            ! Update vapour concentrations to chemistry
            do i = 1, n_cond_org
                if (index_cond(i) /= 0) CH_GAS(index_cond(i)) = conc_vapour(i) *1D-6
            end do
            ! Again sulfuric acid always needs special treatment
            CH_GAS(ind_H2SO4) = conc_vapour(n_cond_tot)*1d-6

          end if ! end of condensation


        ! ..........................................................................................................
        ! COAGULATION
        if (Coagulation) then
            ! Solve particle coagulation
            Call Coagulation_routine(dconc_coag)
            ! Distribute mass
            Call Mass_Number_Change('coagulation')
            ! Update PSD with new concentrations
            current_PSD = new_PSD
        end if ! end of coagulation

    end if
    ! End Aerosol =====================================================================================


    ! SPEED handling, currently not implimented
    IF (error%error_state) THEN
        CALL error_handling(error, speed_up,GTIME%dt,GTIME%sec)
        error%error_state = .false.
    ! If there was no error during the actual timestep
    ELSE
        ! Write printouts to screen
        if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
        ! Save to netcdf
        if (GTIME%savenow) THEN
            CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, VAPOUR_PROP, CURRENT_PSD)
            ! PSD is also saved to txt
            if (Aerosol_flag) THEN
                WRITE(101,*) GTIME%sec, sum(get_conc()*1d-6), get_conc()*1d-6 / LOG10(bin_ratio)
                WRITE(104,*) GTIME%sec, get_conc()*1d-6
            END IF
        END IF

        ! Add main timestep to GTIME'
        GTIME = ADD(GTIME)

    END IF ! SPEED handling

END DO	! Main loop ends
! ==================================================================================================================
! ==================================================================================================================

if (Aerosol_flag) THEN
    CLOSE(101)
    CLOSE(104)
END IF

CALL PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

! Close output file netcdf
CALL CLOSE_FILES()
CALL FINISH


CONTAINS

! =================================================================================================
! Error handling for aerosol dynamics: In case any error appears during the calculation of
! dynamics, this subroutine is called to handle the error. This means adjusting the timestep (either
! via speed_up or dt_box) for the process the error showed up
! Note that this only works if we calculate "future" changes not what happend in the last time step!
! INPUTS:  a) errors; b) speed_up; dt_box
! Outputs: a) updated speed_up or dt_box
! =================================================================================================
SUBROUTINE error_handling(error,speed_up,dt_box,sim_time_sec)

    IMPLICIT none
    INTEGER ::  i    , &
        speed_up(10)
    REAL(DP) :: dt_box  ,&
        sim_time_sec
    type(error_type) :: error

    !Write some error information to file
    write(333,*) 'Precision error at simulation time [s]', sim_time_sec
    write(333,*) '  ERROR TYPE: ', trim(error%error_specification)
    !reduce the speed_up factor or, if necessary, the integration time step:
    IF (speed_up(error%error_process) > 1) THEN
        speed_up(error%error_process) = speed_up(error%error_process) / 2
        write(333,*) '  => reduce speed_up:',speed_up(error%error_process)*2, '->',speed_up(error%error_process)
    ELSE
        dt_box = dt_box / 2.d0
        write(333,*) '  => reduce dt [s]:', dt_box
        DO i = 1,10
            IF (i /= error%error_process) speed_up(i) = speed_up(i) * 2
        END DO
    END IF

END SUBROUTINE error_handling


! =================================================================================================
!
! ACDC Nucleation. See that all values here are in SI-UNITS: CUBIC METERS, KELVINS AND PASCALS.
! Written by Tinja Olenius
! .................................................................................................
! Input for get_acdc_J:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_base:            base (ammonia) concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
! TempK:             Temperature in Kelvins
! IPR:               Ion production rate in ion pairs per second [1/m3/s]. 3d6 is a good guestimate
! dt:                Main time step [s]
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_NH3:        Particle formation rate due to ammonia [1/s/m3]. Sum of J_NH3_BY_IONS
! acdc_cluster_diam: Outgrowing cluster diameter [m]
! J_NH3_BY_IONS:     Particle formation rate by ions (neutral, negative and positive) [1/s/m3].
!                    The outgrowing cluster typically has 5 to 6 H2SO4 in it
! .................................................................................................
! Input for get_acdc_D:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_dma:             DMA concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
! TempK:             Temperature in Kelvins
! dt:                Main time step [s]
! time:              model time; used only for outputting information in cluster
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_DMA:        Particle formation rate due to DMA [1/s/m3]
! acdc_cluster_diam: Outgrowing cluster diameter [m]. The cluster typically has 5 to 6 H2SO4 in it
!
! =================================================================================================

SUBROUTINE ACDC_J(C,iters)
    implicit none
    REAL(dp) :: c_org = 0d0
    INTEGER :: i
    INTEGER, INTENT(IN)  :: iters
    REAL(dp), intent(in) :: C(:)
    LOGICAL, save        :: first_time = .true., ss_handle = .true.
    REAL(dp)             :: H2SO4=0,NH3=0,DMA=0,IPR=0 ! these are created to make the unit conversion

    ! NUCLEATION BY S-ACID AND NH3 - NOTE: in ACDC, ingoing concentrations are assumed to be in 1/m3!!
    IF (inm_H2SO4 /= 0) H2SO4 = C(inm_H2SO4)*1d6
    IF (inm_NH3   /= 0) NH3   = C(inm_NH3)*1d6
    IF (inm_DMA   /= 0) DMA   = C(inm_DMA)*1d6
    IF (inm_IPR   /= 0) IPR   = C(inm_IPR)*1d6
    ! Speed up program by ignoring nucleation when there is none
    if ((NH3 > 1d12 .and. H2SO4>1d9) .or. (.not. skip_acdc)) THEN
        ! Idea of iteration is to not necessarily reach steady state but still get a "better" estimate of instantenous formation rate
        Do i=1,iters
            CALL get_acdc_J(H2SO4,NH3,c_org,C(inm_CS),C(inm_TEMPK),IPR,GTIME,&
                ss_handle,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
        end do
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
    END IF
    ! Speed up program by ignoring nucleation when there is none
    if ((DMA > 1d6 .and. H2SO4>1d9) .or. (.not. skip_acdc)) THEN
        ! Idea of iteration is to not necessarily reach steady state but still get a "better" estimate of instantenous formation rate
        Do i=1,iters
            CALL get_acdc_D(H2SO4,DMA,c_org,C(inm_CS),C(inm_TEMPK),GTIME,ss_handle,J_ACDC_DMA,acdc_cluster_diam)
        end do
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'DMA IGNORED'
    END IF

    ! The first time ACDC is run, it is run to steady state, after that the user supplied option is used
    if (first_time) THEN
        ss_handle = ACDC_solve_ss
        first_time = .false.
    end if
END SUBROUTINE ACDC_J  ! END ACDC Nucleation



! ================================================================================================
! PRINT_KEY_INFORMATION for the user. Can be called for example upon GTIME%printnow with
! if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
! ================================================================================================
SUBROUTINE PRINT_KEY_INFORMATION(C)
    IMPLICIT NONE
    REAL(dp), intent(in) :: C(:)
    print FMT10_2CVU,'ACID C: ', C(inm_H2SO4), ' [1/cm3]','sum(An)/A1',clusteracid,' []'
    ! print FMT10_2CVU,'APINE C: ', C(IndexFromName('APINENE'))
    print FMT10_3CVU,'Temp:', C(inm_TempK), ' [K]','Pressure: ', C(inm_pres), ' [Pa]', 'Air_conc', C_AIR_cc(C(inm_TempK), C(inm_pres)), ' [1/cm3]'
    IF (inm_NH3   /= 0) print FMT10_3CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]','sum(Nn)/N1',clusterbase,' []'
    IF (inm_DMA   /= 0) print FMT10_3CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]', 'J-total', J_TOTAL*1e-6,' [1/cm3]'
    print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', C(inm_CS) , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
    if ((GTIME%sec)>0) print '("| ",a,i0,a,i0.2,a,i0,a,i0.2,t65,a,f6.2,t100,"|")', 'Elapsed time (m:s) ',int(cpu2 - cpu1)/60,':',modulo(int(cpu2 - cpu1),60) ,' Est. time to finish (m:s) ',&
                            int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec))/60,':', MODULO(int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec)),60),&
                            'Realtime/Modeltime: ', (GTIME%sec-start_time_s)/(cpu2 - cpu1)

    print FMT_LEND,
END SUBROUTINE PRINT_KEY_INFORMATION


! =================================================================================================
! If user had opted a save or print interval which does not concide with last timestep, we print
! and save the values at last timesteps
! =================================================================================================
SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY
    implicit none
    GTIME = ADD(GTIME, -GTIME%dt)
    print*,''
    print FMT_HDR, 'Main loop ended'

    if (.not. GTIME%printnow) THEN
        print FMT_MSG, 'Values from the last timestep:'
        print FMT_TIME, GTIME%hms
        CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
    END IF
    if (.not. GTIME%savenow)  CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, VAPOUR_PROP, CURRENT_PSD)
END SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY


! ================================================================================================
! Subroutine checks if the user is sending in some variables that are not in the chemistry. In case of
! inorganics, a warning is issued if it is not in the chemistry, but if the input is an organic and
! not found in the chemistry, the program will give an error and halt.
! ================================================================================================
SUBROUTINE CHECK_INPUT_AGAINST_KPP
    implicit none
    integer :: i,j, check
    print FMT_HDR, 'Checking against KPP for chemicals'
    do i=1,N_VARS
        check = 0
        IF (MODS(I)%col > 0 .or. MODS(I)%MODE > 0 .or. ABS(MODS(I)%SHIFT) > 1d-100) THEN
            DO j=1,size(SPC_NAMES)
                IF (MODS(i)%NAME == TRIM(SPC_NAMES(j))) THEN
                    check = 1
                    print FMT_MSG, 'Found '//TRIM(SPC_NAMES(j))//' from chemistry'
                    ! store the 'key->value map' for input and chemistry
                    INDRELAY_CH(I) = J
                    IF ((I==inm_H2SO4) .and. model_H2SO4 .and. Chemistry_flag) print FMT_SUB, 'Replacing HSO4 input with modelled chemistry'
                    exit
                end if
            END DO
            IF (check == 0 .and. I>LENV) THEN
                print FMT_FAT0, 'You are using an (organic?) compound which does not exist in chemistry: '//TRIM(MODS(i)%NAME)//' '
                IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                IF (ABS(MODS(I)%SHIFT) > 1d-100) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
                print FMT_MSG, 'Good bye.'
                STOP
            END IF
            IF (check == 0 .and. I<LENV .and. &
                ((I==inm_H2SO4) .or. (I==inm_NH3) .or. (I==inm_DMA) .or. (I==inm_SO2) &
                .or. (I==inm_NO) .or. (I==inm_NO2) .or. (I==inm_CO) .or. (I==inm_H2) .or. (I==inm_O3))  ) THEN
                print FMT_WARN0, 'A compound which does not exist in chemistry but could be used elsewhere: '//TRIM(MODS(i)%NAME)//' '
                IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                IF (ABS(MODS(I)%SHIFT) > 1d-100) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
            END IF
        END IF
    END DO
    print FMT_LEND,

END SUBROUTINE CHECK_INPUT_AGAINST_KPP




! ============================================================================================================
! Calculate water saturation vapour pressure (ES), actual vapour pressure (EW) (in Pa) and concentration (CW)
! ============================================================================================================
SUBROUTINE WATER(ES,EW,CW)
    IMPLICIT NONE
    REAL(dp), INTENT(INOUT) :: ES, EW, CW
    REAL(dp)                :: TEMPC
    TEMPC = GTEMPK-273.15d0

    ! Saturation vapour pressure over liquid water; using parametrisation from ???, a0,a1...a6 are in constants.f90
    ES = (a0 + a1 * TEMPC**1 + a2 * TEMPC**2 + a3 * TEMPC**3    &
             + a4 * TEMPC**4 + a5 * TEMPC**5 + a6 * TEMPC**6)* 100

    ! Water vapour pressure
    EW  = GRH * ES / 100d0

    ! Water vapour concentration in molecules per cm3
    CW = EW/GPRES * GC_AIR_NOW

END SUBROUTINE WATER


! ================================================================================================
! Calculate Beta parameter (solar angle above horizon), parameterization by Henderson-Sellers
! ================================================================================================
SUBROUTINE BETA(CH_Beta)
    IMPLICIT NONE
    REAL(dp), INTENT(INOUT) :: CH_Beta
    REAL(dp) :: D_zsd,HRANG,ZSR

    D_zsd  = 0.409d0 * sin(2d0 * pi * (GTIME%JD-79.8) / 365.24)
    HRANG  = pi * (GTIME%sec / 43200. - 1.)
    ZSR    = sin(lat * pi / 180.) * sin(D_zsd) + cos(lat * pi / 180.) * cos(D_zsd) * cos(HRANG)
    CH_Beta   = 90d0 - acos(ZSR) * 180d0/PI
    IF (CH_Beta < 0d0) CH_Beta = 0d0
END SUBROUTINE BETA


! ================================================================================================
! Give the user time to read through the input stuff and continue upon Enter, or wait for amount
! defined in Wait_for (in INITFILE)
! ================================================================================================
SUBROUTINE PAUSE_FOR_WHILE(for)
    INTEGER :: for, rof
    CHARACTER(3) :: secs
    rof = for
    if (for == -1) THEN
        do while(rof == for)
            print *,
            print '(a,20(" "),a)', ACHAR(10),'Press any key to start the run or Ctrl-C to abort'
            read(*,*)
            print FMT_LEND,
            for = 0
        end do
    ELSE IF (for > 0) THEN
        write(secs,'(i0)') for
        print FMT_LEND,
        write(*, FMT_MSG) 'Starting in '//TRIM(secs)//' seconds'
        write(*, '("    ")', advance = 'no')
        do while (for>0)
            CALL SLEEP(1)
            for = for - 1
            write(*, '("..", i0)', advance='no') for
            if (MODULO(rof-for, 20) == 0) write(*,'(a,"    ")', advance='no') ACHAR(10)
        END DO
        write(*,*)
        print FMT_LEND,
    END IF
END SUBROUTINE PAUSE_FOR_WHILE

SUBROUTINE FINISH
    IMPLICIT NONE

    write(*,*)
    write(*,'(a,1(" "),a)', advance='no') 'SIMULATION HAS ENDED. '
    write(*, '(a)') 'SO LONG!'
    write(*,*)
    write(*,*)
    write(333,*) '---------- SIMULATION REACHED END SUCCESFULLY ------------'
    close(333)

END SUBROUTINE FINISH

END PROGRAM SUPERMODEL
