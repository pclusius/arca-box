PROGRAM ARCA_main

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
USE PSD_scheme
USE aerosol_dynamics
USE custom_functions

IMPLICIT NONE

! ==================================================================================================================
! Note about file unit handles: Numbers between 100-499 are not used anywhere in the model, so use these in you need
! to open additional files. Generally 1-99 can be used temporarily in KPP, 600-series is used in ARCA_main.f90,
! 700-series is used in output.f90, 800-series is used in input.f90, 900-series in Chemistry.f90
! ==================================================================================================================

! ==================================================================================================================
! VARIABLE DECLARATION. MOST OF THE GLOBAL VARIABLES ARE DEFINED IN INPUT.F90 and CONSTANTS.F90
! ==================================================================================================================

REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)    ! Array to hold input variables for the current timestep
REAL(DP), ALLOCATABLE :: CH_GAS(:), CH_GAS_old(:)        ! Array to hold all chemistry compounds; old: to restore in case of an error related to timestep handling
REAL(dp), ALLOCATABLE :: conc_fit(:)      ! An array giving particle conc independant of PSD_style [m⁻³]
REAL(dp), ALLOCATABLE :: losses_fit(:),losses_fit0(:),losses_fit1(:), intrp_losses(:) ! An array giving particle conc independant of PSD_style [m-3]
REAL(dp), ALLOCATABLE :: save_measured(:)      ! An aray for saving measurements

INTEGER               :: dmps_ln = 0      ! line number from where background particles are read from
INTEGER               :: dmps_sp_min = 0, dmps_sp_max = 0 ! Indices for dmps_special

! Variables related to chemistry module
REAL(dp) :: CH_RO2,CH_RO2_old      ! RO2 concentration in [molecules / cm^3]
REAL(dp) :: CH_H2O      ! H20 concentration in [molecules / cm^3]
REAL(dp) :: CH_Beta     ! solar zenit angle
REAL(dp) :: EW          ! Water content in Pa
REAL(dp) :: ES          ! Saturation vapour pressure

! Variables related to aerosol module
REAL(dp), dimension(:), allocatable:: conc_vapour

! Transient variables
CHARACTER(:), allocatable:: RUN_OUTPUT_DIR ! Saves the output directory relative to this executable
CHARACTER(1000) :: inibuf       ! Buffer to save backp from the INITfile that was called
INTEGER         :: I,II,J,JJ    ! Loop indices
INTEGER         :: ioi          ! iostat variable
REAL(dp)        :: cpu1, cpu2   ! CPU time in seconds

! Temporary variables -> will be replaced
REAL(dp), ALLOCATABLE :: nominal_dp(:)    ! array with nominal diameters. Stays constant independent of PSD_style

! speed_up: factor for increasing integration time step for individual prosesses
! (1): chemistry; (2):Condensation; (3): Coagulation; (4): Deposition
TYPE(error_type)        :: PRCION
LOGICAL                 :: Handbrake_on = .false.
INTEGER                 :: speed_up(size(PRCION%pr_name,1)) = 1
INTEGER                 :: n_of_Rounds = 0
INTEGER                 :: optis_in_use(2) = [2,3]
REAL(dp), ALLOCATABLE   :: d_dpar(:)  ! array reporting relative changes to the diameter array within a single timestep
REAL(dp), ALLOCATABLE   :: d_npar(:)  ! array reporting relative changes to the particle number array within a single timestep
REAL(dp), ALLOCATABLE   :: d_vap(:)   ! array reporting relative changes to the vapour concentration array within a single timestep
REAL(dp)                :: DT_0

REAL(dp) :: Av
REAL(dp) :: Au
REAL(dp) :: Ad
REAL(dp) :: V_chamber
REAL(dp) :: E_field = 0d0
REAL(dp), ALLOCATABLE :: corgwallTeflon(:)


! Welcoming message
print'(a,t35,a)', achar(10),  '--~:| ARCA BOX MODEL 0.9 |:~--'//achar(10)

! ==================================================================================================================
! Declare most variables and read user input and options in input.f90
CALL READ_INPUT_DATA
! ==================================================================================================================

Av = CHAMBER_CIRCUMFENCE*CHAMBER_HEIGHT
Au = Chamber_floor_area
Ad = Chamber_floor_area
V_chamber = Chamber_floor_area*CHAMBER_HEIGHT


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
    ! if (USE_OPENMP) write(*,FMT_SUB) 'Using openmp'

    ! Initialize the Particle representation
    CALL INITIALIZE_PSD

    ! Initialize the nominal diameter vector
    ALLOCATE(nominal_dp(n_bins_par))
    nominal_dp = get_dp()

    ! Initialize the dummy for dmps fitting
    ALLOCATE(conc_fit(n_bins_par))
    save_measured = conc_fit
    losses_fit     = conc_fit
    losses_fit0    = conc_fit
    losses_fit1    = conc_fit
    ALLOCATE(intrp_losses(size(PAR_LOSSES%sections)))

    ! Allocate the change vectors for integration timestep control
    ALLOCATE(d_dpar(n_bins_par))
    ALLOCATE(d_npar(n_bins_par))

    if (use_dmps_partial .and. use_dmps) THEN
        write(*, FMT_MSG) 'Using dmps_special:'

        ! both limits have to be far enough for the smallest and largest diameters
        if (dmps_lowband_upper_limit > nominal_dp(2)) THEN
            dmps_sp_min = minloc(abs(get_dp()-dmps_lowband_upper_limit),1)
            print FMT_SUB, 'Lower bins replaced from 1 to '//i2chr(dmps_sp_min)
        ELSE IF (dmps_lowband_upper_limit > 0d0) THEN
            print FMT_FAT0, 'dmps_lowband_upper_limit is is too small and will not be used '
            STOP
        END IF

        if (dmps_highband_lower_limit < nominal_dp(n_bins_par-1) &
         .and. dmps_highband_lower_limit > nominal_dp(2)) THEN
            dmps_sp_max = minloc(abs(get_dp()-dmps_highband_lower_limit),1)
            print FMT_SUB, 'Upper bins replaced from indices '//TRIM(i2chr(dmps_sp_max))//' to '//TRIM(i2chr(n_bins_par))
        ELSE IF (dmps_highband_lower_limit > 0d0) THEN
            print FMT_FAT0, 'dmps_highband_lower_limit is is too large and will not be used'
            STOP
        END IF

    END IF



    ALLOCATE(conc_vapour(n_cond_tot))
    ALLOCATE(corgwallTeflon(VAPOUR_PROP%n_condorg))
    corgwallTeflon = 0d0

    ! Allocate the change vectors for integration timestep control
    ALLOCATE(d_vap(size(SPC_NAMES)))

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

    ! Save all condensables for easier access
    OPEN(600,file=RUN_OUTPUT_DIR//"/CondensingVapours.txt",status='replace',action='write')
    do ii = 1,size(VAPOUR_PROP%vapour_names)
        WRITE(600,'(a)') TRIM(VAPOUR_PROP%vapour_names(ii))
    end do
    CLOSE(600)


    ! Open text files to save particle size ditribution in easily accessible format
    ! Sumfile that uses same format as SMEAR sumfiles, (with exception of time resolution; here model save interval is used)
    ! Format is:
    ! 000------------------000------------bin_1 diameter---bin_2 diameter---bin_3 diameter . bin_n diameter
    ! Time 0 (s)--Total particle number---dN/log10(dDp)---dN/log10(dDp)-----dN/log10(dDp) . . dN/log10(dDp)
    ! Time 1 (s)--Total particle number---dN/log10(dDp)---dN/log10(dDp)-----dN/log10(dDp) . . dN/log10(dDp)
    OPEN(601,file=RUN_OUTPUT_DIR//"/particle_conc.sum",status='replace',action='write')
    WRITE(601,*) 0d0,0d0,get_dp()

    ! Datfile where PSD is in linear form. Format is:
    ! 000----------bin_1 diameter---bin_2 diameter---bin_3 diameter . bin_n diameter
    ! Time 0 (s)---------dN--------------dN----------------dN . . . . . . . dN
    ! Time 1 (s)---------dN--------------dN----------------dN . . . . . . . dN
    OPEN(604,file=RUN_OUTPUT_DIR//"/particle_conc.dat",status='replace',action='write')
    WRITE(604,*) 0d0,get_dp()

END IF

! Save backup from the initfile
OPEN(UNIT=606, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=ioi)
open(unit=607, file=RUN_OUTPUT_DIR//'/InitBackup.txt')
DO while (ioi == 0)
    read(606,'(a)', iostat = ioi) inibuf
    write(607,'(a)') TRIM(inibuf)
end do
close(606)
close(607)


! Open netCDF files
CALL OPEN_FILES( RUN_OUTPUT_DIR, Description, MODS, CH_GAS, VAPOUR_PROP)

! If wait_for was defined in user options, wait for a sec
CALL PAUSE_FOR_WHILE(wait_for)

if (TRIM(INITIALIZE_WITH) /= '') THEN
    if (current_PSD%psd_style/=2) THEN
        CALL INITIALIZE_WITH_LAST(CURRENT_PSD%composition_fs, CURRENT_PSD%conc_fs,CH_GAS)
    ELSE! if (current_PSD%psd_style==2) THEN
        CALL INITIALIZE_WITH_LAST(CURRENT_PSD%composition_ma, CURRENT_PSD%conc_ma,CH_GAS)
    END IF
END IF

write(*,*) ''

open(unit=608, file=RUN_OUTPUT_DIR//'/optimization.txt',status='replace',action='write')

if (Use_speed) THEN
    DT_0 = GTIME%dt
    print FMT_HDR, 'Simulation time step will be optimized for speed and precision'
    ! print '("| ",a,6(f7.2, " % "),t100,"|")', 'precision limits: ', change_range(1:3,:)*100
    ! print '("| ",a,4(f7.2, " % "),t100,"|")', '                  ', change_range(4:5,:)*100
    print FMT_HDR, 'First running a few n_of_Rounds without speeding....'
    change_range = TRANSPOSE(RESHAPE([diameter_prec_def,pnumber_prec_def,vapour_prec_def], [2,3]))
    DO i=1,size(change_range,1)
        print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for '//range_names(i), change_range(i,:)
    END DO
    print FMT_LEND
    change_range = change_range *1d-2
    Handbrake_on = .true.
    Use_speed = .false.
else
    print FMT_HDR, 'Beginning simulation with constant timestep'
    write(608, *) 'Time step was not optimized, dt: ',GTIME%dt
    flush(608)
end if

if (start_time_s>0) then
    write(*, '(a)', advance='no') 'Starting simulation at: ', ACHAR(10)
    print FMT_TIME, GTIME%hms
end if

call cpu_time(cpu1) ! For efficiency calculation


! =================================================================================================
! =================================================================================================
DO WHILE (GTIME%SIM_TIME_S - GTIME%sec > -1d-12) ! MAIN LOOP STARTS HERE
! =================================================================================================
    ! =================================================================================================
    ! Store the current state of the aerosol
    ! (i.e. everything that is potentially changed by aerosol dynamics)
    ! =================================================================================================
    ! if (gtime%sec == 26290) GTIME%dt = 5d0
    ! if (gtime%sec == 26290) print*, 'Time step is now 5 sec'
    old_PSD = current_PSD
    CH_GAS_old = CH_GAS
    CH_RO2_old = CH_RO2
    ! =================================================================================================
    ! =================================================================================================
    ! PREPARE AND PRINT OUT SOME FUNDAMENTAL STUFF
    ! =================================================================================================
    ! =================================================================================================

    ! Print time header in a nice way, start with empty row
    if (GTIME%printnow) THEN
        WRITE(*,*) ! Print time
        print FMT_TIME, GTIME%hms
        ! To compare real time vs simulation time, timer is stopped in the beginning
        call cpu_time(cpu2)
    ! if --gui flag was used, print a dot and EOL so the STDOUT-reading in Python GUI will be smoother.
    ELSE
      if (ingui) print'(a)', '.'
    END IF

    ! GC_AIR_NOW, TEMPK, PRES and RH are calculated as global variables and are available everywhere.
    ! H2SO4 CS is based on the input, if provided, otherwise it is calculated from the aerosol population
    ! USE GC_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3
    GTEMPK = interp(timevec, CONC_MAT(:,inm_tempK)) .mod. MODS(inm_tempK)
    GPRES = interp(timevec, CONC_MAT(:,inm_PRES)) .mod. MODS(inm_PRES)
    GRH = interp(timevec, CONC_MAT(:,inm_RH)) .mod. MODS(inm_RH)
    GC_AIR_NOW = C_AIR_cc(GTEMPK, GPRES)

    IF (MODS(inm_CS)%ISPROVIDED) THEN
        GCS = interp(timevec, CONC_MAT(:,inm_CS)) .mod. MODS(inm_CS)
    ELSE
        IF (equal(GTIME%sec, 0d0)) THEN
            GCS = 0d-3
        ELSEIF(Aerosol_flag) THEN
            CONTINUE ! Now GCS will be updated in the aerosol module
        ELSEIF (ACDC) THEN
            PRINT FMT_WARN0, 'Condensation sink is not provided as input, and aerosol module is off -> NO CS'
        ELSE
            IF (GTIME%printnow) print FMT_WARN0, 'NO CS provided or calculated'
        END IF
    END IF
    ! Calculate Water vapour pressure and concentration
    CALL WATER(ES,EW,CH_H2O)

    ! Assign values to input variables. N_VARS will cycle through all variables that user can provide
    ! or tamper, and leave zero if no input or mod was provided
    DO I = 1, N_VARS
        IF (MODS(i)%ISPROVIDED) TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I)) .mod. MODS(I)
        IF (I == inm_CS) TSTEP_CONC(I) = GCS
    END DO
    ! =================================================================================================

    ! =================================================================================================
    ! =================================================================================================
    ! Chemistry
    ! =================================================================================================
    ! =================================================================================================
    IF (Chemistry_flag) THEN

        if (GTIME%hrs<FLOAT_CHEMISTRY_AFTER_HRS) THEN
            DO I = 1, N_VARS ! <-- N_VARS will cycle through all input variables
                IF (INDRELAY_CH(I)>0) THEN ! <-- this will pick those that were paired in CHECK_INPUT_AGAINST_KPP
                    IF (I == inm_H2SO4 .and. .not. model_H2SO4) THEN
                        CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                    ELSE IF (I /= inm_H2SO4) THEN
                        CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                    END IF
                END IF
            END DO
        END IF

        ! Determine the SWR direction
        IF (GTIME%JD==0) THEN
            CH_Beta = 90d0
        ELSE
            ! Solar angle above horizon. For this to properly work, lat, lon and Date need to be defined in INIT_FILE
            call BETA(CH_Beta)
        END IF

        Call CHEMCALC(CH_GAS, GTIME%sec, (GTIME%sec + GTIME%dt), GTEMPK, TSTEP_CONC(inm_swr), CH_Beta,  &
                    CH_H2O, GC_AIR_NOW, GCS, TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2)

        if (model_H2SO4) TSTEP_CONC(inm_H2SO4) = CH_GAS(ind_H2SO4)

        IF (EXTRAFUNCS) CALL AFTER_CHEM(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
    END IF ! IF (Chemistry_flag)
    ! =================================================================================================

    ! =================================================================================================
    ! =================================================================================================
    ! NUCLEATION
    ! Formation rates can be calculated using ACDC or the organic nucleation parametrisation.In
    ! addition to these, formation rate (in 1/s/cm³) can be sent in with NUC_RATE_IN.  These
    ! will be added together, so if one or more has to be omitted, turn them off in options
    ! (like all input, NUC_RATE_IN is zero if it is not sent in). When RESOLVE_BASE is used, NUC_RATE_IN
    ! has to be provided as a reference formation rate. If only one of the bases is to be iterated, the
    ! other rate is coming from ACDC. Using the Ratio in Solve Bases (GUI-> ADVANCED) the ration between
    !
    ! =================================================================================================
    IF (.not. PRCION%err) THEN
        if (ACDC) THEN
            CALL ACDC_J(TSTEP_CONC, acdc_iterations)
            J_TOTAL_M3 = J_ACDC_DMA_M3 + J_ACDC_NH3_M3 ! [particles/s/m^3]
            IF (RESOLVE_BASE) THEN
                J_TOTAL_M3 = TSTEP_CONC(inm_JIN) * 1d6 ! [particles/s/m^3]
            ELSE
                J_TOTAL_M3 = J_TOTAL_M3 + TSTEP_CONC(inm_JIN) * 1d6 ! [particles/s/m^3]
            END IF
        else
          ! Only use input formation rate
          J_TOTAL_M3 = TSTEP_CONC(inm_JIN) * 1d6
        END if

        IF (ORG_NUCL) CALL ORGANIC_NUCL(J_TOTAL_M3)

        IF (EXTRAFUNCS) CALL AFTER_NUCL(TSTEP_CONC,CH_GAS,J_TOTAL_M3)

    END if

    if (GTIME%savenow .and. RESOLVE_BASE) CALL Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)
    ! =================================================================================================


    ! =================================================================================================
    ! =================================================================================================
    ! AEROSOL PART STARTS HERE
    ! =================================================================================================
    ! =================================================================================================
    IF (Aerosol_flag.and.(.not. PRCION%err)) THEN

        ! Read in background particles
        IF (N_MODAL>0) THEN
            call Multimodal(MMODES, get_dp(), conc_fit, N_MODAL)
            CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit*dmps_multi)
            do i = 1, n_bins_par
                CALL set_composition(i,nominal_dp(i))
            end do
            if (GTIME%hrs >= DMPS_read_in_time) N_MODAL = -1d0
        END IF


        if (use_dmps .and. GTIME%min >= (dmps_ln*dmps_tres_min)) THEN
            ! Particles are read in from measuremensts throughout the simulation and saved to Particles.nc for comparison
            CALL GeneratePSDfromInput( BG_PAR%sections,  BG_PAR%conc_matrix(min(dmps_ln+1, size(BG_PAR%time, 1)),:), conc_fit )
            conc_fit = conc_fit*dmps_multi

            ! If the mdodel is still in initialization mode, replace the model particles with measured
            if (GTIME%hrs < DMPS_read_in_time) THEN
                ! Print user info
                if (gtime%printnow) print FMT_SUB, 'Replacing PSD with background particles from '//TRIM(i2chr(dmps_ln+1))&
                                                    //'. measurement (row '//TRIM(i2chr(dmps_ln+2))//' in file)'
                CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit)
                do i = 1, n_bins_par
                    CALL set_composition(i,nominal_dp(i))
                end do
            END IF

            ! If the mdodel is still in initialization mode, replace the model particles with measured
            if (use_dmps_partial &
            .and. (GTIME%hrs .ge. DMPS_read_in_time)&
            .and. (GTIME%hrs .lt. END_DMPS_SPECIAL)&
            ) THEN
                ! Print user info
                if (gtime%printnow) print FMT_MSG, 'Replacing PSD partially with background particles from '//TRIM(i2chr(dmps_ln+1))&
                                                    //'. measurement (row '//TRIM(i2chr(dmps_ln+2))//' in file)'
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
                    do i = dmps_sp_max, n_bins_par
                        CALL set_composition(i,nominal_dp(i))
                    end do

                END IF

            END IF

            ! Next time next line from measurement etc. is read
            dmps_ln = dmps_ln + 1
        END IF

        ! ..........................................................................................................
        ! ADD NUCLEATED PARTICLES TO PSD, IN THE 1st BIN. This using the mixing method, where two particle distros
        ! are mixed together. NOTE that this step is always done if Aerosol_flag is on, even if ACDC_flag is off. In
        ! case no NPF is wanted, turn off ACDC and make sure NUC_RATE_IN is 0
        dmass = 0d0

        dmass(1,VAPOUR_PROP%ind_GENERIC) = nominal_dp(1)**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)
        dconc_dep_mix = 0d0
        dconc_dep_mix(1) = J_TOTAL_M3*GTIME%dt
        ! Negative mixing ratio makes this nucleation
        mix_ratio = -1d0
        CALL Mass_Number_Change('mixing')
        ! Update current_psd
        current_PSD = new_PSD


        ! ..........................................................................................................
        ! CONDENSATION
        if (Condensation .and.(.not. PRCION%err) .and. mod(int(GTIME%sec/GTIME%dt),speed_up(PRCION%con)) == 0) THEN
! print*,'in cond at ',GTIME%sec,Gtime%dt
! IF (n_of_Rounds > 0) PRINT*, 'nO OF ROUNDS:', n_of_Rounds
            ! Pick the condensables from chemistry and change units from #/cm^3 to #/m^3
            conc_vapour = 0d0
            conc_vapour(1:VAPOUR_PROP%n_condorg) =  CH_GAS(index_cond)*1D6 ! mol/m3

            ! Poor sulfuric acid always wants special treatment
            conc_vapour(n_cond_tot) = CH_GAS(ind_H2SO4)*1d6

            ! Update vapour pressures for organics
            VAPOUR_PROP%c_sat(1:VAPOUR_PROP%n_condorg) =  saturation_conc_m3( &
                                VAPOUR_PROP%psat_a(1:VAPOUR_PROP%n_condorg),  &
                                VAPOUR_PROP%psat_b(1:VAPOUR_PROP%n_condorg),  &
                                GTEMPK) * VP_MULTI

            CALL UPDATE_MOLECULAR_DIFF_AND_CSPEED(VAPOUR_PROP)

            IF (CHEM_DEPOSITION) CALL CALCULATE_CHEMICAL_WALL_LOSS(conc_vapour(1:VAPOUR_PROP%n_condorg))

            dmass = 0d0

            CALL Condensation_apc(VAPOUR_PROP,conc_vapour,dmass, GTIME%dt*speed_up(PRCION%con),d_dpar,d_vap)

            ! ERROR HANDLING
            IF ((maxval(ABS(d_dpar)) > change_range(1,2) .or. sum(ABS(d_vap))/vapour_prop%n_condtot > change_range(3,2)) .and. use_speed) THEN   !if the changes in diameter are too big
                ! IF (n_of_Rounds>=0) THEN
                print*, 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                PRCION%err = .true.
                PRCION%proc = PRCION%con
                ! end if
                print*, 'dpar ja höy', maxval(ABS(d_dpar)) > change_range(1,2), sum(ABS(d_vap))/vapour_prop%n_condtot > change_range(3,2)
                print*, 'dpar ja höy', maxval(ABS(d_dpar)),  change_range(1,2), sum(ABS(d_vap))/vapour_prop%n_condtot,  change_range(3,2)
                print*, maxloc(ABS(d_vap)), VAPOUR_PROP%vapour_names( maxloc(ABS(d_vap)))

                IF (maxval(ABS(d_dpar)) > change_range(1,2)) THEN
                    PRCION%err_text = 'Too large diameter change: '//f2chr(maxval(abs(d_dpar)))   !d_dpar(maxloc(abs(d_dpar)))
                ELSE IF (sum(ABS(d_vap))/vapour_prop%n_condtot > change_range(3,2)) THEN
                    PRCION%err_text = 'Too large vapour concentration change: '//f2chr(sum(ABS(d_vap))/vapour_prop%n_condtot)//', UL: '//f2chr(change_range(3,2))   !d_vap(maxloc(abs(d_vap)))
                ELSE
                    PRCION%err_text = 'Too large diameter and vapour concentration change: '//f2chr(maxval(abs(d_dpar)))//', '//f2chr(sum(ABS(d_vap))/vapour_prop%n_condtot)   !d_vap(maxloc(abs(d_vap)))
                END IF
                PRINT*,'Precision error in condensation at '//GTIME%hms//', '//PRCION%err_text
                ! PRINT*,'what, where:', d_dpar(maxloc(abs(d_dpar))), maxloc(abs(d_dpar))
                dmass = 0.d0

            ELSE  ! everything fine -> apply changes
                ! Calculate growth rates
                IF (GTIME%printnow .and. CALC_GR) THEN
                    CALL PRINT_GROWTH_RATE
                END IF
! if (sum(ABS(d_vap))/vapour_prop%n_condtot>0.1) print*, gtime%hms,GTIME%sec, 'max höyrynmuutos', sum(ABS(d_vap))/vapour_prop%n_condtot
! print*, gtime%hms,GTIME%sec, 'max höyrynmuutos', sum(ABS(d_vap))/vapour_prop%n_condtot
                ! Distribute mass
                CALL Mass_Number_Change('condensation')

                ! Update current_psd
                current_PSD = new_PSD

                ! Update vapour concentrations to chemistry
                CH_GAS(index_cond) = conc_vapour(1:VAPOUR_PROP%n_condorg) *1D-6
                ! Again sulfuric acid always needs special treatment
                CH_GAS(ind_H2SO4) = conc_vapour(n_cond_tot)*1d-6
                ! Check whether timestep can be increased:
                IF (maxval(ABS(d_dpar)) < change_range(1,1) .and. sum(ABS(d_vap))/vapour_prop%n_condtot < change_range(3,1) .and. use_speed) THEN
                    speed_up(PRCION%con) = speed_up(PRCION%con) * 2
                    Print*,GTIME%hms//'-> Cond. speed multiplier now:', speed_up(PRCION%con)
                    WRITE(608,*) GTIME%hms//'-> Cond. speed multiplier now:', speed_up(PRCION%con)
                END IF
            END IF

        end if ! end of condensation


        ! ..........................................................................................................
        ! COAGULATION
        if (Coagulation .and.(.not. PRCION%err) .and. mod(int(GTIME%sec/GTIME%dt),int(speed_up(PRCION%coa))) == 0) then
            ! Solve particle coagulation
            dconc_coag = 0.d0
            Call Coagulation_routine(dconc_coag, GTIME%dt*speed_up(PRCION%coa),d_npar)

            ! ERROR HANDLING
            IF (maxval(ABS(d_npar)) > change_range(2,2) .and. use_speed) THEN   ! if the changes in particle numbers are too big
                PRCION%err = .true.
                PRCION%proc = PRCION%coa
                !PRINT*,'ERROR',maxloc(abs(d_dpar)), d_dpar(maxloc(abs(d_dpar))),maxloc(abs(d_vap)), d_vap(maxloc(abs(d_vap)))
                PRCION%err_text = 'Too large particle number change'   !d_dpar(maxloc(abs(d_dpar))
                PRINT*,'Precision error in coagulation at '//GTIME%hms//', '//PRCION%err_text
                dconc_coag = 0.d0
            ELSE  ! everything fine -> apply changes
                ! Distribute mass
                Call Mass_Number_Change('coagulation')
                ! Update PSD with new concentrations
                current_PSD = new_PSD
                !Check whether timestep can be increased:
                IF (maxval(ABS(d_npar)) < change_range(2,1) .and. use_speed) THEN
                  speed_up(PRCION%coa) = speed_up(PRCION%coa) * 2
                  Print*,GTIME%hms//'-> Coag. speed multiplier now:', speed_up(PRCION%coa)
                  WRITE(608,*) GTIME%hms//'-> Coag. speed multiplier now:', speed_up(PRCION%coa)
                END IF
            END IF
        end if ! end of coagulation


        ! ..........................................................................................................
        ! Deposition
        if ((Deposition .or. LOSSES_FILE /= '') .and.(.not. PRCION%err) .and. mod(int(GTIME%sec/GTIME%dt),int(speed_up(PRCION%dep))) == 0) then
            ! Solve particle coagulation
            dconc_dep_mix = 0.d0

            IF (LOSSES_FILE /= '') THEN

                ! The losses file is interpolated spatially and temporally to fit the current bin structure and time
                losses_fit = [((INTERP( PAR_LOSSES%sections,                                                            &
                             [((INTERP(PAR_LOSSES%time, PAR_LOSSES%conc_matrix(:,i) )), i=1,size(PAR_LOSSES%time, 1))], &
                             timein=nominal_dp(j))                                                                      &
                             ), j=1,n_bins_par)]

                ! Deposited concentratios calculated here
                dconc_dep_mix = get_conc() * (1 - EXP(-losses_fit*GTIME%dt))


            ELSE
                call deposition_velocity(nominal_dp,ustar,CHAMBER_CIRCUMFENCE*CHAMBER_HEIGHT,Chamber_floor_area,Chamber_floor_area,&
                                        Chamber_floor_area*CHAMBER_HEIGHT,GTEMPK,GPRES,E_field,dconc_dep_mix,GTIME%dt*speed_up(PRCION%dep))
            END IF
            ! ERROR HANDLING; Check whether changes are within limits:
            d_npar = 0.d0
            WHERE (get_conc()>0d0) d_npar = dconc_dep_mix/get_conc()

            IF (maxval(ABS(d_npar)) > change_range(2,2) .and. use_speed) THEN   ! if the changes in diameter are too big
                PRCION%err = .true.
                PRCION%proc = PRCION%dep

                PRCION%err_text = 'Too large particle number change' !d_dpar(maxloc(abs(d_dpar))
                PRINT*,'Precision error in deposition at '//GTIME%hms//', '//PRCION%err_text
                dconc_dep_mix = 0.d0
            ELSE  ! everything fine -> apply changes
                ! Distribute mass
                Call Mass_Number_Change('deposition')
                ! Update PSD with new concentrations
                current_PSD = new_PSD
                !Check whether timestep can be increased:
                IF (maxval(ABS(d_npar)) < change_range(2,1) .and. use_speed) THEN
                  speed_up(PRCION%dep) = speed_up(PRCION%dep) * 2
                  Print*,GTIME%hms//'-> Depos. speed multiplier now:', speed_up(PRCION%dep)
                  WRITE(608,*) GTIME%hms//'-> Depos. speed multiplier now:', speed_up(PRCION%dep)
                END IF
            END IF


            ! Update PSD with new concentrations
            current_PSD = new_PSD
        end if ! end of Deposition

    end if
    ! End Aerosol =====================================================================================

    if (gtime%printnow) print*, gtime%hms,GTIME%sec, 'max höyrynmuutos', sum(ABS(d_vap))/vapour_prop%n_condtot
    ! SPEED handling
    IF (PRCION%err) THEN  ! In case of a timestep error (i.e. too large changes in aerosol dynamics)
        print*, n_of_Rounds

        PRINT*,'Precision error is active', PRCION%err, TRIM(PRCION%pr_name(PRCION%proc))
        CALL error_handling(PRCION, speed_up)
        ! RESET ERROR
        PRCION%err = .false.
        PRCION%proc = 0
        PRCION%err_text = ''
        ! Reset the system to the state at previous timestep
        current_PSD = old_PSD
        CH_GAS = CH_GAS_old
        CH_RO2 = CH_RO2_old
        ! Reset relative change vectors documenting the precision:
        d_dpar = 0.d0
        d_npar = 0.d0
        d_vap  = 0.d0

        PRINT*,'Error handled:', PRCION%err, PRCION%proc
        write(608,'("error state: ",l,", error process: ",i0)') PRCION%err, PRCION%proc
        flush(608)
        print*,''

    ! If there was no error during the actual timestep
    ELSE
        ! PRCION%err = .false.
        ! PRCION%proc = 0
        ! PRCION%err_text = ''
        ! IF (Use_speed .and. n_of_Rounds<0) n_of_Rounds = n_of_Rounds+1
        ! IF (Use_speed .and. n_of_Rounds>0) then
        !     n_of_Rounds = -20
        !     Use_speed = .FALSE.
        ! end if
        ! Write printouts to screen
        if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
        ! if (GTIME%printnow) print*, GTIME%dt, speed_up

        ! Save to netcdf
        if (GTIME%savenow) THEN

            if (Aerosol_flag) THEN
                ! PSD is also saved to txt
                WRITE(601,*) GTIME%sec, sum(get_conc()*1d-6), get_conc()*1d-6 / LOG10(bin_ratio)
                WRITE(604,*) GTIME%sec, get_conc()*1d-6
                save_measured = conc_fit/dmps_multi

            END IF

            CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3_M3, J_ACDC_DMA_M3, VAPOUR_PROP, save_measured,&
                            1d9*3600/(GTIME%dt*speed_up(PRCION%con))*get_dp()*d_dpar)

        END IF

        if (GTIME%dt<DT_0 .and. minval(speed_up(optis_in_use))>1 .and. MODULO(int(GTIME%sec), int(DT_0)*2)==0) THEN
            print*, GTIME%hms//': Upping the main dt from ', GTIME%dt, 'to', GTIME%dt*2
            GTIME%dt = GTIME%dt * 2
            speed_up = speed_up / 2
        END IF

        ! Add main timestep to GTIME'
        GTIME = ADD(GTIME)


    END IF ! SPEED handling

    if (Handbrake_on) THEN
        n_of_Rounds = n_of_Rounds +1
        IF (n_of_Rounds > 30) THEN
            Handbrake_on = .false.
            Use_speed = .true.
            n_of_Rounds = 0
            print FMT_MSG, 'Engaging speed...'
        END IF
    ! else
    !     IF (n_of_Rounds > 2) THEN
    !         print FMT_MSG, 'Putting on the handbrake...'
    !         Handbrake_on = .true.
    !         Use_speed = .false.
    !         n_of_Rounds = 0
    !     END IF
    END IF

END DO	! Main loop ends
! ==================================================================================================================


CALL PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

if (Aerosol_flag) THEN
    CLOSE(601)
    CLOSE(604)
END IF


! Close output file netcdf
CALL CLOSE_FILES(RUN_OUTPUT_DIR)
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
SUBROUTINE error_handling(PRCION,speed_up)

    IMPLICIT none
    INTEGER ::  i, speed_up(:)
    type(error_type) :: PRCION

    ! Write some error information to file
    write(*,FMT_WARN0) 'Precision error at time '//GTIME%hms//' ('//i2chr(int(GTIME%sec))//' sec)'
    write(608,*) 'Precision error at time '//GTIME%hms//' ('//i2chr(int(GTIME%sec))//' sec)'
    write(*,FMT_SUB) 'current simulation timestep [s]'//f2chr(GTIME%dt)
    write(608,*) 'current simulation timestep [s]', GTIME%dt
    write(*,FMT_SUB) 'Process causing troubles: '//TRIM(PRCION%pr_name(PRCION%proc))
    write(608,*) 'Cause of problem: ',  TRIM(PRCION%pr_name(PRCION%proc))
    write(*,FMT_SUB) 'Cause of problem: '//trim(PRCION%err_text)
    write(608,*) '  ERROR TYPE: ', trim(PRCION%err_text)
    ! reduce the speed_up factor or, if necessary, the integration time step:
    IF (speed_up(PRCION%proc) > 1) THEN
        speed_up(PRCION%proc) = speed_up(PRCION%proc) / 2
        write(608,*) '  => reduce speed_up:',speed_up(PRCION%proc)*2, '->',speed_up(PRCION%proc)
        write(*,FMT_SUB) '  => reduce speed_up:'//i2chr(speed_up(PRCION%proc)*2)//'->'//i2chr(speed_up(PRCION%proc))
    ELSE
        if (GTIME%dt/5>=DT_0) THEN
            GTIME%dt = GTIME%dt / 2.d0
            write(608,*) '  => reduce dt [s]:', GTIME%dt
            write(*,FMT_WARN0) '  => reduce dt [s]:'//f2chr(GTIME%dt)
        ELSE
            write(*,FMT_WARN0) '  Can not reduce dt any more'//f2chr(GTIME%dt)
        END IF
        print*, gtime%hms,GTIME%sec, 'max höyrynmuutos', sum(ABS(d_vap))/vapour_prop%n_condtot
        ! DO i = 1,size(speed_up)
        !     IF (i /= PRCION%proc) speed_up(i) = speed_up(i) * 2
        ! END DO
        n_of_Rounds = n_of_Rounds + 1
    END IF
    write(608,*) ''
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
! J_ACDC_NH3_M3:        Particle formation rate due to ammonia [1/s/m3]. Sum of J_NH3_BY_IONS
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
! J_ACDC_DMA_M3:        Particle formation rate due to DMA [1/s/m3]
! acdc_cluster_diam: Outgrowing cluster diameter [m]. The cluster typically has 5 to 6 H2SO4 in it
!
! =================================================================================================

SUBROUTINE ACDC_J(C,iters)
    implicit none
    REAL(dp) :: c_org = 0d0
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
        CALL get_acdc_J(H2SO4,NH3,c_org,GCS,C(inm_TEMPK),IPR,GTIME,&
            ss_handle,J_ACDC_NH3_M3,acdc_cluster_diam, J_NH3_BY_IONS, iters)
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
    END IF
    ! Speed up program by ignoring nucleation when there is none
    if ((DMA > 1d6 .and. H2SO4>1d9) .or. (.not. skip_acdc)) THEN
        ! Idea of iteration is to not necessarily reach steady state but still get a "better" estimate of instantenous formation rate
        CALL get_acdc_D(H2SO4,DMA,c_org,GCS,C(inm_TEMPK),GTIME,ss_handle,J_ACDC_DMA_M3,acdc_cluster_diam, iters)
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'DMA IGNORED'
    END IF

    ! The first time ACDC is run, it is run to steady state, after that the user supplied option is used
    if (first_time) THEN
        ss_handle = ACDC_solve_ss
        if (RESOLVE_BASE) ss_handle = .TRUE.
        first_time = .false.
    end if

END SUBROUTINE ACDC_J  ! END ACDC Nucleation


! =================================================================================================
!
! Organic nucleation based on Roldin et al., Nat Comm. 2019
! .................................................................................................
SUBROUTINE ORGANIC_NUCL(J_TOTAL_M3)
    real(dp) :: J_TOTAL_M3
    real(dp) :: ORGS,kJ,dH,J15
    real(dp), parameter :: dG = -15.1, dS = -61.1
    integer :: i,jj,n
    integer, save, allocatable  :: inds(:)
    logical, save  :: first_run = .True.
    character(25)  :: name

    if (first_run) THEN
        OPEN(UNIT=609, FILE='ModelLib/nucl_homs.txt', STATUS='OLD', ACTION='READ', iostat=jj)
        n = rowcount(609)
        print FMT_MSG, 'Using parametrisation for organic nucleation with '//i2chr(n)//' nucleating compounds'
        allocate(inds(n))
        DO i=1,n
            read(609, *) name
            if (IndexFromName(TRIM(name), SPC_NAMES)>0) inds(i) = IndexFromName(TRIM(name), SPC_NAMES)
        END DO
        first_run = .False.
    END IF

    ORGS = sum(CH_GAS(inds))

    dH = dG - dS*GTEMPK
    kJ = 5d-13*EXP(-dH/Rg * (1/GTEMPK - 1/298d0))

    J15 = kJ * ORGS*TSTEP_CONC(inm_H2SO4)
    if (GTIME%printnow) print FMT_MSG, 'Organic formation rate: '//TRIM(ADJUSTL(f2chr(J15)))//' [/s/cm3]'
    J_TOTAL_M3 = J_TOTAL_M3 + J15*1d6

END SUBROUTINE ORGANIC_NUCL

SUBROUTINE CALCULATE_CHEMICAL_WALL_LOSS(conc)
    IMPLICIT NONE
    real(dp) :: conc(:)
    real(dp) :: conc_dummy(VAPOUR_PROP%n_condorg),corg_old(VAPOUR_PROP%n_condorg)
    ! Reversible wall losses of condensable VOCs:
    real(dp) :: kwallTeflon(VAPOUR_PROP%n_condorg),kwallbTeflon(VAPOUR_PROP%n_condorg)
    real(dp) :: ke_chamber
    integer  :: nc ! short for number of condensible organic vapours

    nc = VAPOUR_PROP%n_condorg

    ! kwallTeflon based on McMurry PH & Grosjean D (1985) Gas and aerosol wall losses in Teflon film smog chambers. Environ. Sci. Technol. 19(12):1176-1182.
    ke_chamber=EDDYK    ! Coefficient of eddy diffusion (s^-1), usually treated as an unknown fitting parameter

    kwallTeflon=((Av+Au+Ad)/V_chamber)*(ALPHAWALL*VAPOUR_PROP%c_speed(1:nc)/4D0)/(1D0+(pi/2D0) &
                * (ALPHAWALL*VAPOUR_PROP%c_speed(1:nc)/(4D0*sqrt(ke_chamber*VAPOUR_PROP%diff(1:nc)))))

    kwallbTeflon=kwallTeflon/(Na/VAPOUR_PROP%c_sat(1:nc)*50D-6) ! Cw=40 mumol/m^3, => ~10 mg/m^3 as in Zang et al., 2014
    conc_dummy = conc

    ! Use a small internal time step of 0.01 s to update the gas and particle wall pool of condensable organic molec:
    DO j=1,INT(GTIME%dt/0.01D0)
        corg_old=conc_dummy
        conc_dummy=conc_dummy*EXP(-kwallTeflon*0.01D0)+kwallbTeflon*corgwallTeflon*0.01D0
        corgwallTeflon=corgwallTeflon+(corg_old-conc_dummy)
    END DO

    if (GTIME%printnow .and. SUM(conc)>0d0) print*, 'mass lost to walls: ', &
                100-SUM(VAPOUR_PROP%molec_mass(1:nc)*conc_dummy)&
                /SUM(VAPOUR_PROP%molec_mass(1:nc)*conc)*1d2
    conc = conc_dummy  ! Vector with concentrations of condensable organic compounds


END SUBROUTINE CALCULATE_CHEMICAL_WALL_LOSS


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
    IF (inm_NH3   /= 0) print FMT10_3CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3_M3*1d-6, ' [1/s/cm3]','sum(Nn)/N1',clusterbase,' []'
    IF (inm_DMA   /= 0) print FMT10_3CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA_M3*1d-6, ' [1/s/cm3]', 'J-total', J_TOTAL_M3*1e-6,' [1/s/cm3]'
    print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', GCS , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
    if ((GTIME%sec)>0 .and. (cpu2 - cpu1 > 0d0)) print '("| ",a,i0,a,i0.2,a,i0,a,i0.2,t65,a,f7.2,t100,"|")', 'Elapsed time (m:s) ',int(cpu2 - cpu1)/60,':',modulo(int(cpu2 - cpu1),60) ,' Est. time to finish (m:s) ',&
                            int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec))/60,':', MODULO(int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec)),60),&
                            'Realtime/Modeltime: ', (GTIME%sec-start_time_s)/(cpu2 - cpu1)
    ! print FMT_, 'current integration time step', GTIME%dt

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
    if (.not. GTIME%savenow)  THEN
        if (Aerosol_flag) THEN
            ! PSD is also saved to txt
            WRITE(601,*) GTIME%sec, sum(get_conc()*1d-6), get_conc()*1d-6 / LOG10(bin_ratio)
            WRITE(604,*) GTIME%sec, get_conc()*1d-6
            save_measured = conc_fit/dmps_multi
        END IF

        CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3_M3, J_ACDC_DMA_M3, VAPOUR_PROP, save_measured,&
                        1d9*3600/(GTIME%dt*speed_up(PRCION%con))*get_dp()*d_dpar)


    END IF
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
        IF (MODS(I)%ISPROVIDED) THEN
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
                IF (ABS(MODS(I)%SHIFT) > 0d0) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
                print FMT_MSG, 'Good bye.'
                STOP
            END IF
            IF (check == 0 .and. I<LENV .and. &
                ((I==inm_H2SO4) .or. (I==inm_NH3) .or. (I==inm_DMA) .or. (I==inm_SO2) &
                .or. (I==inm_NO) .or. (I==inm_NO2) .or. (I==inm_CO) .or. (I==inm_H2) .or. (I==inm_O3))  ) THEN
                print FMT_WARN0, 'A compound which does not exist in chemistry but could be used elsewhere: '//TRIM(MODS(i)%NAME)//' '
                IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                IF (ABS(MODS(I)%SHIFT) > 0d0) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
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
            WRITE(*,*)
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
    write(608,*) '---------- SIMULATION REACHED END SUCCESFULLY ------------'
    close(608)

END SUBROUTINE FINISH


SUBROUTINE PRINT_GROWTH_RATE
    IMPLICIT NONE
    CHARACTER(len=60) :: fmt
    CHARACTER(len=90) :: fmt2
    INTEGER :: nb
    nb = size(GGR)

    ii = MINLOC(ABS(nominal_dp-GR_bins(1)),1)
    GGR(1) = 1d9 * 3600/(GTIME%dt*speed_up(PRCION%con)) * (sum((nominal_dp(1:ii) * d_dpar(1:ii)))/ii)
    if (nb>1) THEN
        DO i=2,size(GGR)
            jj = MINLOC(ABS(nominal_dp-GR_bins(i)),1)
            GGR(i) = 1d9 * 3600/(GTIME%dt*speed_up(PRCION%con)) * (sum((nominal_dp(ii:jj) * d_dpar(ii:jj)))/(jj-ii))
            ii = jj
        END DO
    END IF
    print FMT_MSG, 'Growth rates (nm/h) between sizes (nm, starting from smallest bin)'
    fmt2 = '("| ",a,f6.1,t18,'//TRIM(i2chr(MIN(SIZE(GGR,1), 3)))//'("   -->",f6.1,),a,t100,"|")'
    fmt  = '("| ",a,t18,'//TRIM(i2chr(MIN(SIZE(GGR,1), 3)))//'(es9.2,3(" ")),a,t100,"|")'

    print fmt2, 'Sizes:  ',min_particle_diam*1d9, GR_bins(1:MIN(SIZE(GGR,1), 3))*1d9, ' [nm]'
    print fmt, 'GR:', GGR(1:MIN(SIZE(GGR,1), 3)), ' [nm/h]'

END SUBROUTINE PRINT_GROWTH_RATE


END PROGRAM ARCA_main
