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

INTEGER               :: dmps_ln = 0, dmps_ln_old      ! line number from where background particles are read from
INTEGER               :: dmps_sp_min = 0, dmps_sp_max = 0 ! Indices for dmps_special

! Variables related to chemistry module
REAL(dp) :: CH_RO2,CH_RO2_old      ! RO2 concentration in [molecules / cm^3]
REAL(dp) :: CH_H2O      ! H20 concentration in [molecules / cm^3]
REAL(dp) :: CH_Beta     ! solar zenit angle
REAL(dp) :: EW          ! Water content in Pa
REAL(dp) :: ES          ! Saturation vapour pressure
LOGICAL  :: acdc_goback = .false.
! Variables related to aerosol module
REAL(dp), dimension(:), allocatable:: conc_vapour

! Transient variables
CHARACTER(:), allocatable:: RUN_OUTPUT_DIR ! Saves the output directory relative to this executable
CHARACTER(1000) :: inibuf       ! Buffer to save backp from the INITfile that was called
INTEGER         :: I,II,J,JJ    ! Loop indices
INTEGER         :: ioi          ! iostat variable
INTEGER(dint)   :: add_rounds=1          ! iostat variable
REAL(dp)        :: cpu1, cpu2,cpu3,cpu4 ! CPU time in seconds
REAL(dp)        :: errortime(3) = -9999d0 ! CPU time in seconds
! LOGICAL         :: in_turn(4) = .true.

! Temporary variables -> will be replaced
REAL(dp), ALLOCATABLE :: nominal_dp(:)    ! array with nominal diameters. Stays constant independent of PSD_style

LOGICAL                 :: Handbrake_on = .false.
INTEGER(dint)           :: n_of_Rounds = 0_dint
REAL(dp), ALLOCATABLE   :: J_distr(:)
REAL(dp), ALLOCATABLE   :: d_dpar(:)  ! array reporting relative changes to the diameter array within a single timestep
REAL(dp), ALLOCATABLE   :: d_npar(:)  ! array reporting relative changes to the particle number array within a single timestep
REAL(dp), ALLOCATABLE   :: d_vap(:)   ! array reporting relative changes to the vapour concentration array within a single timestep
REAL(dp) :: Av
REAL(dp) :: Au
REAL(dp) :: Ad
REAL(dp) :: V_chamber
REAL(dp) :: E_field = 0d0
REAL(dp), ALLOCATABLE :: corgwallTeflon(:)
! CHARACTER(len=8) :: CHEM
CHARACTER(len=64) :: CurrentChemistry


#ifdef CHEM
CurrentChemistry = CHEM
#endif

! Welcoming message
print'(a,t35,a)', achar(10),  '--~:| ARCA BOX MODEL 0.9 |:~--'//achar(10)
print FMT_HDR, 'Compiled with "'//TRIM(currentChemistry)//'" chemistry module'
print*, ''
! ==================================================================================================================
! Declare most variables and read user input and options in input.f90
CALL READ_INPUT_DATA
! ==================================================================================================================

Av = CHAMBER_CIRCUMFENCE*CHAMBER_HEIGHT
Au = Chamber_floor_area
Ad = Chamber_floor_area
V_chamber = Chamber_floor_area*CHAMBER_HEIGHT

CALL LINK_VARIABLES

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

    ALLOCATE(J_distr(   MINLOC(abs(nominal_dp-nominal_dp(1)*1.15),1)   ))
    do ii=1,size(J_distr)
        J_distr(ii) = 0.5*exp(-0.6d0*ii)
    end do
    J_distr(1) = J_distr(1) + 1-sum(J_distr)
    print FMT_MSG, 'Distributing nucleated particles over '//TRIM(i2chr(size(J_distr, 1)))//' bins.'

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
    d_dpar = 0d0
    d_npar = 0d0

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
    ALLOCATE(d_vap(VAPOUR_PROP%n_condtot))
    d_vap = 0d0

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
CALL OPEN_FILES( RUN_OUTPUT_DIR, Description,CurrentChemistry, MODS, CH_GAS, VAPOUR_PROP)

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
open(unit=610, file=RUN_OUTPUT_DIR//'/optimChanges.txt',status='replace',action='write')
write(610,*) '# time_sec           max_d_vap           max_d_npar           max_d_dpar'

if (Use_speed) THEN

    print FMT_HDR, 'Simulation time step will be optimized for precision'
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for particle diameter',      DDIAM_RANGE
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for particle concentration', DPNUM_RANGE
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for vapour concentration',   DVAPO_RANGE

    ! Input was in percetages, here we change them to fractional
    DDIAM_RANGE = DDIAM_RANGE * 1d-2
    DPNUM_RANGE = DPNUM_RANGE * 1d-2
    DVAPO_RANGE = DVAPO_RANGE * 1d-2
    ! This is used to control the start (or halt) of optimization
    ! Handbrake_on = .true.
    Use_speed = .true.
else
    print FMT_HDR, 'Beginning simulation with constant timestep'
    write(608, *) 'Time step not optimized, dt: ',GTIME%dt
    flush(608)
end if

if (ENABLE_END_FROM_OUTSIDE) THEN
    print FMT_SUB, "Simulation can be stopped using file called ENDNOW.INIT with text 'STOP' in it, saved"
    print FMT_SUB, "in output directory. The program stops in orderly fashion upon next PRINTNOW time and"
    print FMT_SUB, "closes the output files."
end if

print FMT_LEND

if (start_time_s>0) then
    write(*, '(a)', advance='no') 'Starting simulation at: ', ACHAR(10)
    print FMT_TIME, GTIME%hms
end if

call cpu_time(cpu1) ! For efficiency calculation

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@        @@@@@@@@@@@@@@                ,@@@@@@@@@@@@@@@@@@,            @@@@@@@@@@@@@@@@.       ,@@@@@@@@@@@@
!@@@@@@@@@@@          @@@@@@@@@@@@@                     @@@@@@@@@@.                    @@@@@@@@@@@/         *@@@@@@@@@@@
!@@@@@@@@@@            @@@@@@@@@@@@       #@@@@%         &@@@@@@           (%%(       @@@@@@@@@@@#           #@@@@@@@@@@
!@@@@@@@@@      @       @@@@@@@@@@@       #@@@@@@@&       @@@@&        @@@@@@@@@@@@@@@@@@@@@@@@@#      @      #@@@@@@@@@
!@@@@@@@@      &@@       @@@@@@@@@@       #@@@@@@@@       @@@@       &@@@@@@@@@@@@@@@@@@@@@@@@@%      @@@      %@@@@@@@@
!@@@@@@@      .@@@@       @@@@@@@@@       #@@@@@@@       @@@@        @@@@@@@@@@@@@@@@@@@@@@@@@@      @@@@@      @@@@@@@@
!@@@@@@       @@@@@@      ,@@@@@@@@                    ,@@@@@       .@@@@   MAIN LOOP   @@@@@@      #@@@@@#      @@@@@@@
!@@@@@       @@@@@@@@      *@@@@@@@                 @@@@@@@@@        @@@@@ STARTS HERE @@@@@@       @@@@@@@.      @@@@@@
!@@@@                       /@@@@@@       #@@@        @@@@@@@@       #@@@@@@@@@@@@@@@@@@@@@@                       @@@@@
!@@@                         #@@@@@       #@@@@%        @@@@@@#        @@@@@@@@@@@@@@@@@@@@                         @@@@
!@@       &@@@@@@@@@@@@       &@@@@       #@@@@@@        @@@@@@@          .%@@&/       @@@       @@@@@@@@@@@@@       @@@
!@       ,@@@@@@@@@@@@@@       &@@@       #@@@@@@@*        @@@@@@@                     @@       @@@@@@@@@@@@@@@       @@
!       *@@@@@@@@@@@@@@@@       @@@       #@@@@@@@@@        @@@@@@@@@@             @@@@@       @@@@@@@@@@@@@@@@@       @
!=======================================================================================================================

MAINLOOP: DO ! The main loop, runs until time is out. For particular reasons the time is checked at the end of the loop
    ! Here we check which parts of the loop are entered, based on the timestep optimization
    if (Use_speed) THEN
        PRC%in_turn(PRC%cch) = MODULO(n_of_Rounds,speed_up(PRC%cch))==0_dint
        PRC%in_turn(PRC%coa) = MODULO(n_of_Rounds,speed_up(PRC%coa))==0_dint
        PRC%in_turn(PRC%dep) = MODULO(n_of_Rounds,speed_up(PRC%dep))==0_dint
        PRC%in_turn(4) = (PRC%in_turn(PRC%cch).or.(PRC%in_turn(PRC%coa).or.PRC%in_turn(PRC%dep)))

        if (PRC%increase(1).and.PRC%in_turn(1)) THEN
            if (MODULO(n_of_Rounds,speed_up(1)*2)==0_dint) CALL increase_speed(1)
        end if
        if (PRC%increase(2).and.PRC%in_turn(2)) THEN
            if (MODULO(n_of_Rounds,speed_up(2)*2)==0_dint) CALL increase_speed(2)
        end if
        if (PRC%increase(3).and.PRC%in_turn(3)) THEN
            if (MODULO(n_of_Rounds,speed_up(3)*2)==0_dint) CALL increase_speed(3)
        end if

    END IF

    if (PRC%in_turn(4)) THEN
        ! =================================================================================================
        ! =================================================================================================
        ! Store the current state of the aerosol
        ! (i.e. everything that is potentially changed by aerosol dynamics)
        ! =================================================================================================
        old_PSD = current_PSD
        CH_GAS_old = CH_GAS
        CH_RO2_old = CH_RO2
        dmps_ln_old = dmps_ln
        ! =================================================================================================
    ELSE
        print*, 'Entered an unnecessary loop, should not have happened.'
    END IF


    ! Print time header in a nice way, start with empty row
    if (GTIME%printnow) THEN
        WRITE(*,*) ! Print time
        print FMT_TIME, GTIME%hms//' ('//TRIM(i2chr(int(GTIME%sec*1000d0)))//' msec)'
        ! To compare real time vs simulation time, timer is stopped in the beginning
        call cpu_time(cpu2)
        cpu4=cpu2
    ! if --gui flag was used, print a dot and EOL so the STDOUT-reading in Python GUI will be smoother.
    ELSE if (PRC%in_turn(4)) THEN
        if (ingui) print'(a)', '.'
        call cpu_time(cpu3)
        if (cpu3-cpu4>15d0) THEN
            print*, '..Alive at ',GTIME%hms
            cpu4=cpu3
        END IF
    END IF

    ! Assign values to input variables. N_VARS will cycle through all variables that user can provide
    ! or alter, and value leave zero if no input was provided.
    ! GC_AIR_NOW, TEMPK, PRES and RH are saved as global variables and are available everywhere.
    ! H2SO4 CS is based on the input, if provided, otherwise it is calculated from the aerosol population
    ! USE GC_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3

if (PRC%in_turn(4)) THEN
    DO I = 1, N_VARS
        IF (MODS(i)%ISPROVIDED) TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I)) .mod. MODS(I)
        IF (I == 2) THEN
            GTEMPK = TSTEP_CONC(inm_TempK)
            GPRES  = TSTEP_CONC(inm_pres)
            GC_AIR_NOW = C_AIR_cc(GTEMPK, GPRES)
        END IF
    END DO
    DO I = 1, N_VARS
        IF ((MODS(i)%ISPROVIDED).and.(MODS(i)%TIED /= '')) THEN
            TSTEP_CONC(I) = TSTEP_CONC(INDRELAY_TIED(I)) * MODS(I)%MULTI + UCONV(MODS(I)%SHIFT, MODS(I))
        END IF
    END DO

    GRH = TSTEP_CONC(inm_RH)

    IF (MODS(inm_CS)%ISPROVIDED) THEN
        GCS = TSTEP_CONC(inm_CS)
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

END IF

in_turn_cch_1: if (PRC%in_turn(PRC%cch)) THEN

    ! =================================================================================================

    ! =================================================================================================
    ! =================================================================================================
    ! Chemistry
    ! =================================================================================================
    ! =================================================================================================
    CHEMISTRY_ROUTINES: IF (Chemistry_flag) THEN
        ! Calculate Water vapour pressure and concentration
        CALL WATER(ES,EW,CH_H2O)

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

        Call CHEMCALC(CH_GAS, GTIME%sec, (GTIME%sec + GTIME%dt*speed_up(PRC%cch)), GTEMPK, max(0d0,TSTEP_CONC(inm_swr)),&
                    CH_Beta,CH_H2O, GC_AIR_NOW, GCS, TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2)

        ! if ( ( minval(CH_GAS)<-1d0 ).and.GTIME%sec>100d0) print*,'Negative values from chemistry, setting to zero: ', MINVAL(CH_GAS)
        WHERE (CH_GAS<0d0) CH_GAS = 0d0

        if (model_H2SO4) TSTEP_CONC(inm_H2SO4) = CH_GAS(ind_H2SO4)

        IF (AFTER_CHEM_ON) CALL AFTER_CHEM(TSTEP_CONC,CH_GAS,CH_GAS_old,CH_RO2,CH_RO2_old,J_TOTAL_M3)
    END IF CHEMISTRY_ROUTINES
    ! IF (Chemistry_flag)
    ! =================================================================================================
END if in_turn_cch_1

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
in_turn_acdc: if (PRC%in_turn(PRC%coa)) THEN
    ! NUCLEATION_ROUTINES: IF (.not. PRC%err) THEN
        if (ACDC) THEN
            CALL ACDC_J(TSTEP_CONC)
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

        IF (AFTER_NUCL_ON) CALL AFTER_NUCL(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
        if (GTIME%savenow .and. RESOLVE_BASE) CALL Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)
    ! END if NUCLEATION_ROUTINES

    ! =================================================================================================
END if in_turn_acdc




in_turn_any: if (PRC%in_turn(4)) THEN

    ! =================================================================================================
    ! =================================================================================================
    ! AEROSOL PART STARTS HERE
    ! =================================================================================================
    ! =================================================================================================
    ! AEROSOL_ROUTINES: IF (Aerosol_flag.and.(.not. PRC%err)) THEN
    AEROSOL_ROUTINES: IF (Aerosol_flag) THEN

        ! Read in background particles
        IF (N_MODAL>0) THEN
            call Multimodal(MMODES, get_dp(), conc_fit, N_MODAL)
            CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit*dmps_multi)
            do i = 1, n_bins_par
                CALL set_composition(i,nominal_dp(i))
            end do
            if (GTIME%hrs >= DMPS_read_in_time) N_MODAL = -1d0
        END IF

        PARTICLE_INIT: if (use_dmps .and. GTIME%min >= (dmps_ln*dmps_tres_min)) THEN
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
        END IF PARTICLE_INIT

    in_turn_cch_2: if (PRC%in_turn(PRC%cch)) THEN

        ! ..........................................................................................................
        ! ADD NUCLEATED PARTICLES TO PSD, IN THE 1st BIN. This using the mixing method, where two particle distros
        ! are mixed together. NOTE that this step is always done if Aerosol_flag is on, even if ACDC_flag is off. In
        ! case no NPF is wanted, turn off ACDC and make sure NUC_RATE_IN is 0
        dmass = 0d0
        dconc_dep_mix = 0d0

        mix_ratio = -1d0
        dmass(1:size(J_distr,1),VAPOUR_PROP%ind_GENERIC) = nominal_dp(1:size(J_distr,1))**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)
        ! do ii=1,size(J_distr,1)
        !     dconc_dep_mix(ii) = J_TOTAL_M3*GTIME%dt*speed_up(PRC%cch)*J_distr(ii)
        ! end do
        dconc_dep_mix(1:size(J_distr,1)) = J_TOTAL_M3*GTIME%dt*speed_up(PRC%cch)*J_distr

        CALL Mass_Number_Change('mixing')
        ! Negative mixing ratio makes this nucleation

        ! Update current_psd
        current_PSD = new_PSD


        ! ..........................................................................................................
        ! CONDENSATION
        onlyIfCondIsUsed: if (Condensation .and.(.not. PRC%err)) THEN

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

            d_vap = 0
            d_dpar = 0

            CALL Condensation_apc(VAPOUR_PROP,conc_vapour,dmass, GTIME%dt*speed_up(PRC%cch),d_dpar,d_vap)

            ! ERROR HANDLING
            IF (     (maxval(ABS(d_dpar)) > DDIAM_RANGE(2) .or. maxval(abs(d_vap)) > DVAPO_RANGE(2))  &
            .and. use_speed .and. speed_up(PRC%cch)>1) THEN   ! if the changes in diameter are too big
                IF (.not.PRC%err) THEN
                    IF (maxval(ABS(d_dpar)) > DDIAM_RANGE(2)) THEN
                        call SET_ERROR(PRC%cch, 'Too large diameter change: '//f2chr(maxval(abs(d_dpar))))

                    ELSE IF (MAXVAL(abs(d_vap)) > DVAPO_RANGE(2)) THEN
                        call SET_ERROR(PRC%cch,'Too large vapour concentration change: '//f2chr(MAXVAL(abs(d_vap)))//ACHAR(10)&
                        //'Up. lim. of vap: '//f2chr(DVAPO_RANGE(2))//' Troublevapour '//VAPOUR_PROP%vapour_names(i)//f2chr(conc_vapour(i)) )

                    ELSE
                        call SET_ERROR(PRC%cch,'Too large diameter and vapour concentration change: '//f2chr(maxval(abs(d_dpar)))//', '//f2chr(MAXVAL(abs(d_vap))))
                    END IF
                END IF
                dmass = 0.d0

            ELSE IF (.not.PRC%err) THEN ! everything fine -> apply changes
                ! PRC%err = .false.
                ! Calculate growth rates
                IF (GTIME%printnow .and. CALC_GR) THEN
                    CALL PRINT_GROWTH_RATE
                END IF

                ! Distribute mass
                CALL Mass_Number_Change('condensation')

                current_PSD = new_PSD

                ! Update vapour concentrations to chemistry
                CH_GAS(index_cond) = conc_vapour(1:VAPOUR_PROP%n_condorg) *1D-6
                ! Again sulfuric acid always needs special treatment
                CH_GAS(ind_H2SO4) = conc_vapour(n_cond_tot)*1d-6

                ! Check whether timestep can be increased:
                IF (maxval(ABS(d_dpar)) < product(DDIAM_RANGE)**(0.5d0) .and. MAXVAL(abs(d_vap)) < product(DVAPO_RANGE)**(0.5d0) .and. use_speed) THEN
                    if (speed_up(PRC%cch) * 2 * Gtime%dt < speed_dt_limit(PRC%cch+1)) THEN
                        PRC%increase(PRC%cch) = .true.
                    END IF
                END IF
            END IF
        end if onlyIfCondIsUsed
        ! end of condensation
    END IF in_turn_cch_2


    in_turncoa: if (PRC%in_turn(PRC%coa)) THEN

        ! ..........................................................................................................
        ! COAGULATION
        onlyIfCoagIsUsed: if (Coagulation .and.(.not. PRC%err)) then
            ! Solve particle coagulation
            dconc_coag = 0.d0

            Call Coagulation_routine(dconc_coag, GTIME%dt*speed_up(PRC%coa),d_npar)

            ! ERROR HANDLING

            IF ((MINVAL(d_npar) < 0 .or. maxval(ABS(d_npar)) > DPNUM_RANGE(2)) .and. use_speed .and. speed_up(PRC%coa)>1) THEN   ! if the changes in particle numbers are too big
                call set_error(PRC%coa, 'Too large particle number change')
                dconc_coag = 0.d0
            END IF
            ! Distribute mass
            if (.not. PRC%err) Call Mass_Number_Change('coagulation')


            if (.not. PRC%err) THEN                ! Update PSD with new concentrations
                current_PSD = new_PSD
                ! Check whether timestep can be increased:
                IF (maxval(ABS(d_npar)) < product(DPNUM_RANGE)**(0.5d0) .and. use_speed) THEN
                    if (speed_up(PRC%coa) * 2 * Gtime%dt < speed_dt_limit(PRC%coa+1)) THEN
                        PRC%increase(PRC%coa) = .true.
                        ! call increase_speed(PRC%coa)
                    END IF
                END IF
            END IF
        end if onlyIfCoagIsUsed
        ! end of coagulation
    END IF in_turncoa


    in_turn_dep: if (PRC%in_turn(PRC%dep)) THEN
        ! ..........................................................................................................
        ! Deposition
        OnlyIfDepoIsUsed: if ((Deposition .or. LOSSES_FILE /= '') .and.(.not. PRC%err)) then
            ! Solve particle coagulation
            dconc_dep_mix = 0.d0

            IF (LOSSES_FILE /= '') THEN

                ! The losses file is interpolated spatially and temporally to fit the current bin structure and time
                losses_fit = [((INTERP( PAR_LOSSES%sections,                                                            &
                             [((INTERP(PAR_LOSSES%time, PAR_LOSSES%conc_matrix(:,i) )), i=1,size(PAR_LOSSES%time, 1))], &
                             timein=nominal_dp(j))                                                                      &
                             ), j=1,n_bins_par)]

                ! Deposited concentratios calculated here
                dconc_dep_mix = get_conc() * (1 - EXP(-losses_fit*GTIME%dt*speed_up(PRC%dep)))


            ELSE
                call deposition_velocity(nominal_dp,ustar,CHAMBER_CIRCUMFENCE*CHAMBER_HEIGHT,Chamber_floor_area,Chamber_floor_area,&
                                        Chamber_floor_area*CHAMBER_HEIGHT,GTEMPK,GPRES,E_field,dconc_dep_mix,GTIME%dt*speed_up(PRC%dep))
            END IF
            ! ERROR HANDLING; Check whether changes are within limits:

            d_npar = 0.d0
            WHERE (get_conc()>0d0) d_npar = dconc_dep_mix/get_conc()
            IF (maxval(ABS(d_npar)) > DPNUM_RANGE(2) .and. use_speed .and. speed_up(PRC%dep)>1) THEN   ! if the changes in diameter are too big
                call SET_ERROR(PRC%dep,'Too large particle number change')
                dconc_dep_mix = 0.d0
            ELSE  ! everything fine -> apply changes
                ! PRC%err = .false.

                ! Distribute mass
                Call Mass_Number_Change('deposition')
                ! Update PSD with new concentrations
                current_PSD = new_PSD
                !Check whether timestep can be increased:
                IF (maxval(ABS(d_npar)) < product(DPNUM_RANGE)**(0.5d0) .and. use_speed) THEN
                    if (speed_up(PRC%dep) * 2 * Gtime%dt < speed_dt_limit(PRC%dep+1)) THEN
                        PRC%increase(PRC%dep) = .true.
                    END IF
                END IF

            END IF


            ! Update PSD with new concentrations
            current_PSD = new_PSD
        end if OnlyIfDepoIsUsed
    END IF in_turn_dep
        ! end of Deposition
    end if AEROSOL_ROUTINES

END IF in_turn_any
    ! End Aerosol =====================================================================================

    ! SPEED handling
    CHECK_PRECISION: IF (PRC%err) THEN  ! In case of a timestep error (i.e. too large changes in aerosol dynamics)

        ! We need to decrease dt of some process, so we cancel any other speed increases;
        ! they can be done on the next loop if all is okay
        PRC%increase = .false.

        CALL error_handling!(PRCION, speed_up)
        ! print*, 'hyppy tuntemattomaan', add_rounds
        ! Reset the system to the state at previous timestep
        current_PSD = old_PSD
        CH_GAS = CH_GAS_old
        CH_RO2 = CH_RO2_old
        acdc_goback = .true.
        dmps_ln = dmps_ln_old

        ! Reset relative change vectors documenting the precision:
        d_dpar = 0.d0
        d_npar = 0.d0
        d_vap  = 0.d0

        ! If there was no error during the actual timestep
    ELSE

        ! Write printouts to screen
        if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
        if (GTIME%printnow .and. ENABLE_END_FROM_OUTSIDE) CALL CHECK_IF_END_CMD_GIVEN
        ! if (GTIME%printnow) print*, GTIME%dt, speed_up

        ! Save to netcdf
        if (GTIME%savenow) THEN

            ! PSD is also saved to txt
            if (Aerosol_flag) THEN
                WRITE(601,*) GTIME%sec, sum(get_conc()*1d-6), get_conc()*1d-6 / LOG10(bin_ratio)
                WRITE(604,*) GTIME%sec, get_conc()*1d-6
                WRITE(610,*) GTIME%sec, maxval(abs(d_vap)), maxval(abs(d_npar)), maxval(abs(d_dpar))
                FLUSH(610)
                save_measured = conc_fit/dmps_multi
            END IF

            CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3_M3, J_ACDC_DMA_M3, VAPOUR_PROP, save_measured,&
                            1d9*3600/(GTIME%dt*speed_up(PRC%cch))*get_dp()*d_dpar)

        END IF

        ! Add main timestep to GTIME
        if (GTIME%sec + GTIME%dt <= GTIME%SIM_TIME_S) THEN

            add_rounds = minval(pack(speed_up, speed_up > 0))! - modulo(n_of_Rounds,add_rounds)

            GTIME = ADD(GTIME, GTIME%dt*add_rounds)

            n_of_Rounds = n_of_Rounds + add_rounds

        ELSE ! Exit the main loop when time is up
            EXIT
        END IF


    END IF CHECK_PRECISION
    ! SPEED handling


END DO MAINLOOP
!=======================================================================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!       *@@ MAIN LOOP  @@       @@@       #@@@@@@@@@        @@@@@@@@@@             @@@@@       @@@@@@@@@@@@@@@@@       @
!@       ,@ ENDS HERE  @       &@@@       #@@@@@@@*        @@@@@@@                     @@       @@@@@@@@@@@@@@@       @@
!@@       &@@@@@@@@@@@@       &@@@@       #@@@@@@        @@@@@@@          .%@@&/       @@@       @@@@@@@@@@@@@       @@@
!@@@                         #@@@@@       #@@@@%        @@@@@@#        @@@@@@@@@@@@@@@@@@@@                         @@@@
!@@@@                       /@@@@@@       #@@@        @@@@@@@@       #@@@@@@@@@@@@@@@@@@@@@@                       @@@@@
!@@@@@       @@@@@@@@      *@@@@@@@                 @@@@@@@@@        @@@@@@@@@@@@@@@@@@@@@@@@       @@@@@@@.      @@@@@@
!@@@@@@       @@@@@@      ,@@@@@@@@                    ,@@@@@       .@@@@@@@@@@@@@@@@@@@@@@@@@      #@@@@@#      @@@@@@@
!@@@@@@@      @@@@@       @@@@@@@@@       #@@@@@@@       @@@@        @@@@@@@@@@@@@@@@@@@@@@@@@@      @@@@@      @@@@@@@@
!@@@@@@@@      &@@       @@@@@@@@@@       #@@@@@@@@       @@@@       &@@@@@@@@@@@@@@@@@@@@@@@@@%      @@@      %@@@@@@@@
!@@@@@@@@@      @       @@@@@@@@@@@       #@@@@@@@&       @@@@&        @@@@@@@@@@@@@@@@@@@@@@@@@#      @      #@@@@@@@@@
!@@@@@@@@@@            @@@@@@@@@@@@       #@@@@%         &@@@@@@           (%%(       @@@@@@@@@@@#           #@@@@@@@@@@
!@@@@@@@@@@@          @@@@@@@@@@@@@                     @@@@@@@@@@.                    @@@@@@@@@@@/         *@@@@@@@@@@@
!@@@@@@@@@@@@        @@@@@@@@@@@@@@                ,@@@@@@@@@@@@@@@@@@,            @@@@@@@@@@@@@@@@.       ,@@@@@@@@@@@@
!=======================================================================================================================

CALL PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

if (Aerosol_flag) THEN
    CLOSE(601)
    CLOSE(604)
    CLOSE(610)
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
SUBROUTINE error_handling
    IMPLICIT none

    ! Write some error information to file
    CALL OUTPUTBOTH(608,FMT_NOTE0,'Precision tolerance exceeded at time '//GTIME%hms//' ('//TRIM(f2chr((GTIME%sec)))//' sec)')
    CALL OUTPUTBOTH(608,FMT_SUB,  'Process exceeding tolerance: '//TRIM(PRC%pr_name(PRC%proc)))
    CALL OUTPUTBOTH(608,FMT_SUB,  'Message from process: '//trim(PRC%err_text))
    CALL OUTPUTBOTH(608,FMT_SUB,  'Current simulation and process timestep [s]'//f2chr(GTIME%dt)//', '//TRIM(f2chr(GTIME%dt*speed_up(PRC%proc))) )

    ! reduce the speed_up factor or, if necessary, the integration time step:
    IF (speed_up(PRC%proc) > 1) THEN
        speed_up(PRC%proc) = speed_up(PRC%proc) / 2
        CALL OUTPUTBOTH(608,FMT_SUB, '  => reduce speed_up:'//di2chr(speed_up(PRC%proc)*2)//'->'//di2chr(speed_up(PRC%proc)))
        errortime(PRC%proc) = GTIME%sec
    ELSE
        CALL OUTPUTBOTH(608,FMT_NOTE0, '  => SHOULD REDUCE DT')
    END IF

    CALL OUTPUTBOTH(608, FMT_MSG,'Timesteps adjusted.')
    print FMT_LEND,

    flush(608)
    print*,''

    ! RESET ERROR
    PRC%err = .false.
    PRC%proc = 0
    PRC%err_text = ''
    ! add_rounds = minval(pack(speed_up, speed_up > 0))
    write(608,*) ''


END SUBROUTINE error_handling


! Routine to increase speed so that unused speed factors stay as
! the highest number and are not hindering main timestep changes
Subroutine INCREASE_SPEED(ii)
    IMPLICIT NONE
    INTEGER :: ii

    ! do ii=1,SIZE(PRC%pr_name)
        if (PRC%increase(ii) .and. GTIME%sec-errortime(ii)>60d0) THEN
            speed_up(ii) = speed_up(ii) * 2
            CALL OUTPUTBOTH(608,'*','  '//GTIME%hms//' --> Speeding '//TRIM(PRC%pr_name(ii))//' to '//TRIM(di2chr(speed_up(ii))) )
            print*, ''
        END IF
    ! end do
    PRC%increase(ii) = .false.

    if ((.not. Deposition).and.(LOSSES_FILE == '')) speed_up(PRC%dep) = -999_dint
    if (.not. Coagulation) speed_up(PRC%coa) = -999_dint

END Subroutine INCREASE_SPEED


! =================================================================================================
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

SUBROUTINE ACDC_J(C)
    implicit none
    REAL(dp) :: c_org = 0d0
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
        CALL get_acdc_J(H2SO4,NH3,c_org,GCS,C(inm_TEMPK),IPR,GTIME,&
            ss_handle,J_ACDC_NH3_M3,acdc_cluster_diam, J_NH3_BY_IONS, GTIME%dt*speed_up(PRC%cch), acdc_goback)
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
    END IF
    ! Speed up program by ignoring nucleation when there is none
    if ((DMA > 1d6 .and. H2SO4>1d9) .or. (.not. skip_acdc)) THEN
        CALL get_acdc_D(H2SO4,DMA,c_org,GCS,C(inm_TEMPK),GTIME,ss_handle,J_ACDC_DMA_M3,acdc_cluster_diam, GTIME%dt*speed_up(PRC%cch), acdc_goback)
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'DMA IGNORED'
    END IF
    acdc_goback = .false.
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
        OPEN(UNIT=609, FILE='ModelLib/required/nucl_homs.txt', STATUS='OLD', ACTION='READ', iostat=jj)
        n = rowcount(609)
        print FMT_MSG, 'Using parametrisation for organic nucleation with '//i2chr(n)//' nucleating compounds'
        allocate(inds(n))
        inds = 0
        DO i=1,n
            read(609, *) name
            if (IndexFromName(TRIM(name), SPC_NAMES)>0) inds(i) = IndexFromName(TRIM(name), SPC_NAMES)
        END DO
        if (PRODUCT(inds) == 0) THEN
            print FMT_LEND,
            print FMT_FAT0, "List of nucleating organic compounds has gases which are not in the chemistry."
            print FMT_SUB, "Options: Edit the file 'ModelLib/required/nucl_homs.txt' or turn of organic nucleation."
            print FMT_SUB, "Bye."
            print FMT_LEND,
            STOP
        END IF
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
    DO j=1,INT(GTIME%dt*speed_up(PRC%cch)/0.01D0)
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

    if (Aerosol_flag) print FMT_MSG, 'Max change in vapours '&
            //TRIM(f2chr(1d2*d_vap(maxloc(d_vap,1))))//'% for '//VAPOUR_PROP%vapour_names(maxloc((d_vap)))
    if (Aerosol_flag) print FMT_MSG, 'Max change in par conc. '&
            //TRIM(f2chr(1d2*d_npar(maxloc(d_npar,1))))//'% in bin # '//i2chr((maxloc(d_npar,1)))
    if (Aerosol_flag) print FMT_MSG, 'Max change in par diam. '&
            //TRIM(f2chr(1d2*d_dpar(maxloc(d_dpar,1))))//'% in bin # '//i2chr((maxloc(d_dpar,1)))

    print FMT10_2CVU,'ACID C: ', C(inm_H2SO4), ' [1/cm3]','sum(An)/A1',clusteracid,' []'
    print FMT10_3CVU,'Temp:', C(inm_TempK), ' [K]','Pressure: ', C(inm_pres), ' [Pa]', 'Air_conc', C_AIR_cc(C(inm_TempK), C(inm_pres)), ' [1/cm3]'
    IF (inm_NH3   /= 0) print FMT10_3CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3_M3*1d-6, ' [1/s/cm3]','sum(Nn)/N1',clusterbase,' []'
    IF (inm_DMA   /= 0) print FMT10_3CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA_M3*1d-6, ' [1/s/cm3]', 'J-total', J_TOTAL_M3*1e-6,' [1/s/cm3]'
    print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', GCS , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
    if ((GTIME%sec)>0 .and. (cpu2 - cpu1 > 0d0)) print '("| ",a,i0,a,i0.2,a,i0,a,i0.2,t65,a,f7.2,t100,"|")', 'Elapsed time (m:s) ',int(cpu2 - cpu1)/60,':',modulo(int(cpu2 - cpu1),60) ,' Est. time to finish (m:s) ',&
                            int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec))/60,':', MODULO(int((cpu2 - cpu1)/((GTIME%sec))*(GTIME%SIM_TIME_S-GTIME%sec)),60),&
                            'Realtime/Modeltime: ', (GTIME%sec-start_time_s)/(cpu2 - cpu1)
    print '("| ",a,3(", ",f6.2)," (",3(" ",i0),")"t100,"|")', 'Time steps (s) and multipliers: '//TRIM(f2chr(GTIME%dt))//' (main)', GTIME%dt*speed_up, speed_up
    print FMT_LOOPEND,

END SUBROUTINE PRINT_KEY_INFORMATION


! =================================================================================================
! If user had opted a save or print interval which does not concide with last timestep, we print
! and save the values at last timesteps
! =================================================================================================
SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY
    implicit none
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
                        1d9*3600/(GTIME%dt*speed_up(PRC%cch))*get_dp()*d_dpar)


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

    close(610)

END SUBROUTINE FINISH


SUBROUTINE PRINT_GROWTH_RATE
    IMPLICIT NONE
    CHARACTER(len=60) :: fmt
    CHARACTER(len=90) :: fmt2
    INTEGER :: nb
    nb = size(GGR)

    ii = MINLOC(ABS(nominal_dp-GR_bins(1)),1)
    GGR(1) = 1d9 * 3600/(GTIME%dt*speed_up(PRC%cch)) * (sum((nominal_dp(1:ii) * d_dpar(1:ii)))/ii)
    if (nb>1) THEN
        DO i=2,size(GGR)
            jj = MINLOC(ABS(nominal_dp-GR_bins(i)),1)
            GGR(i) = 1d9 * 3600/(GTIME%dt*speed_up(PRC%cch)) * (sum((nominal_dp(ii:jj) * d_dpar(ii:jj)))/(jj-ii))
            ii = jj
        END DO
    END IF
    print FMT_MSG, 'Growth rates (nm/h) between sizes (nm, starting from smallest bin)'
    fmt2 = '("| ",a,f6.1,t18,'//TRIM(i2chr(MIN(SIZE(GGR,1), 3)))//'("   -->",f6.1,),a,t100,"|")'
    fmt  = '("| ",a,t18,'//TRIM(i2chr(MIN(SIZE(GGR,1), 3)))//'(es9.2,3(" ")),a,t100,"|")'

    print fmt2, 'Sizes:  ',min_particle_diam*1d9, GR_bins(1:MIN(SIZE(GGR,1), 3))*1d9, ' [nm]'
    print fmt, 'GR:', GGR(1:MIN(SIZE(GGR,1), 3)), ' [nm/h]'

END SUBROUTINE PRINT_GROWTH_RATE

SUBROUTINE LINK_VARIABLES
    IMPLICIT NONE
    INTEGER :: i
    print FMT_HDR, 'Checking for linked variables'
    DO i=1,N_VARS
        IF ((MODS(i)%ISPROVIDED).and.(MODS(i)%TIED /= '')) THEN
            if (I < 3) THEN
                PRINT FMT_FAT0, "Tying temperature or pressure to other variables is not possible."
                stop
            END IF
            ii = IndexFromName(MODS(i)%TIED)
            IF (ii == 0) THEN
                print FMT_FAT0, TRIM(MODS(i)%NAME)//' is tied to "'//TRIM(MODS(i)%TIED)//'" which is unavailable.'
                stop
            ELSE
                print FMT_MSG, 'Tying '//TRIM(MODS(i)%NAME)//' to '//TRIM(MODS(i)%TIED)
                INDRELAY_TIED(I) = ii
            END IF
        END IF
    END DO
    print FMT_LEND,

END SUBROUTINE LINK_VARIABLES

SUBROUTINE OUTPUTBOTH(unit, format, text)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: unit
    CHARACTER(len=*), INTENT(in) :: format
    CHARACTER(len=*), INTENT(in) :: text

    if (FORMAT == '*') THEN
        write(*,*) text
        write(unit,*) text
    ELSE
        write(*,format) text
        write(unit,Format) text
    END IF

END SUBROUTINE OUTPUTBOTH

! This subroutine makes it possible to end a simulation from outside using a file called ENDNOW.INIT
SUBROUTINE CHECK_IF_END_CMD_GIVEN
    implicit none
    integer :: ioi
    CHARACTER(len=256) :: command, date, time
    OPEN(620, file=RUN_OUTPUT_DIR//"/ENDNOW.INIT", iostat=ioi, STATUS='OLD',action='read')
    if (ioi == 0) THEN
        read(620,*,iostat=ioi) command
        close(620)
        if ((ioi == 0).and.(TRIM(command)=='STOP')) THEN
            print FMT_MSG, 'ENDNOW.INIT WAS USED, STOPPING THE SIMULATION FROM OUTSIDE. GOOD BYE.'
            GTIME%SIM_TIME_S = GTIME%sec
            call date_and_time(date,time)
            OPEN(620, file=RUN_OUTPUT_DIR//"/ENDNOW.INIT", iostat=ioi, STATUS='REPLACE',action='WRITE')
            write(620,*) 'Run was terminated by the user at simulation time '//GTIME%hms//'. Real date is: '//TRIM(date)//', time is '//TRIM(time(1:6))
            close(620)
            ioi = RENAME(RUN_OUTPUT_DIR//"/ENDNOW.INIT", RUN_OUTPUT_DIR//"/ENDNOW._OLD")
            if (ioi /= 0) print*, 'Error while renaming the end cmd file.'
            !
        END IF
    END IF

END SUBROUTINE CHECK_IF_END_CMD_GIVEN

END PROGRAM ARCA_main
