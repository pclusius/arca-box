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

PROGRAM ARCA_main

USE second_MAIN                            ! Main second file
USE second_PARAMETERS                      ! CH_NSPEC (originally NSPEC) and chemical indices, ind_xxxx, come from here
USE second_Precision,  ONLY : dp           ! KPP Numerical type
USE second_Monitor,    ONLY : SPC_NAMES    ! Names of chemicals from KPP
USE SECOND_REACTIVITY, ONLY : NREACTIVITY  !
USE Chemistry
USE constants
USE AUXILLARIES
USE INPUT
use OUTPUT
USE PSD_scheme
USE aerosol_dynamics
USE custom_functions

IMPLICIT NONE

REAL(dp), ALLOCATABLE  :: conc_pp(:,:)      ! [#/m^3] Particle phase concentrations, DIM(n_bins_par,n_cond_tot)
real(dp) :: T_prev = 0d0
real(dp) :: P_prev = 0d0

! ==================================================================================================================
! Note about file unit handles: Numbers between 100-499 are not used anywhere in the model, so use these in you need
! to open additional files. Generally 1-99 can be used temporarily in KPP, 600-series is used in ARCA_main.f90,
! 700-series is used in output.f90, 800-series is used in input.f90, 900-series in Chemistry.f90
! ==================================================================================================================

! ==================================================================================================================
! VARIABLE DECLARATION. MOST OF THE GLOBAL VARIABLES ARE DEFINED IN INPUT.F90 and CONSTANTS.F90
! ==================================================================================================================

REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)    ! Array to hold input variables for the current timestep
REAL(DP), ALLOCATABLE :: CH_GAS(:), CH_GAS_old(:),d_chem(:), reactivities(:) ! Array to hold all chemistry compounds; old: to restore in case of an error related to timestep handling
REAL(dp), ALLOCATABLE :: conc_fit(:)      ! An array giving particle conc independant of PSD_style [m⁻³]
REAL(dp), ALLOCATABLE :: losses_fit(:)    ! interpolated loss rates [/s]
REAL(dp), ALLOCATABLE :: losses_fit0(:),losses_fit1(:), intrp_losses(:), inv_loss(:) ! Experimental stuff, deprecated [m-3]
REAL(dp), ALLOCATABLE :: save_measured(:)      ! An array for saving measurements

INTEGER               :: dmps_ln = 0, dmps_ln_old      ! line number from where background particles are read from
INTEGER               :: dmps_sp_min = 0, dmps_sp_max = 0 ! Indices for dmps_special

! Variables related to chemistry module
REAL(dp) :: CH_RO2,CH_RO2_old       ! RO2 concentration in [molecules / cm^3]
REAL(dp) :: CH_H2O                  ! H20 concentration in [molecules / cm^3]
REAL(dp) :: CH_Beta                 ! solar zenit angle
REAL(dp) :: EW                      ! Water content in Pa
REAL(dp) :: ES                      ! Saturation vapour pressure
REAL(dp) :: STW                     ! Water surface energy density / surface tension

! Variables related to aerosol module
REAL(dp), dimension(:), allocatable:: conc_vapour

! Transient variables
CHARACTER(:), allocatable:: RUN_OUTPUT_DIR      ! Saves the output directory relative to this executable
CHARACTER(1000) :: inibuf                       ! Buffer to save backp from the INITfile that was called
INTEGER         :: I,II,J,JJ                    ! Loop indices
INTEGER         :: ioi                          ! iostat variable
INTEGER(dint)   :: will_stall=0                 ! transient variable
INTEGER(dint)   :: add_rounds=1                 ! transient variable
INTEGER(dint)   :: bookkeeping(3,2) = 0_dint    ! save total rounds and latest point of error
REAL(dp)        :: cpu1,cpu2,cpu3,cpu4          ! CPU times in seconds
REAL(dp)        :: errortime(3) = -9999d0       ! CPU time in seconds
REAL(dp)        :: H2SO4                        ! Transient keeper for H2SO4 variable
LOGICAL         :: acdc_goback = .false.        ! Used for time step optimization

REAL(dp), ALLOCATABLE   :: nominal_dp(:)        ! [m] array with nominal diameters. Stays constant independent of PSD_style
INTEGER(dint)           :: n_of_Rounds = 0_dint ! Stores the total rounds, most of whic are skipped if OPTIMIZE_DT
REAL(dp), ALLOCATABLE   :: J_distr(:)           ! Transient variable for distributing ACDC particles
REAL(dp), ALLOCATABLE   :: d_dpar(:)            ! array reporting relative changes to the diameter array within a single timestep
REAL(dp), ALLOCATABLE   :: d_npar(:)            ! array reporting relative changes to the particle number array within a single timestep
REAL(dp), ALLOCATABLE   :: d_npdep(:)           ! array reporting relative changes to the particle number array within a single timestep
REAL(dp), ALLOCATABLE   :: d_vap(:)             ! array reporting relative changes to the vapour concentration array within a single timestep
REAL(dp)                :: refdvap  = 0d0       ! holds the reference maximum of d_vap, which is the average of three largest or three
integer                 :: irefdvap = 1         !   :... smallest (if negative) in d_vap, whichever is larger in magnitude. irefvap is the index of
logical                 :: logdvap  = .false.   !   :... the chosen maximum or minimum value. logdvap is true if abs(refdvap) is over the precision limits
logical                 :: update_K = .true.    ! is coagulation coefficient to be updated? In the main loop TRUE only if Temp or Press differs enough from the prev time step.
REAL(dp) :: A_vert                              ! chamber vertical wall area [m²]
REAL(dp) :: Vol_chamber                         ! chamber volume [m³]
REAL(dp) :: E_field = 0d0                       ! chamber Electric field (not implemented yet) [V/m]
! REAL(dp) :: wl_rates(2)                         ! Transient vector for storing reversible chemical loss rates for chamber [1/s]
CHARACTER(len=64) :: CurrentChem, CurrentVers ! Name of the chemistry module and current compiled version


! This block is handled by C preprocessor --------------------------!
#ifdef CHEM
CurrentChem = CHEM
#endif

! This block is handled by C preprocessor --------------------------!
#ifdef VERSION
CurrentVers = VERSION
#endif
! end of C preprocessor commands -----------------------------------!

! Welcome message
print'(a,t35,a)', achar(10),  '--~:| ARCA BOX MODEL '//TRIM(CurrentVers)//' |:~--'//achar(10)
print FMT_HDR, 'Compiled with "'//TRIM(CurrentChem)//'" chemistry module'
print*, ''
! ==================================================================================================================

! Declare most variables and read user input and options in input.f90
CALL READ_INPUT_DATA
! ==================================================================================================================

! ==================================================================================================================
! If some variables were linked together, they will be handled here
CALL LINK_VARIABLES
! ==================================================================================================================

! ==================================================================================================================
IF (Chemistry_flag.or.Condensation) THEN
    ! Check that the input exists in chemistry, or if not, print warning
    CALL CHECK_INPUT_AGAINST_KPP
    ! This only called once for KPP in the beginning
    call READ_rates(CurrentChem)
    IF (Chemistry_flag) CALL KPP_SetUp(factorsForReactionRates)
ENDIF
! ==================================================================================================================

! ==================================================================================================================
! Calculate chamber dimensions, assume square floor shape.
! Considering the simplicity of the parametrization, this is good enough. ==========================================
A_vert = 4d0*SQRT(CHAMBER_FLOOR_AREA)*CHAMBER_HEIGHT
Vol_chamber = CHAMBER_FLOOR_AREA*CHAMBER_HEIGHT
if (((Deposition.and.LOSSES_FILE=='').or.Chem_Deposition).and.(EQUAL(Vol_chamber, 0d0))) &
STOP '   ----------- CHAMBER VOLUME IS ZERO! -------------'
! ==================================================================================================================

! ==================================================================================================================
! If particles are considered -> initialize a particle representation, set initial PSD and determine composition
INIT_AEROSOL: IF (Aerosol_flag) THEN

    ! Initialize the Particle representation
    CALL INITIALIZE_PSD

    ! Initialize the nominal diameter vector
    ALLOCATE(nominal_dp(n_bins_par))
    nominal_dp = get_dp()

    ! J_Distr is a vector containing the weighing factor for newly nucleated particles. Majority (~75-90%) of new clustes
    ! are put in the first model bin (whatever that is), and the rest are distributed in the next couple of bins in
    ! diminishing numbers (exp decay), but not further than NPF_DIST (default=1.15) times the smallest diameter.
    ALLOCATE(J_distr(   MINLOC(abs(nominal_dp-nominal_dp(1)*NPF_DIST),1)   ))
    do ii=1,size(J_distr)
        J_distr(ii) = 0.5*exp(-0.7d0*ii)
    end do
    J_distr(1) = J_distr(1) + 1-sum(J_distr)
    print FMT_MSG, 'Distributing nucleated particles over '//TRIM(i2chr(size(J_distr, 1)))//' bins.'

    ! Initialize the dummy for dmps fitting
    ALLOCATE(conc_fit(n_bins_par))
    save_measured = conc_fit
    losses_fit    = conc_fit

    ! Allocate the change vectors for integration timestep control
    ALLOCATE(d_dpar(n_bins_par))
    ALLOCATE(d_npar(n_bins_par))
    ALLOCATE(d_npdep(n_bins_par))
    d_dpar = 0d0
    d_npar = 0d0
    d_npdep = 0d0

    if (use_dmps_partial .and. use_dmps) THEN
        ! write(*, FMT_MSG) 'Using dmps_special:'

        ! both limits have to be far enough from the smallest and largest diameters
        if (dmps_lowband_upper_limit > nominal_dp(2)) THEN
            dmps_sp_min = minloc(abs(get_dp()-dmps_lowband_upper_limit),1)
            print FMT_SUB, 'Lower bins replaced from 1 to '//i2chr(dmps_sp_min)
        ELSE IF (dmps_lowband_upper_limit > 0d0) THEN
            STOP '     ---- dmps_lowband_upper_limit is is too small and will not be used ----'
        END IF

        if (dmps_highband_lower_limit < nominal_dp(n_bins_par-1) &
         .and. dmps_highband_lower_limit > nominal_dp(2)) THEN
            dmps_sp_max = minloc(abs(get_dp()-dmps_highband_lower_limit),1)
            print FMT_SUB, 'Upper bins replaced from indices '//TRIM(i2chr(dmps_sp_max))//' to '//TRIM(i2chr(n_bins_par))
        ELSE IF (dmps_highband_lower_limit > 0d0) THEN
            STOP '     ---- dmps_highband_lower_limit is is too large and will not be used ----'
        END IF

    END IF


    ALLOCATE(conc_vapour(n_cond_tot))
    ALLOCATE(CONDENSED_MASS(n_cond_tot))
    CONDENSED_MASS = 0d0

    ! Allocate the change vectors for integration timestep control
    ALLOCATE(d_vap(max(VAPOUR_PROP%n_cond_tot,1)))
    d_vap = 0d0

    if (Deposition) THEN
        ALLOCATE(Depos_composition(max(VAPOUR_PROP%n_cond_tot,1)))
        Depos_composition = 0d0
    END if
    if (Chem_Deposition) THEN
        ALLOCATE(c_org_wall(VAPOUR_PROP%n_cond_org))
        c_org_wall = 0d0
        c_org_wall_old = c_org_wall
        ALLOCATE(VAP_DEP_MASS_WALLS(n_cond_org))
        VAP_DEP_MASS_WALLS = 0d0
    END IF

    do ii=1,n_bins_par
        if (Kelvin_taylor) THEN
            pre_Kelvin(ii,:) = 2D0*VAPOUR_PROP%surf_tension * VAPOUR_PROP%molar_mass / ( Rg * 300d0 * VAPOUR_PROP%density * nominal_dp(ii)/2d0)
        else
            pre_Kelvin(ii,:) = EXP(4D0*VAPOUR_PROP%surf_tension * VAPOUR_PROP%molar_mass / ( Rg * 300d0 * VAPOUR_PROP%density * nominal_dp(ii)) )
        end if

    end do

    print FMT_LEND,

else ! This is run even if aerosol module is off

    ALLOCATE(conc_vapour(1))
    ALLOCATE(d_vap(1))
    ALLOCATE(d_dpar(1))
    ALLOCATE(d_npar(1))
    ALLOCATE(d_npdep(1))
    d_dpar  = 0d0
    d_npar  = 0d0
    d_npdep = 0d0

END IF INIT_AEROSOL

! If Condensation is off, the changes is gases is monitored from chemistry
if (Chemistry_flag.and..not.Condensation) ALLOCATE(d_chem(NSPEC))

! Allocate and initiate the timestep vector for interpolated input
ALLOCATE(TSTEP_CONC(N_VARS))
TSTEP_CONC = 0

! Allocate and initiate the timestep vector for chemical compounds
ALLOCATE(CH_GAS(NSPEC))
ALLOCATE(reactivities(NREACTIVITY))
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

    ! Save all condensibles for easier access
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
open(unit=607, file=RUN_OUTPUT_DIR//'/InitBackup_temp.txt')
DO while (ioi == 0)
    read(606,'(a)', iostat = ioi) inibuf
    write(607,'(a)') TRIM(inibuf)
end do
close(606)
close(607)
call system('mv '//RUN_OUTPUT_DIR//'/InitBackup_temp.txt '//RUN_OUTPUT_DIR//'/InitBackup.txt', I)
if (I/=0) print FMT_WARN0, "Could not save InitBackup.txt, is the file locked?"

! Open netCDF files
CALL OPEN_FILES( RUN_OUTPUT_DIR, Description,CurrentChem,CurrentVers, MODS, CH_GAS, reactivities, VAPOUR_PROP,Vol_chamber)

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

open(unit=610, file=RUN_OUTPUT_DIR//'/Changes.txt',status='replace',action='write')

if (OPTIMIZE_DT) THEN

    write(610,'("# ddiam, dpnum, dvap: ", 3(es10.3," ",es10.3))') DDIAM_RANGE*1d-2,DPNUM_RANGE*1d-2,DVAPO_RANGE*1d-2
    print FMT_HDR, 'Simulation time step will be optimized for precision'
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for particle diameter',      DDIAM_RANGE
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for particle concentration', DPNUM_RANGE
    print'("| ",a,": ",t47,2(f7.2, " % "),t100,"|")', 'precision limits for vapour concentration',   DVAPO_RANGE

    ! Input was in percetages, here we change them to fractional
    DDIAM_RANGE = DDIAM_RANGE * 1d-2
    DPNUM_RANGE = DPNUM_RANGE * 1d-2
    DVAPO_RANGE = DVAPO_RANGE * 1d-2
else
    write(610,'("# ddiam, dpnum, dvap: ", 3(es10.3," ",es10.3))') 0d0,0d0,0d0,0d0,0d0,0d0

    print FMT_HDR, 'Beginning simulation with constant timestep'
    write(608, *) 'Time step not optimized, dt: ',GTIME%dt
    flush(608)
end if

write(610,'(4(a,"                "),a)') '#  time_sec  ','max_d_diam','max_d_npar','max_d_vap ','max_d_npdep'

print*, ''
print FMT_LEND
print FMT_MSG, 'Starting main loop'

if (ENABLE_END_FROM_OUTSIDE) THEN
    print FMT_SUB, "Simulation can be stopped using file called ENDNOW.INIT with text 'STOP' in it, saved"
    print FMT_SUB, "in the output directory. The program stops in orderly fashion upon next PRINTNOW time"
    print FMT_SUB, "and closes the output files."
end if

print FMT_LEND

if (start_time_s>0) then
    write(*, '(a)', advance='no') 'Starting simulation at: ', ACHAR(10)
    print FMT_TIME, GTIME%hms
end if

! Sanity checks for bad input, add more when needed
if (EQUAL(Cw_eqv,0d0)) STOP 'Effective wall concentration should not be zero.'
! End sanity checks, entering bonkers mode

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
    if (OPTIMIZE_DT) THEN
        PRC%in_turn(PRC%cch) = MODULO(n_of_Rounds,speed_up(PRC%cch)) == 0_dint.and.(Condensation.or.Chemistry_flag)
        PRC%in_turn(PRC%coa) = MODULO(n_of_Rounds,speed_up(PRC%coa)) == 0_dint.and.Coagulation
        PRC%in_turn(PRC%dep) = MODULO(n_of_Rounds,speed_up(PRC%dep)) == 0_dint.and.Deposition
        PRC%in_turn(4)       = (PRC%in_turn(PRC%cch).or.(PRC%in_turn(PRC%coa).or.PRC%in_turn(PRC%dep)))

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
        if (Chem_Deposition) c_org_wall_old = c_org_wall
        ! =================================================================================================
    ELSE
        print FMT_INTRMT, 'Entered an unnecessary loop, should not have happened.'
    END IF


    ! Assign values to input variables. N_VARS will cycle through all variables that user can provide
    ! or alter, and value leave zero if no input was provided.
    ! GC_AIR_NOW, TEMPK, PRES and RH are saved as global variables and are available everywhere.
    ! H2SO4 CS is based on the input, if provided, otherwise it is calculated from the aerosol population
    ! USE GC_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3

if (PRC%in_turn(4)) THEN
    DO I = 1, N_VARS
        IF (MODS(i)%ISPROVIDED) TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I),unit=FILE_TIME_UNIT) .mod. MODS(I)
        IF (I == 2) THEN
            GTEMPK     = TSTEP_CONC(inm_TempK)
            GPRES      = TSTEP_CONC(inm_pres)
            GC_AIR_NOW = C_AIR_cc(GTEMPK, GPRES)
        END IF
    END DO
    ! Here the tied variables are handled.
    DO I = 1, N_VARS
        IF ((MODS(i)%ISPROVIDED).and.(MODS(i)%TIED /= '')) THEN
            TSTEP_CONC(I) = TSTEP_CONC(INDRELAY_TIED(I)) * MODS(I)%MULTI + UCONV(MODS(I)%SHIFT, MODS(I))
        END IF
    END DO

    GRH = TSTEP_CONC(inm_RH)
    H2SO4 = TSTEP_CONC(inm_H2SO4)

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


    CHEMISTRY_ROUTINES: IF (Chemistry_flag .or. CONDENSATION) THEN
        ! Calculate Water vapour pressure and concentration
        CALL WATER(ES,EW,CH_H2O,STW)

        if (swr_is_time_dependent) THEN
            swr_spectrum = [(INTERP(swr_times, swr_temporal_data(:,ii)), ii=1,84)]
        end if

        if (GTIME%hrs<=FLOAT_CHEMISTRY_AFTER_HRS) THEN
            DO I = 1, N_VARS ! <-- N_VARS will cycle through all input variables
                ! This will pick those that were paired in CHECK_INPUT_AGAINST_KPP
                IF (INDRELAY_CH(I)>0) THEN
                    ! In concentrations
                    IF (I<=N_CONCENTRATIONS.and.GTIME%hrs<=FLOAT_CONC_AFTER_HRS) THEN
                        CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                    ! In emissions
                    ELSE IF (I>N_CONCENTRATIONS.and.GTIME%hrs<=FLOAT_EMIS_AFTER_HRS) THEN
                        CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                    END IF
                END IF
            END DO
        END IF

        ! Determine the SWR direction. GTIME%JD==0 means not date is set (using index, like in chamber)
        IF (GTIME%JD==0) THEN
            CH_Beta = 90d0
        ELSE
            ! Solar angle above horizon. For this to properly work, lat, lon and Date need to be defined in INIT_FILE
            call BETA(CH_Beta)
        END IF

        IF (Chemistry_flag) Call CHEMCALC(CH_GAS, GTIME%sec, (GTIME%sec + GTIME%dt*speed_up(PRC%cch)), GTEMPK, max(0d0,TSTEP_CONC(inm_swr)),&
                    CH_Beta,CH_H2O, GC_AIR_NOW, GCS, TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2, reactivities, swr_spectrum, SWR_IS_ACTINICFLUX)

        if ( ( minval(CH_GAS)<-1d2 ).and.GTIME%sec>100d0) print FMT_WARN0,&
        'Negative values from chemistry, setting to zero: '//SPC_NAMES(MINLOC(CH_GAS))//', '//TRIM(f2chr(MINVAL(CH_GAS)))
        WHERE (CH_GAS<0d0) CH_GAS = 0d0

        ! Normally changes in chemical compounds is monitored from condensation scheme
        if (.not. Condensation) THEN
            d_chem = 0d0
            where (CH_GAS_OLD>MIN_CONCTOT_CC_FOR_DVAP) d_chem = ((CH_GAS-CH_GAS_OLD)/CH_GAS_OLD)

            if (OPTIMIZE_DT) THEN
                if (maxval(abs(d_chem)) > DVAPO_RANGE(2)) THEN
                    CALL SET_ERROR(PRC%cch, 'Too large change in chemistry: '//f2chr(1d2* d_chem(maxloc(abs(d_chem),1)))//'%' )
                else if (maxval(abs(d_chem)) < DVAPO_RANGE(1)) THEN
                    if (speed_up(PRC%cch) * 2 * Gtime%dt <= DT_UPPER_LIMIT(PRC%cch)) THEN
                        PRC%increase(PRC%cch) = .true.
                    END IF
                END IF
            END IF
        END IF


        if (H2SO4_ind_in_chemistry>0) H2SO4 = CH_GAS(H2SO4_ind_in_chemistry)

        IF (AFTER_CHEM_ON) CALL AFTER_CHEM(TSTEP_CONC,CH_GAS,CH_GAS_old,CH_RO2,CH_Beta,CH_H2O,CH_RO2_old,J_TOTAL_M3)

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
    ! NOTE: Since ARCA 1.2, RESOLVE BASE has been disabled
    ! =================================================================================================

    in_turn_acdc: if (PRC%in_turn(4).and..not.PRC%err) THEN
        if (ACDC) THEN
            CALL ACDC_J(TSTEP_CONC, GTIME%dt*minval(pack(speed_up, speed_up > 0)))
            J_TOTAL_M3 = sum(G_ACDC(:)%J_OUT_M3(1)) ! [particles/s/m^3]
            J_TOTAL_M3 = J_TOTAL_M3 + TSTEP_CONC(inm_JIN) * 1d6 ! [particles/s/m^3]
        else
          ! Only use input formation rate
          J_TOTAL_M3 = TSTEP_CONC(inm_JIN) * 1d6
        END if

        IF (ORG_NUCL) CALL ORGANIC_NUCL(J_TOTAL_M3)

        IF (AFTER_NUCL_ON) CALL AFTER_NUCL(TSTEP_CONC,CH_GAS,J_TOTAL_M3)
        ! if (GTIME%savenow .and. RESOLVE_BASE .and..not.OPTIMIZE_DT) CALL Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)

    ! =================================================================================================
    END if in_turn_acdc


    in_turn_any: if (PRC%in_turn(4).and..not.PRC%err) THEN

    ! =================================================================================================
    ! =================================================================================================
    ! AEROSOL PART STARTS HERE
    ! =================================================================================================
    ! =================================================================================================

    AEROSOL_ROUTINES: IF (Aerosol_flag) THEN

        ! Read in background particles
        IF (N_MODAL>0) THEN
            call Multimodal(MMODES, get_dp(), conc_fit, N_MODAL)
            conc_fit = conc_fit*1d6
            CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit)
            do i = 1, n_bins_par
                if (GTIME%sec<1d0) THEN
                    CALL set_composition(i,nominal_dp(i), .false.)
                else
                    CALL set_composition(i,nominal_dp(i), Use_old_composition)
                end if
            end do
            if (GTIME%hrs >= DMPS_read_in_time) N_MODAL = -1d0
        END IF

        PARTICLE_INIT: if (use_dmps .and. GTIME%min >= (dmps_ln*dmps_tres_min)) THEN
            ! Particles are read in from measuremensts throughout the simulation and saved to Particles.nc for comparison
            CALL GeneratePSDfromInput( BG_PAR%sections,  BG_PAR%conc_matrix(min(dmps_ln+1, size(BG_PAR%time, 1)),:), conc_fit )
            conc_fit = conc_fit*dmps_multi

            ! If the mdodel is still in initialization mode, replace the model particles with measured
            if (GTIME%hrs <= DMPS_read_in_time) THEN
                ! Print user info
                if (gtime%printnow) print*, ''
                if (gtime%printnow) print FMT_INTRMT, 'Replacing PSD with background particles from '//TRIM(i2chr(dmps_ln+1))&
                                                    //'. measurement (row '//TRIM(i2chr(dmps_ln+2))//' in file)'
                CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit)
                do i = 1, n_bins_par
                    if (GTIME%sec<1d0) THEN
                        CALL set_composition(i,nominal_dp(i), .false.)
                    else
                        CALL set_composition(i,nominal_dp(i), Use_old_composition)
                    end if
                end do
            END IF

            ! If the mdodel is still in initialization mode, replace the model particles with measured
            if (use_dmps_partial &
            .and. (GTIME%hrs > DMPS_read_in_time)&
            .and. (GTIME%hrs < END_DMPS_PARTIAL) &
            ) THEN
                ! Print user info
                if (gtime%printnow) print*, ''
                if (gtime%printnow) print FMT_INTRMT, 'Replacing PSD partially with background particles from '//TRIM(i2chr(dmps_ln+1))&
                                                    //'. measurement (row '//TRIM(i2chr(dmps_ln+2))//' in file)'
                ! if small particles are replaced with input
                if (dmps_sp_min>0) THEN
                    ! Print user info
                    ! if (gtime%printnow) print FMT_SUB, 'Lowband upper limit: '//i2chr(dmps_sp_min)
                    ! Sets the concentration and composition to the partially added particles PSD based on the mfractions.
                    CALL send_conc(current_PSD%dp_range(1),nominal_dp(dmps_sp_min),conc_fit)
                    do i = 1, dmps_sp_min
                        if (GTIME%sec<1d0) THEN
                            CALL set_composition(i,nominal_dp(i), .false.)
                        else
                            CALL set_composition(i,nominal_dp(i), Use_old_composition)
                        end if
                    end do
                END IF

                ! If large particles are replaced with input
                if (dmps_sp_max>0) THEN
                    ! Print user info
                    ! if (gtime%printnow) print FMT_SUB, 'Highband lower limit: '//i2chr(dmps_sp_max)
                    ! Sets the concentration and composition to the partially added particles PSD based on the mfractions.
                    CALL send_conc(nominal_dp(dmps_sp_max),current_PSD%dp_range(2),conc_fit)
                    do i = dmps_sp_max, n_bins_par
                        if (GTIME%sec<1d0) THEN
                            CALL set_composition(i,nominal_dp(i), .false.)
                        else
                            CALL set_composition(i,nominal_dp(i), Use_old_composition)
                        end if
                    end do

                END IF

            END IF

            ! Next time next line from measurement etc. is read
            dmps_ln = dmps_ln + 1
        END IF PARTICLE_INIT

    add_particles: if (PRC%in_turn(4).and.J_TOTAL_M3>0d0) THEN
        ! ..........................................................................................................
        ! ADD NUCLEATED PARTICLES TO PSD, IN THE first couple of bins. This using the mixing method, where two particle
        ! distros are mixed together. NOTE that this step is always done if Aerosol_flag is on, even if ACDC_flag is off.
        ! In case no NPF is wanted, turn off ACDC and put NUC_RATE_IN to 0

        dmass = 0d0
        dconc_dep_mix = 0d0

        ! Negative mixing ratio makes this nucleation
        mix_ratio = -1d0

        ! New particles are assigned GENERIC composition.
        dmass(1:size(J_distr,1),VAPOUR_PROP%ind_GENERIC) = nominal_dp(1:size(J_distr,1))**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)
        ! The new particles are distributed over bins spanning bins that are 15% larger in diameter than the first bin, typically 3 bins
        dconc_dep_mix(1:size(J_distr,1)) = J_TOTAL_M3*GTIME%dt*minval(pack(speed_up, speed_up > 0))*J_distr

        CALL Mass_Number_Change('mixing')

        ! Update current_psd
        current_PSD = new_PSD

    END IF add_particles


        in_turn_cch_2: if (PRC%in_turn(PRC%cch)) THEN
        ! ..........................................................................................................
        ! CONDENSATION
        onlyIfCondIsUsed: if (Condensation .and.(.not. PRC%err)) THEN

            ! Pick the condensibles from chemistry and change units from #/cm^3 to #/m^3
            conc_vapour = 0d0
            conc_vapour(1:VAPOUR_PROP%n_cond_org-1) =  CH_GAS(index_cond)*1D6 ! mol/m3

            ! Poor sulfuric acid always needs special treatment
            conc_vapour(VAPOUR_PROP%ind_H2SO4) = H2SO4*1d6

            ! Update vapour pressures for organics
            VAPOUR_PROP%c_sat(1:VAPOUR_PROP%n_cond_org) =  saturation_conc_m3( &
                                VAPOUR_PROP%psat_a(1:VAPOUR_PROP%n_cond_org),  &
                                VAPOUR_PROP%psat_b(1:VAPOUR_PROP%n_cond_org),  &
                                GTEMPK) * VP_MULTI


            CALL UPDATE_MOLECULAR_DIFF_AND_CSPEED(VAPOUR_PROP)

            IF (CHEM_DEPOSITION) CALL CALCULATE_CHEMICAL_WALL_LOSS(conc_vapour(1:VAPOUR_PROP%n_cond_org),c_org_wall)

            dmass  = 0d0
            d_vap  = 0
            d_dpar = 0

            if (Kelvin_exp .and. abs(GTEMPK-T_prev)>5d-2) THEN
                Kelvin_Eff = pre_Kelvin ** (300d0/GTEMPK)
                T_prev = GTEMPK
            END IF

            if (Kelvin_taylor) Kelvin_Eff = 1d0 + pre_Kelvin * 300d0 / GTEMPK + (pre_Kelvin * 300d0 / GTEMPK) **2 / 2d0 + (pre_Kelvin * 300d0 / GTEMPK) **3 / 6d0

            CALL Condensation_apc(VAPOUR_PROP,conc_vapour,dmass, GTIME%dt*speed_up(PRC%cch),d_dpar,d_vap, kelvin_eff)

            ! calculate reference d_vap value, largest in magnitude but with sign
            call avg3(d_vap,refdvap, irefdvap)
            logdvap = abs(refdvap)>DVAPO_RANGE(2)

            ! ERROR HANDLING
            IF ( OPTIMIZE_DT.and.(maxval(ABS(d_dpar(2:))) > DDIAM_RANGE(2).or.d_dpar(1)>DDIAM_RANGE(2).or.logdvap )  &   ! ... if the changes in diam or vap. conc are too big
                .and. speed_up(PRC%cch)>1) THEN                                             ! ... only if OPTIMIZE_DT and we are not yet at minimum timestep
                IF (.not.PRC%err) THEN ! We actually might have error state from Condensation_apc, we don't want to overwrite it here
                    IF (maxval(ABS(d_dpar)) > DDIAM_RANGE(2)) THEN
                        call SET_ERROR(PRC%cch, 'Too large diameter change: '//f2chr(1d2*d_dpar(maxloc(abs(d_dpar),1)))//'%, bin# '//i2chr(maxloc(abs(d_dpar),1)))
                    ELSE IF ( logdvap ) THEN
                        call SET_ERROR(PRC%cch,'Too large d_vap: '//f2chr(refdvap*1d2)//'% in: '//TRIM(VAPOUR_PROP%vapour_names(irefdvap)))
                    ELSE
                        call SET_ERROR(PRC%cch,'Too large diameter and vapour concentration change: '//f2chr(1d2*maxval(abs(d_dpar)))//'%, '//f2chr(MAXVAL(abs(d_vap))))
                    END IF
                END IF
                dmass = 0.d0

            ELSE IF (.not.PRC%err) THEN ! everything fine -> apply changes

                CONDENSED_MASS = CONDENSED_MASS + SUM(dmass * &
                                    (TRANSPOSE(SPREAD(get_conc(current_PSD),1,VAPOUR_PROP%n_cond_tot))), 1)

                ! Distribute mass
                CALL Mass_Number_Change('condensation')

                ! Update current_psd
                current_PSD = new_PSD

                ! NOTE Update vapour concentrations to chemistry is done at the end of the timestep

                ! Check whether timestep can be increased:
                IF (OPTIMIZE_DT &
                    .and. maxval(ABS(d_dpar)) < DDIAM_RANGE(1) &
                    .and. abs(refdvap) < DVAPO_RANGE(1) &
                    ) THEN

                    if (speed_up(PRC%cch) * 2 * Gtime%dt <= DT_UPPER_LIMIT(PRC%cch)) THEN
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

            dconc_coag = 0.d0

            if (abs(GTEMPK-T_prev)>5d-2 .or. abs(GPRES-P_prev)>10_dp) THEN
                update_K = .true.
                T_prev   = GTEMPK
                P_prev   = GPRES
            END IF
            ! call opti() ! # 10

            ! Solve particle coagulation
            Call COAGULATION_ROUTINE(dconc_coag, GTIME%dt*speed_up(PRC%coa),d_npar, update_K)
            ! call opti() ! # 10

            ! ERROR HANDLING
            IF ( OPTIMIZE_DT .and. speed_up(PRC%coa)>1 .and. (MINVAL(d_npar) < 0 .or. maxval(ABS(d_npar)) > DPNUM_RANGE(2))) THEN   ! if the changes in particle numbers are too big
                call set_error(PRC%coa, 'Too large particle number change: '//f2chr(1d2*maxval(abs(d_npar)))//'%')
                dconc_coag = 0.d0
            END IF
            ! Distribute mass
            if (.not. PRC%err) Call Mass_Number_Change('coagulation')

            if (.not. PRC%err) THEN

                ! Update PSD with new concentrations
                current_PSD = new_PSD
                ! Check whether timestep can be increased:
                IF (maxval(ABS(d_npar)) < DPNUM_RANGE(1) .and. OPTIMIZE_DT) THEN
                    if (speed_up(PRC%coa) * 2 * Gtime%dt <= DT_UPPER_LIMIT(PRC%coa)) THEN
                        PRC%increase(PRC%coa) = .true.
                    END IF
                END IF
            END IF
        end if onlyIfCoagIsUsed
        ! end of coagulation
    END IF in_turncoa

    in_turn_dep: if (PRC%in_turn(PRC%dep)) THEN
        ! ..........................................................................................................
        ! Deposition
        OnlyIfDepoIsUsed: if ((Deposition) .and.(.not. PRC%err)) then
            ! Solve particle coagulation
            dconc_dep_mix = 0.d0

            IF (LOSSES_FILE /= '') THEN
                if (CONSTANT_PAR_LOSS_RATE>=0d0) THEN
                    losses_fit=CONSTANT_PAR_LOSS_RATE
                ELSE
                    ! The losses file is interpolated spatially and temporally to fit the current bin structure and time
                    losses_fit = [ ((INTERP( PAR_LOSSES%sections,                                                             &
                                 [((INTERP(PAR_LOSSES%time, PAR_LOSSES%conc_matrix(:,i),unit=LOSSFILE_TIME_UNIT)), i=1,size(PAR_LOSSES%sections, 1))], &
                                  timein=nominal_dp(j)) ), j=1,n_bins_par)]
                end if

                ! Deposited concentratios calculated here
                dconc_dep_mix = get_conc() * (1 - EXP(-losses_fit*GTIME%dt*speed_up(PRC%dep)))

            ELSE
                call deposition_velocity(get_dp(),ustar,A_vert,CHAMBER_FLOOR_AREA,CHAMBER_FLOOR_AREA,&
                                        Vol_chamber,GTEMPK,GPRES,E_field,dconc_dep_mix,GTIME%dt*speed_up(PRC%dep),losses_fit)
            END IF

            ! ERROR HANDLING; Check whether changes are within limits:
            d_npdep = 0.d0
            WHERE (get_conc()>0d0) d_npdep = dconc_dep_mix/get_conc()
            IF (maxval(ABS(d_npdep)) > DPNUM_RANGE(2) .and. OPTIMIZE_DT .and. speed_up(PRC%dep)>1) THEN   ! if the changes in diameter are too big
                call SET_ERROR(PRC%dep,'Too large particle number change'//f2chr(1d2*maxval(abs(d_npdep)))//'%')
                dconc_dep_mix = 0.d0

            ELSE  ! everything fine -> apply changes
                ! Distribute mass
                Call Mass_Number_Change('deposition')

                ! Save lost mass in Depos_composition
                Depos_composition = Depos_composition + SUM(get_composition(current_PSD)* &
                                    (TRANSPOSE(SPREAD((get_conc(current_PSD) - get_conc(new_PSD)),1,VAPOUR_PROP%n_cond_tot))), 1)

                ! Update PSD with new concentrations
                current_PSD = new_PSD
                !Check whether timestep can be increased:
                IF (maxval(ABS(d_npdep)) < DPNUM_RANGE(1) .and. OPTIMIZE_DT) THEN
                    if (speed_up(PRC%dep) * 2 * Gtime%dt <= DT_UPPER_LIMIT(PRC%dep)) THEN
                        PRC%increase(PRC%dep) = .true.
                    END IF
                END IF

            END IF

            ! WHY WAS THIS HERE??! Update PSD with new concentrations
            ! current_PSD = new_PSD
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

        CALL error_handling

        ! Reset the system to the state at previous timestep
        current_PSD = old_PSD
        CH_GAS      = CH_GAS_old
        CH_RO2      = CH_RO2_old
        acdc_goback = .true.
        dmps_ln     = dmps_ln_old
        if (Chem_Deposition) c_org_wall_old = c_org_wall

        ! Reset relative change vectors documenting the precision:
        d_dpar      = 0.d0
        d_npar      = 0.d0
        d_vap       = 0.d0
        d_npdep     = 0.d0

        ! If there was no error during the actual timestep
    ELSE

        ! Write printouts to screen
        if (GTIME%printnow) THEN
            CALL PRINT_HEADER
            ! Calculate growth rates
            IF (CALC_GR) CALL PRINT_GROWTH_RATE
            CALL PRINT_KEY_INFORMATION(TSTEP_CONC)

        ELSE
            ! if --gui flag was used, print a dot and EOL so the STDOUT-reading in Python GUI will be smoother.
            if (ingui) print'(a)', '.'
            call cpu_time(cpu3)
            if (cpu3-cpu4>15d0) THEN
                print*, '... ',GTIME%hms, '......... '
                cpu4=cpu3
            END IF

        END IF

        ! Save to netcdf
        if (GTIME%savenow) THEN

            ! PSD is also saved to txt
            if (Aerosol_flag) THEN
                WRITE(601,*) GTIME%sec, sum(get_conc()*1d-6), get_conc()*1d-6 / LOG10(bin_ratio)
                WRITE(604,*) GTIME%sec, get_conc()*1d-6
                ! WRITE(605,*) GTIME%sec, G_COAG_SINK
                save_measured = conc_fit/1d6
            END IF

            WRITE(610,*) GTIME%sec, d_dpar(maxloc(abs(d_dpar),1)), d_npar(maxloc(abs(d_npar),1)), refdvap, d_npdep(maxloc(abs(d_npdep),1))
            FLUSH(610)

            CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,reactivities,conc_vapour*1d-6,VAPOUR_PROP, save_measured,&
                            1d9*3600/(GTIME%dt*speed_up(PRC%cch))*get_dp()*d_dpar,losses_fit)
            ! At this point the simulation can be stopped if desired
            if (ENABLE_END_FROM_OUTSIDE) CALL CHECK_IF_END_CMD_GIVEN

        END IF

        ! Add main timestep to GTIME
        if (GTIME%sec + GTIME%dt <= GTIME%SIM_TIME_S) THEN

            add_rounds = minval(pack(speed_up, speed_up > 0))! - modulo(n_of_Rounds,add_rounds)

            GTIME = ADD(GTIME, GTIME%dt*add_rounds)

            n_of_Rounds = n_of_Rounds + add_rounds

            if (PRC%in_turn(1)) bookkeeping(1,1) = bookkeeping(1,1) + 1_dint
            if (PRC%in_turn(2)) bookkeeping(2,1) = bookkeeping(2,1) + 1_dint
            if (PRC%in_turn(3)) bookkeeping(3,1) = bookkeeping(3,1) + 1_dint
            will_stall = 0

            ! Update vapour concentrations to chemistry
            if (Condensation) &
            CH_GAS(index_cond) = conc_vapour(1:VAPOUR_PROP%n_cond_org-1) *1D-6
            ! ----------------------------------------

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

! print*, 'total rounds: ',n_of_Rounds
! print*, 'kelvin calculated: ', i_kel

if (Aerosol_flag) THEN
    CLOSE(601)
    CLOSE(604)
    CLOSE(610)
END IF


! Close output file netcdf
CALL CLOSE_FILES(RUN_OUTPUT_DIR)
if (OPTIMIZE_DT) THEN
    print*, 'Total rounds for chemistry and/or condensation: ', bookkeeping(1,1)
    print*, 'Total rounds for nucleation and/or coagulation: ', bookkeeping(2,1)
    print*, 'Total rounds for deposition                   : ', bookkeeping(3,1)
END IF

CALL FINISH


CONTAINS


! Print a pretty header
SUBROUTINE PRINT_HEADER
    IMPLICIT NONE
    WRITE(*,*)
    ! Print time
    if (GTIME%dt*minval(pack(speed_up, speed_up>0))>1d0) THEN
        print FMT_TIME, GTIME%hms//' ('//TRIM(i2chr(int(GTIME%sec)))//' s)'
    ELSE
        print FMT_TIME, GTIME%hms//' ('//TRIM(i2chr(int(GTIME%sec*1000d0)))//' ms)'
    END IF
    ! To compare real time vs simulation time, timer is stopped in the beginning
    call cpu_time(cpu2)
    cpu4=cpu2
END SUBROUTINE PRINT_HEADER


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
    ! This header is not printed/saved if it was already showed in this times step
    if (bookkeeping(PRC%proc,2) < n_of_Rounds) THEN
        CALL OUTPUTBOTH(608,FMT_NOTE0,'Precision tolerance exceeded in '//TRIM(UCASE(PRC%pr_name(PRC%proc)))//&
                                      ' at time '//GTIME%hms//' ('//TRIM(f2chr((GTIME%sec)))//' sec)',-1)
    END IF
    CALL OUTPUTBOTH(608,FMT_SUB,  ' What happened: '//trim(PRC%err_text),-1)
    ! reduce the speed_up factor or, if necessary, the integration time step:
    IF (speed_up(PRC%proc) > 2) THEN

        speed_up(PRC%proc) = speed_up(PRC%proc) / 2

        CALL OUTPUTBOTH(608,'(" - ",a)', GTIME%hms//' --> Reducing '//TRIM(PRC%pr_name(PRC%proc))//' to '&
                                        //di2chr(speed_up(PRC%proc))//' (new process dt: '//TRIM(f2chr(GTIME%dt*speed_up(PRC%proc)))//' s)')
        errortime(PRC%proc) = GTIME%sec

    ELSE IF (speed_up(PRC%proc) == 2) THEN
            print*,'minimum dt would be reached in ', TRIM(PRC%pr_name(minloc(abs(speed_up),1))), '-> halving all other timesteps'
            WHERE (speed_up > 2) speed_up = speed_up / 2
            if (maxval(speed_up) == 2) THEN
                WHERE (speed_up == 2) speed_up = 1
            END IF
            print*, 'New multiplyers:', speed_up
    ELSE
        will_stall = will_stall+1
        if (PRC%proc==PRC%cch.and.Chemistry_flag.and..not.Condensation) THEN
            CALL OUTPUTBOTH(608,FMT_WARN0, ' => Reduce time step or increase dC tolerance, chemistry can usually handle it ')
        else
            CALL OUTPUTBOTH(608,FMT_WARN0, ' => YOU SHOULD REDUCE MAIN TIMESTEP ')
        END IF
        if (will_stall>1000) THEN
            print *,''
            print *, 'Unfortunately it seems the simulation is getting nowhere, so stopping now. Your options are:'
            print *, '  1) Increase tolerances'
            print *, '  2) Reduce minimum time step'
            print *, '  3) Use constant timestep and accept the outcome'
            print *, '  4) If nothing helps and you think you have reasonable input, contact support'
            print *, ''
            CALL CLOSE_FILES(RUN_OUTPUT_DIR)
            CALL FINISH
            STOP 'Saved files, bye.'
        END IF
    END IF

    ! CALL OUTPUTBOTH(608, FMT_MSG,'Timesteps adjusted.')
    ! print FMT_LEND,

    bookkeeping(PRC%proc,2) = n_of_Rounds
    CALL OUTPUTBOTH(608,'*','',-1)
    flush(608)

    ! RESET ERROR STATE
    PRC%err = .false.
    PRC%proc = 0
    PRC%err_text = ''

    ! if (minval(abs(speed_up))==1) THEN
    !     print*,'minimum dt reached in ', TRIM(PRC%pr_name(minloc(abs(speed_up),1))), '-> halving all timesteps'
    !     WHERE (speed_up > 1) speed_up = speed_up / 2
    !     print*, 'New multiplyers:', speed_up
    ! END IF

END SUBROUTINE error_handling


! Routine to increase speed so that unused speed factors stay as
! the highest number and are not hindering main timestep changes
Subroutine INCREASE_SPEED(ii)
    IMPLICIT NONE
    INTEGER :: ii

    if (PRC%increase(ii) .and. GTIME%sec-errortime(ii)>(2.1d0*GTIME%dt * speed_up(ii))) THEN
        speed_up(ii) = speed_up(ii) * 2
        CALL OUTPUTBOTH(608,'*','+ '//GTIME%hms//' --> Speeding '//TRIM(PRC%pr_name(ii))//' to '&
                        //TRIM(di2chr(speed_up(ii)))//' (new process dt: '//TRIM(f2chr(GTIME%dt*speed_up(ii)))//' s)' )
    END IF
    PRC%increase(ii) = .false.
    errortime(ii) = -999d0
    if (.not. Deposition ) speed_up(PRC%dep) = -9999_dint
    if (.not. Coagulation) speed_up(PRC%coa) = -9999_dint
    if (.not. Condensation .and..not. Chemistry_flag) speed_up(PRC%cch) = -9999_dint

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
! HOW TO ADD MORE ACDC SUBMODULES?
! - Duplicate any of the src/ACDC/ACDC_* directories and rename it with the next available number.
!   If, for example one adds a sixth module, the new directory is named ACDC_06
! - In the new directory ACDC_06, delete files acdc_equations_0x*.f90 and acdc_system_0x*.f90
! - In the new directory ACDC_06, name all the remaining .f90 files with *_0x#.f90 -> *_0x6.f90
! - In the new system directory, find and replace 0X* with 0X6 (case insensitive search/replace) INSIDE
!   the following files:
!      - acdc_simulation_setup_0x*.f90
!      - driver_acdc_J_0x*.f90 and
!      - get_acdc_J_0x*.f90
! - In this subroutine, update the USE statements and SELECT_PROPER_ACDC_SUBMODULE using the same logic
! - In the GUI->Cluster formation, update the new system settings and press Update ACDC_06 and recompile ARCA
! - The makefile should not need any updating, but in case compiling error, run make from the terminal:
! make cleanmain && make
! =================================================================================================


SUBROUTINE ACDC_J(C, dt)
    use get_acdc_J_0X1, only : acdc_plugin_0X1=>acdc_plugin
    use get_acdc_J_0X2, only : acdc_plugin_0X2=>acdc_plugin
    use get_acdc_J_0X3, only : acdc_plugin_0X3=>acdc_plugin
    use get_acdc_J_0X4, only : acdc_plugin_0X4=>acdc_plugin
    use get_acdc_J_0X5, only : acdc_plugin_0X5=>acdc_plugin

    implicit none

    REAL(dp), intent(in) :: C(:)
    REAL(dp), intent(in) :: dt
    LOGICAL, save        :: first_time = .true., ss_handle
    REAL(dp)             :: IPR = 0 ! Ion pairs / second / cm3, converted to Ion pairs / second / m3
    INTEGER :: ii,jj

    if (first_time) THEN
        ss_handle = ACDC_solve_ss
        ! if (RESOLVE_BASE) ss_handle = .TRUE.
        ! Now set to false in the end of the routine !
        ! first_time = .false.
    end if

    IF (inm_IPR   /= 0) IPR   = C(inm_IPR)*1d6

    ! -------- NEW ACDC -------------
    do ii=1,size(G_ACDC)

        if (G_ACDC(ii)%inuse) THEN
            do jj=1,size(G_ACDC(ii)%ACDC_MONOMER_NAMES)
                ! by default, G_ACDC(ii)%ACDC_monConc(jj) = 0d0
                If (G_ACDC(ii)%ACDC_LINK_IND(jj)>0) THEN
                    if (G_ACDC(ii)%ACDC_LINK_IND(jj) == inm_H2SO4) THEN
                        G_ACDC(ii)%ACDC_monConc(jj) = H2SO4*1d6
                    ELSE
                        G_ACDC(ii)%ACDC_monConc(jj) = C(G_ACDC(ii)%ACDC_LINK_IND(jj))*1d6
                    END IF
                ELSE if (G_ACDC(ii)%ACDC_LINK_IND(jj)<0) THEN
                    G_ACDC(ii)%ACDC_monConc(jj) = CH_GAS(-1*G_ACDC(ii)%ACDC_LINK_IND(jj))*1d6
                END IF
            end do
            SELECT_PROPER_ACDC_SUBMODULE: select case (ii)
                case (1)
                    CALL acdc_plugin_0X1(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,GCS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, acdc_goback, (GTIME%printnow.and.PRINT_ACDC))
                case (2)
                    CALL acdc_plugin_0X2(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,GCS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, acdc_goback, (GTIME%printnow.and.PRINT_ACDC))
                case (3)
                    CALL acdc_plugin_0X3(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,GCS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, acdc_goback, (GTIME%printnow.and.PRINT_ACDC))
                case (4)
                    CALL acdc_plugin_0X4(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,GCS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, acdc_goback, (GTIME%printnow.and.PRINT_ACDC))
                case (5)
                    CALL acdc_plugin_0X5(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,GCS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, acdc_goback, (GTIME%printnow.and.PRINT_ACDC))
            end select SELECT_PROPER_ACDC_SUBMODULE

            G_ACDC(ii)%J_OUT_CM3(:) = 1d-6 * G_ACDC(ii)%J_OUT_M3(:)
            ! Missing values are substituted w -1
            WHERE (G_ACDC(ii)%J_OUT_CM3(:)<-9.999d-7 ) G_ACDC(ii)%J_OUT_CM3(:) = -1d0

            if (first_time) print FMT_SUB, 'ACDC system '//TRIM(i2chr(ii))//', Outgrowing cluster size: '//TRIM(f2chr(G_ACDC(ii)%Cl_diam))
            if (first_time.and.Aerosol_flag) THEN
                if (nominal_dp(1)/G_ACDC(ii)%Cl_diam > 1d0) THEN
                    print FMT_FAT0, 'System nr: '//TRIM(i2chr(ii))
                    print *, 'Cluster diameter ', G_ACDC(ii)%Cl_diam
                    STOP 'Smallest particle bin is larger than outgrowing cluster'
                else if (nominal_dp(1)/G_ACDC(ii)%Cl_diam < 0.9d0) THEN
                    print FMT_WARN0, 'System nr: '//TRIM(i2chr(ii))//', outgrowing cluster diameter:'//TRIM(f2chr(G_ACDC(ii)%Cl_diam))
                    print FMT_WARN0, 'Outgrowing cluster diameter is '//TRIM(f2chr(G_ACDC(ii)%Cl_diam/nominal_dp(1)))//' times larger than the bin where it is put'
                    print FMT_WARN0, 'solution: increase minimum Dp'
                end if
            end if
        END IF
    END DO

    if (first_time) first_time = .false.
    acdc_goback = .false.



END SUBROUTINE ACDC_J  ! END ACDC Nucleation


! =================================================================================================
! Parametrisation of organic nucleation as used in Roldin et al., Nat Comm. 2019
! .................................................................................................
SUBROUTINE ORGANIC_NUCL(J_TOTAL_M3)
    real(dp) :: J_TOTAL_M3
    real(dp) :: ORGS,kJ,dH,J15
    real(dp), parameter :: dG = -15.1d3   ! [cal/mol]
    real(dp), parameter :: dS = 61.1      ! [cal/mol/K] NOTE this value had a typo in Nat. Comm. and should be positive!
    real(dp), PARAMETER :: Rcal = 1.98720 ! gas constant in calories/mol/K
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
            print FMT_FAT0, "List of nucleating organic compounds has vapours which are not in the chemistry."
            print FMT_SUB, "Options: Edit the file 'ModelLib/required/nucl_homs.txt' or turn of organic nucleation."
            print FMT_LEND,
            STOP "  ----------- Bye -------------"
        END IF
        first_run = .False.
    END IF
    ORGS = sum(CH_GAS(inds))

    dH = dG - dS*GTEMPK
    kJ = 5d-13*EXP(-dH/(Rcal) * (1/GTEMPK - 1/298d0))

    J15 = kJ * ORGS*H2SO4
    if (GTIME%printnow) print FMT_MSG, 'Organic formation rate: '//TRIM(ADJUSTL(f2chr(J15)))//' [/s/cm3]'
    J_TOTAL_M3 = J_TOTAL_M3 + J15*1d6

END SUBROUTINE ORGANIC_NUCL


! ================================================================================================
! Chemical wall loss - only for condensible vapours, based on model for Teflon coated chamber/bag
! k_g2w and k_w2g based on McMurry PH & Grosjean D (1985) Gas and aerosol wall losses in
! Teflon film smog chambers. Environ. Sci. Technol. 19(12):1176-1182. See also Zang et al., 2014
! EDDYK = Coefficient of eddy diffusion (s^-1), usually treated as an unknown fitting parameter
! ================================================================================================
SUBROUTINE CALCULATE_CHEMICAL_WALL_LOSS(conc, c_org_wall)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: conc(:),c_org_wall(:)
    INTEGER,PARAMETER  :: NEQ = 2
    REAL(kind(1.d0))   :: Y(NEQ), RWORK(52),RTOL=1D-4,ATOL(NEQ)=[1.D-8,1.D-14],RPAR(NEQ),T, TOUT
    INTEGER            :: IWORK(32),IOUT,LRW=52,LIW=32,MF=21,IOPT=0,ISTATE,ITASK=1,ITOL=2
    real(dp)           :: old_gas_mass
    real(dp)           :: k_g2w(VAPOUR_PROP%n_cond_org) ! Reversible wall loss rate gas->wall
    real(dp)           :: k_w2g(VAPOUR_PROP%n_cond_org) ! Reversible wall loss rate wall->gas
    integer            :: nc,ii ! nc=short for number of condensible organic vapours, loop indices

    nc = VAPOUR_PROP%n_cond_org

    k_g2w = ((A_vert + 2d0*Chamber_floor_area)/Vol_chamber)*(VAPOUR_PROP%alphawall(1:nc)*VAPOUR_PROP%c_speed(1:nc)/4D0)/(1D0+(pi/2D0) &
                 * (VAPOUR_PROP%alphawall(1:nc)*VAPOUR_PROP%c_speed(1:nc)/(4D0*sqrt(EDDYK*VAPOUR_PROP%diff(1:nc)))))
    k_w2g = k_g2w/(Na/VAPOUR_PROP%c_sat(1:nc) * Cw_eqv) ! Cw_eqv=40 mumol/m^3, => ~10 mg/m^3 as in Zhang et al., 2014

    old_gas_mass = SUM(VAPOUR_PROP%molec_mass(1:nc)*conc)

    do ii = 1,VAPOUR_PROP%n_cond_org
        Y           = [conc(ii),c_org_wall(ii)]
        RPAR        = [k_g2w(ii),k_w2g(ii)]
        T           = 0d0
        TOUT        = GTIME%dt*speed_up(PRC%cch)
        ISTATE      = 1
        CALL DVODE(WLOSSEQN,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,WLJAC,MF,RPAR)
        VAP_DEP_MASS_WALLS(ii) = VAP_DEP_MASS_WALLS(ii) + VAPOUR_PROP%molec_mass(ii) * (Y(2) - c_org_wall(ii))
        conc(ii)    = Y(1)
        c_org_wall(ii) = Y(2)
    end do

    if (GTIME%printnow .and. SUM(conc)>0d0) THEN
        print'(" * ",es10.2,"% of condensible vapour mass lost to walls in one timestep")', &
            1d2*(1d0-SUM(VAPOUR_PROP%molec_mass(1:nc)*conc)/old_gas_mass)
        print '(a,es12.3)', '   Organic vapours on walls [µg]  :', sum(VAPOUR_PROP%molec_mass(1:nc)*c_org_wall)*1d9
        print '(a,es12.3)', '   Organic vapours in chamber [µg]:', Vol_chamber * sum(VAPOUR_PROP%molec_mass(1:nc)*conc)*1d9
    END IF

END SUBROUTINE CALCULATE_CHEMICAL_WALL_LOSS

SUBROUTINE WLOSSEQN (NEQ, T, Y, YDOT, RPAR)
    INTEGER :: NEQ
    REAL(kind(1.d0)):: RPAR(NEQ), T,Y(NEQ), YDOT(NEQ)
    YDOT(1) = -RPAR(1)*Y(1) + RPAR(2)*Y(2)
    YDOT(2) = -YDOT(1)
end subroutine WLOSSEQN

SUBROUTINE WLJAC (NEQ, T, Y, ML, MU, PD, NRPD, RPAR)
    INTEGER :: NEQ, ML, MU, NRPD
    REAL(kind(1.d0)):: RPAR(NEQ), T, PD(NRPD,NEQ),Y(NEQ), YDOT(NEQ)
    PD(1,1) = -RPAR(1)
    PD(1,2) =  RPAR(2)
    PD(2,1) = -PD(1,1)
    PD(2,2) = -PD(1,2)
END SUBROUTINE WLJAC

! ================================================================================================
! PRINT_KEY_INFORMATION for the user. Can be called for example upon GTIME%printnow with
! if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
! ================================================================================================
SUBROUTINE PRINT_KEY_INFORMATION(C)
    IMPLICIT NONE
    REAL(dp), intent(in) :: C(:)

    if (Condensation) print FMT_MSG, 'Max change in vapours '&
            //TRIM(f2chr(1d2* d_vap(irefdvap) ))//'% for '//VAPOUR_PROP%vapour_names(irefdvap)
    if (.not.Condensation.and.Chemistry_flag) print FMT_MSG, 'Max change in chemistry '&
            //TRIM(f2chr(1d2* d_chem(maxloc(abs(d_chem),1)) ))//'% for '//SPC_NAMES(maxloc(abs(d_chem),1))
    if (Aerosol_flag) print FMT_MSG, 'Max change in par conc. '&
            //TRIM(f2chr(1d2*d_npar(maxloc(d_npar,1))))//'% in bin # '//i2chr((maxloc(d_npar,1)))
    if (Aerosol_flag) print FMT_MSG, 'Max change in par diam. '&
            //TRIM(f2chr(1d2*d_dpar(maxloc(d_dpar,1))))//'% in bin # '//i2chr((maxloc(d_dpar,1)))

    print FMT10_3CVU,'Temp:', C(inm_TempK), ' [K]','Pressure: ', C(inm_pres), ' [Pa]', 'Air_conc', C_AIR_cc(C(inm_TempK), C(inm_pres)), ' [1/cm3]'
    if (OH_ind_in_chemistry>0) THEN
        print FMT10_2CVU, 'H2SO4: ', H2SO4, ' [1/cm3]','OH: ',CH_GAS(ind_OH),' [1/cm3]'
    else
        print FMT10_CVU, 'H2SO4: ', H2SO4, ' [1/cm3]'
    end if
    IF (inm_NH3   /= 0) print FMT10_2CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', G_ACDC(1)%J_OUT_M3(1)*1d-6, ' [1/s/cm3]'!,'sum(Nn)/N1',clusterbase,' []'
    IF (inm_DMA   /= 0) print FMT10_3CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', G_ACDC(2)%J_OUT_M3(1)*1d-6, ' [1/s/cm3]', 'J-total', J_TOTAL_M3*1e-6,' [1/s/cm3]'
    print FMT10_3CVU, 'Jion neut:', sum(G_ACDC(1:2)%J_OUT_CM3(2)) , ' [/s/cm3]','Jion pos:', SUM(G_ACDC(1:2)%J_OUT_CM3(3)) , ' [/s/cm3]','Jion neg:', SUM(G_ACDC(1:2)%J_OUT_CM3(4)) , ' [/s/cm3]'
    print FMT10_2CVU, 'C-sink:', GCS , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
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
            save_measured = conc_fit/1d6
        END IF

        CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,reactivities,conc_vapour*1d-6, VAPOUR_PROP, save_measured,&
                        1d9*3600/(GTIME%dt*speed_up(PRC%cch))*get_dp()*d_dpar,losses_fit)

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

    if (model_H2SO4) THEN
        if (H2SO4_ind_in_chemistry<1) THEN
            print FMT_FAT0, 'H2SO4 should explicitly come from chemistry but it does not exist there. Either '
            print FMT_SUB, '  uncheck "Replace any input H2SO4 with modelled..." and include H2SO4 to input, or '
            print FMT_SUB, '  use chemistry that has H2SO4 included.'
            stop ' '
        ELSE
            MODS(inm_H2SO4)%ISPROVIDED = .false.
            print FMT_SUB, 'Replacing HSO4 input with modelled chemistry'
        END IF
    END IF

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
                    exit
                end if
            END DO
            IF (Chemistry_flag.and.check == 0 .and. I>LENV) THEN
                print FMT_FAT0, 'You are using an (organic?) compound which does not exist in chemistry: '//TRIM(MODS(i)%NAME)//' '
                IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                IF (ABS(MODS(I)%SHIFT) > 0d0) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
                print FMT_MSG, 'Good bye.'
                STOP
            END IF
            IF (Chemistry_flag.and.check == 0 .and. I<LENV .and. &
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
! Calculate water saturation vapour pressure (ES), current vapour pressure (EW, in Pa), concentration in
! molec/cm³ (CW) and surface tension STW [n/m²]
! ============================================================================================================
SUBROUTINE WATER(ES,EW,CW,STW)
    IMPLICIT NONE
    REAL(dp), INTENT(INOUT) :: ES, EW, CW, STW
    REAL(dp)                :: TEMPC
    ! Saturation vapour pressure of water in Pa
    REAL, PARAMETER        :: a0 = 6.107799961,     & ! Parameters to calculate the saturation vapour pressure for water
                              a1 = 4.436518524E-1,  &
                              a2 = 1.428945805E-2,  &
                              a3 = 2.650648471E-4,  &
                              a4 = 3.031240396E-6,  &
                              a5 = 2.034080948E-8,  &
                              a6 = 6.136820929E-11


    TEMPC = GTEMPK-273.15d0

    ! Saturation vapour pressure over liquid water; using parametrisation from ???, a0,a1...a6 are in constants.f90
    ES = (a0 + a1 * TEMPC**1 + a2 * TEMPC**2 + a3 * TEMPC**3    &
             + a4 * TEMPC**4 + a5 * TEMPC**5 + a6 * TEMPC**6)* 100

    ! Water vapour pressure
    EW  = GRH * ES / 100d0

    ! Water vapour concentration in molecules per cm3
    CW = EW/GPRES * GC_AIR_NOW

    ! Water surface tension N/m²
    STW =  0.117296d0 - 0.152362d-3 * GTEMPK

END SUBROUTINE WATER


! ================================================================================================
! Calculate Beta parameter (solar angle [°] above horizon), parameterization by Henderson-Sellers
! ================================================================================================
SUBROUTINE BETA(CH_Beta) ! [degrees]
    IMPLICIT NONE
    REAL(dp), INTENT(INOUT) :: CH_Beta
    REAL(dp)                :: D_zsd,HRANG,ZSR

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
    INTEGER      :: for, rof
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

! Fill vector INDRELAY_TIED with the corresponding source variable indices
SUBROUTINE LINK_VARIABLES
    IMPLICIT NONE
    INTEGER :: i
    print FMT_HDR, 'Checking for linked variables'
    DO i=1,N_VARS
        IF ((MODS(i)%ISPROVIDED).and.(MODS(i)%TIED /= '')) THEN
            if (I < 3) THEN
                stop "  ------- Tying temperature or pressure to other variables is not possible. ---"
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

! Handles output from subroutine Error_handling. Output 'text' to screen and/or 'unit' with 'format'
! Targets:
!    tar == -1: text file (or specifically, to 'unit')
!    tar ==  1: screen
!    tar ==  0: both, default if 'tar' is not provided
SUBROUTINE OUTPUTBOTH(unit, format, text,tar)
    IMPLICIT NONE
    INTEGER, INTENT(in)          :: unit
    CHARACTER(len=*), INTENT(in) :: format
    CHARACTER(len=*), INTENT(in) :: text
    INTEGER, OPTIONAL            :: tar
    INTEGER                      :: w

    if (PRESENT(tar)) then
        w = tar
    else
        w = 0
    end if
    if (FORMAT == '*') THEN
        if (w < 1)  write(unit,*) text
        if (w > -1) write(*,*) text
    ELSE
        if (w < 1)  write(unit,Format) text
        if (w > -1) write(*,format) text
    END IF

END SUBROUTINE OUTPUTBOTH

! This subroutine makes it possible to end a simulation from outside using a file called ENDNOW.INIT
! To save IO time, only the existence of the file is first checked. If file ENDNOW.INIT exists, then
! command STOP is checked. If it is found, the simulation is stopped by setting the maximum simulation
! time to current time.
SUBROUTINE CHECK_IF_END_CMD_GIVEN
    implicit none
    integer            :: ioi
    CHARACTER(len=256) :: command, date, time

    OPEN(666, file=RUN_OUTPUT_DIR//"/ENDNOW.INIT", iostat=ioi, STATUS='OLD',action='read')
    if (ioi == 0) THEN
        read(666,*,iostat=ioi) command
        close(666)
        if ((ioi == 0).and.(TRIM(command)=='STOP')) THEN
            print FMT_MSG, 'ENDNOW.INIT WAS USED, STOPPING THE SIMULATION FROM OUTSIDE. GOOD BYE.'
            GTIME%SIM_TIME_S = GTIME%sec
            call date_and_time(date,time)
            OPEN(666, file=RUN_OUTPUT_DIR//"/ENDNOW.INIT", iostat=ioi, STATUS='REPLACE',action='WRITE')
            write(666,*) 'Run was terminated by the user at simulation time '//GTIME%hms//'. Real date is: '//TRIM(date)//', time is '//TRIM(time(1:6))
            close(666)
            ioi = RENAME(RUN_OUTPUT_DIR//"/ENDNOW.INIT", RUN_OUTPUT_DIR//"/ENDNOW._OLD")
            if (ioi /= 0) print*, 'Error while renaming the end cmd file.'
        END IF
    END IF

END SUBROUTINE CHECK_IF_END_CMD_GIVEN

END PROGRAM ARCA_main
