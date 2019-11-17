PROGRAM Supermodel

    USE second_MAIN                         ! Main second file
    USE second_PARAMETERS                   ! CH_NSPEC (originally NSPEC) and chemical indices, ind_xxxx, come from here
    USE second_Precision,  ONLY : dp        ! KPP Numerical type
    ! USE SECOND_PARAMETERS, ONLY : NSPEC     ! Number of chemical species
    ! USE second_Global,     ONLY : NPHOT     ! Number of photochemical reactions used in KPP
    USE second_Monitor,    ONLY : SPC_NAMES ! Names of chemicals from KPP
    USE Chemistry
    USE constants
    USE AUXILLARIES
    USE INPUT
    use OUTPUT
    USE ACDC_NH3
    USE ACDC_DMA
    USE aerosol_auxillaries
    USE Aerosol


    IMPLICIT NONE

    INTEGER  :: I

    ! MOST OF THE VARIABLES ARE DEFINED IN INPUT.F90
    REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)
    REAL(DP), ALLOCATABLE :: CH_GAS(:)
    REAL(dp) :: J_ACDC_NH3 = 0d0
    REAL(dp) :: J_ACDC_DMA = 0d0
    REAL(dp) :: J_NH3_BY_IONS(3) = 0d0
    REAL(dp) :: acdc_cluster_diam = 2.17d-9

    REAL(dp) ::                                &
        CH_RO2,                                &   ! RO2 concentration in [molecules / cm^3]
        CH_H2O,                                &   ! H20 concentration in [molecules / cm^3]
        CH_TIME_kpp,                           &   ! Start time in KPP
        CH_END_kpp,                            &   ! End time in KPP
        CH_Beta,                               &   ! solar zenit angle
        EW,ES                                      ! Water content in Pa, Saturation vapour pressure

    ! Lukas (and others), these are leftovers from Michael that are defined otherwise but I didn't want to delete them yet (PC)
    ! CH_RES1, CH_RES2, CS_H2SO4, CS_HNO3,   &   ! CS for H2SO4 and HNO3 for KPP
    ! CH_AIR,                                &   ! AIR concentration
    ! CH_Glo,                                &   ! Global shortwave radiation
    ! dt_chem,                               &   ! Chem DT

    ! Carlton, are these yours, I just commented them since they are not in use now (PC)
    ! REAL(dp),DIMENSION(:), ALLOCATABLE :: values !vapours%vapour_number
    ! REAL(dp) :: timestep, end_time


    !speed_up: factor for increasing integration time step for individual prosesses
    !...(1): Photolysis, (2): chemistry; (3):Nucleation; (4):Condensation; (5): Coagulation; (6): Deposition
    INTEGER :: speed_up(10) = 1
    TYPE(error_type) :: error

    TYPE(aerosol_setup)       :: AER_setup
    TYPE(particle_properties) :: AER_par_prop
    TYPE(initial_distribution):: AER_init_dist
    TYPE(ambient_properties)  :: ambient

    CALL read_input_data ! Declare most variables and read user input and options in input.f90

    CALL CHECK_INPUT_AGAINST_KPP ! Check that the input exists in chemistry, or if not, print warning


    ! ???
    if (VAP_logical) then
        call Aerosol_intialization(AER_setup, AER_par_prop, vapours, AER_init_dist,ambient)
        !call condensation_routine(values,timestep,end_time,AER_setup, AER_par_prop,ambient)
        print*, shape(vapours%c_sat)
    end if

    ALLOCATE(TSTEP_CONC(N_VARS))
    TSTEP_CONC = 0
    ALLOCATE(CH_GAS(size(SPC_NAMES)))
    CH_GAS = 0

    ! This only called once for KPP in the beginning
    IF (Chemistry_flag .EQV. .true.) THEN
      CALL KPP_SetUp
    ENDIF

    !Open error output file
    open(unit=333, file='output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'_error_output.txt')

    CALL OPEN_FILES(('output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)), Description, MODS, CH_GAS, VAPOURS)

    !Open output file


    ! do i=1, size(XTRAS)
    !   print*,
    !   print*, XTRAS(I)%name
    !
    !   print'("           ", 38(es9.2))', XTRAS(I)%sections
    !   do k=1, size(XTRAS(I)%time)
    !     print'(39(es9.2))', XTRAS(I)%time(k),XTRAS(I)%binseries(K,:)
    !   end do
    ! end do
    !
    CALL PAUSE_FOR_WHILE(wait_for)

    print*,;print FMT_HDR, 'Beginning simulation' ! Information to user

    ! =================================================================================================
    DO WHILE (MODELTIME%SIM_TIME_S - MODELTIME%sec > -1d-12) ! MAIN LOOP STARTS HERE
        ! =================================================================================================
        ! Print time header in a nice way, start with empty row
        if (MODELTIME%printnow) print *, ! Print time
        if (MODELTIME%printnow) print FMT_TIME, MODELTIME%hms

        ! =================================================================================================
        ! C_AIR_NOW is needed to handle the unit conversions from ppt, ppm etc.
        ! USE C_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3
        C_AIR_NOW = C_AIR_cc(interp(timevec, CONC_MAT(:,inm_tempK))  .mod. MODS(inm_tempK), interp(timevec, CONC_MAT(:,inm_pres))  .mod. MODS(inm_pres))
        ! =================================================================================================


        ! =================================================================================================
        ! Assign values to input variables

        DO I = 1, N_VARS ! <-- N_VARS will cycle through all variables that user can provide or tamper, and leave zero if no input or mod was provided
            IF ((I==inm_TempK) .or. (MODS(I)%col>0) .or. (MODS(I)%MODE > 0) .or. (ABS(MODS(I)%Shift) > 1d-100)) THEN
                TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I)) .mod. MODS(I)
              ! INDRELAY(I)>0 means that user must have provided a column from an input file; MODS(I)%MODE > 0 means NORMALD is in use
            END IF
        END DO

        ! =================================================================================================

        ! =================================================================================================
        ! Calculate Water vapour pressure and concentration
        CALL WATER(ES,EW,CH_H2O)
        ! =================================================================================================

        ! Chemistry
        IF (Chemistry_flag) THEN

            DO I = 1, N_VARS ! <-- N_VARS will cycle through all variables
                IF (INDRELAY_CH(I)>0) THEN ! <-- this will pick those that were paired in CHECK_INPUT_AGAINST_KPP
                    CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
                END IF
            END DO

            CH_TIME_kpp = MODELTIME%sec
            CH_END_kpp = MODELTIME%sec + MODELTIME%dt_chem

            call BETA(CH_BETA) ! this to properly work, lat, lon and Julian Day need to be defined in INIT_FILE

            Call CHEMCALC(CH_GAS, CH_TIME_kpp, CH_END_kpp, TSTEP_CONC(inm_TempK), TSTEP_CONC(inm_swr), CH_Beta,  &
                          CH_H2O, C_AIR_NOW, TSTEP_CONC(inm_CS), TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2)

            if (model_H2SO4) TSTEP_CONC(inm_H2SO4) = CH_GAS(ind_SA)

        END IF ! IF (Chemistry_flag)

        ! =================================================================================================
        ! NUCLEATION
        IF (NUCLEATION .and. (.not. error%error_state)) THEN
            if (ACDC) THEN
                CALL ACDC_J(TSTEP_CONC)
            else
              continue
              ! CALL SOME_OTHER_NUCLEATION_TYPE
            END if
        END if
        ! =================================================================================================

        ! Condensation
        if (VAP_logical) then

      !      call condensation_routine(values,timestep,end_time,AER_setup, AER_par_prop,ambient,vapours)
        !print*, shape(vapours%c_sat)

        end if
        ! Coagulation
        ! Deposition

        ! =================================================================================================
        ! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
          ! if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
          ! if (MODELTIME%savenow) CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, vapours)
        ! =================================================================================================

        IF (error%error_state) THEN !ERROR handling
            CALL error_handling(error, speed_up,MODELTIME%dt,MODELTIME%sec)

            !restore old state
            !status = status_old

            error%error_state = .false.
        ELSE !There was no error during the actual timestep
            ! =================================================================================================

            ! if (MODELTIME%printnow) print FMT_TIME, MODELTIME%hms ! Print time

            ! =================================================================================================

            ! =================================================================================================
            ! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
            if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
            if (MODELTIME%savenow) CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, vapours)
            ! =================================================================================================

            ! =================================================================================================
            ! Add main timestep to modeltime
            MODELTIME = ADD(MODELTIME)
          ! =================================================================================================

        END IF !ERROR hanldling
    ! =================================================================================================
    END DO	! Main loop ends
    ! =================================================================================================

    CALL PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

    CALL CLOSE_FILES() !Close output file netcdf
    CALL FAREWELL ! Ask if *general.nc is plotted

CONTAINS


  !==================================================================================================
  !
  ! Error handling for aerosol dynamics: In case any error appears during the calculation of
  ! dynamics, this subroutine is called to handle the error. This means adjusting the timestep (either
  ! via speed_up or dt_box) for the process the error showed up
  ! Note that this only works if we calculate "future" changes not what happend in the last time step!
  ! INPUTS:  a) errors; b) speed_up; dt_box
  ! Outputs: a) updated speed_up or dt_box
  !
  !==================================================================================================

  SUBROUTINE error_handling(error,speed_up,dt_box,sim_time_sec)

    IMPLICIT none
    INTEGER ::  i    , &
        speed_up(10)
    REAL(DP) :: dt_box  ,&
        sim_time_sec
    type(error_type) :: error

    !Write some error informaiton to file
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

  SUBROUTINE ACDC_J(C)
    implicit none
    real(dp) :: c_org = 0d0
    real(dp), intent(in) :: C(:)
    real(dp)             :: H2SO4=0,NH3=0,DMA=0,IPR=0 ! these are created to make the unit conversion

    ! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!
    IF (inm_H2SO4 /= 0) H2SO4 = C(inm_H2SO4)*1d6
    IF (inm_NH3   /= 0) NH3   = C(inm_NH3)*1d6
    IF (inm_DMA   /= 0) DMA   = C(inm_DMA)*1d6
    IF (inm_IPR   /= 0) IPR   = C(inm_IPR)*1d6
    ! Speed up program by ignoring nucleation when there is none
    if (NH3 > 1d12 .or. J_ACDC_NH3 > 1d-6) THEN
        CALL get_acdc_J(H2SO4,NH3,c_org,C(inm_CS),C(inm_TEMPK),IPR,MODELTIME,&
            ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)

    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (MODELTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
    END IF
    ! Speed up program by ignoring nucleation when there is none
    if (DMA > 1d6 .or. J_ACDC_DMA > 1d-6) THEN
        CALL get_acdc_D(H2SO4,DMA,c_org,C(inm_CS),C(inm_TEMPK),MODELTIME,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (MODELTIME%printnow) print FMT_SUB, 'DMA IGNORED'
    END IF

  END SUBROUTINE ACDC_J  ! END ACDC Nucleation


  ! ! ================================================================================================
  ! ! This is just a placeholder for the time being and has no function
  ! ! ================================================================================================
  ! SUBROUTINE SOME_OTHER_NUCLEATION_TYPE()
  !   IMPLICIT NONE
  !   REAL(dp) :: J1, J2
  !   J1 = 1d0
  !   J2 = 1d0
  ! END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE


  ! ================================================================================================
  ! PRINT_KEY_INFORMATION for the user. Can be called for example upon MODELTIME%printnow with
  ! if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
  ! ================================================================================================
  SUBROUTINE PRINT_KEY_INFORMATION(C)
    IMPLICIT NONE
    real(dp), intent(in) :: C(:)
    print FMT10_2CVU,'ACID C: ', C(inm_H2SO4), ' [1/cm3]', 'Temp:', C(inm_TempK), 'Kelvin'
    print FMT10_CVU,'APINE C: ', C(IndexFromName('APINENE')), ' [1/cm3]'
    print FMT10_2CVU,'Pressure: ', C(inm_pres), ' [Pa]', 'Air_conc', C_AIR_cc(C(inm_TempK), C(inm_pres)), ' [1/cm3]'
    IF (inm_NH3   /= 0) print FMT10_2CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
    IF (inm_DMA   /= 0) print FMT10_2CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
    print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', C(inm_CS) , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
    print FMT_LEND,
  END SUBROUTINE PRINT_KEY_INFORMATION




  ! =================================================================================================
  ! If user had opted a save or print interval which does not concide with last timestep, we print
  ! and save the values at last timesteps
  ! =================================================================================================
  SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY
      implicit none
      MODELTIME = ADD(MODELTIME, -MODELTIME%dt)
      print*
      print FMT_HDR, 'Main loop ended'

      if (.not. MODELTIME%printnow) THEN
          print FMT_MSG, 'Values from the last timestep:'
          print FMT_TIME, MODELTIME%hms
          CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
      END IF
      if (.not. MODELTIME%savenow)  CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, vapours)
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
                      INDRELAY_CH(I) = J !store the key->value map for input and chemistry
                      exit
                  end if
              END DO
              IF (check == 0 .and. I>LENV) THEN
                  print FMT_FAT0, 'You are using an (organic?) compound which does not exist in chemistry: '//TRIM(MODS(i)%NAME)//' '
                  IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                  IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                  IF (ABS(MODS(I)%SHIFT - 0) > 1d-100) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
                  print FMT_MSG, 'Good bye.'
                  STOP
              END IF
              IF (check == 0 .and. I<LENV .and. &
                  ((I==inm_H2SO4) .or. (I==inm_NH3) .or. (I==inm_DMA) .or. (I==inm_SO2) &
                  .or. (I==inm_NO) .or. (I==inm_NO2) .or. (I==inm_CO) .or. (I==inm_H2) .or. (I==inm_O3))  ) THEN
                  print FMT_WARN0, 'A compound which does not exist in chemistry but could be used elsewhere: '//TRIM(MODS(i)%NAME)//' '
                  IF (MODS(I)%COL > 0) print FMT_SUB, 'In INITFILE; &NML_MODS (col) <- input from file.'
                  IF (MODS(I)%MODE > 0) print FMT_SUB, 'In INITFILE; &NML_MODS - a function for input in use.'
                  IF (ABS(MODS(I)%SHIFT - 0) > 1d-100) print FMT_SUB, 'In INITFILE; &NML_MODS (shift) <- modification of value.'
              END IF
          END IF
      END DO

  END SUBROUTINE CHECK_INPUT_AGAINST_KPP




  ! ============================================================================================================
  ! Calculate water saturation vapour pressure (ES), actual vapour pressure (EW) (in Pa) and concentration (CW)
  ! ============================================================================================================
  SUBROUTINE WATER(ES,EW,CW)
    IMPLICIT NONE
    REAL(dp), INTENT(INOUT) :: ES, EW, CW
    REAL(dp)                :: TEMPC
    TEMPC = TSTEP_CONC(inm_TempK)-273.15d0

    ! Saturation vapour pressure over liquid water; using parametrisation from ???, a0,a1...a6 are in constants.f90
    ES = (a0 + a1 * TEMPC**1 + a2 * TEMPC**2 + a3 * TEMPC**3    &
             + a4 * TEMPC**4 + a5 * TEMPC**5 + a6 * TEMPC**6)* 100

    ! Water vapour pressure
    EW  = TSTEP_CONC(inm_RH) * ES / 100d0

    ! Water vapour concentration in molecules per cm3
    CW = EW/TSTEP_CONC(inm_pres) * C_AIR_NOW

  END SUBROUTINE WATER



  ! ================================================================================================
  ! Calculate Beta parameter (solar angle above horizon), parameterization by Henderson-Sellers
  ! ================================================================================================
  SUBROUTINE BETA(CH_BETA)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: CH_BETA
    real(dp) :: D_zsd,HRANG,ZSR

    D_zsd  = 0.409d0 * sin(2d0 * pi * (MODELTIME%JD-79.8) / 365.24)
    HRANG  = pi * (MODELTIME%sec / 43200. - 1.)
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
    if (for == 0) THEN
      do while(rof == for)
        print FMT_LEND,
        print '(a,20(" "),a)', ACHAR(10),'Press any key to start the run or Ctrl-C to abort'
        read(*,*)
        print FMT_LEND,
        for = -1
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

  SUBROUTINE FAREWELL
    IMPLICIT NONE
    character(1) :: buf

    write(*,*)
    write(*,'(a,1(" "),a)', advance='no') 'SIMULATION HAS ENDED. Plot general output (requires Python3). y? '
    read(*,*) buf
    if (UCASE(buf) == 'Y') CALL EXECUTE_COMMAND_LINE('python3 Scripts/PlotNetCDF.py '//'output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'_general.nc')
    write(*, '(a)') ' SO LONG!'
    write(*,*)
    write(*,*)

  END SUBROUTINE FAREWELL

END PROGRAM SUPERMODEL
