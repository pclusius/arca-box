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
    ! use omp_lib


    IMPLICIT NONE

    ! Transient variables
    CHARACTER(90) :: buf
    CHARACTER(1000) :: inibk
    CHARACTER(:), allocatable:: RUN_OUTPUT_DIR
    INTEGER       :: in_channels ! number of diameters -> passed over to the fitting subroutine
    INTEGER       :: I, J, val1, ioi

    ! MOST OF THE VARIABLES ARE DEFINED IN INPUT.F90
    REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)
    REAL(DP), ALLOCATABLE :: CH_GAS(:)
    REAL(DP), ALLOCATABLE :: shitarray(:)

    REAL(dp) :: CH_RO2      ! RO2 concentration in [molecules / cm^3]
    REAL(dp) :: CH_H2O      ! H20 concentration in [molecules / cm^3]
    REAL(dp) :: CH_Beta     ! solar zenit angle
    REAL(dp) :: EW          ! Water content in Pa
    REAL(dp) :: ES          ! Saturation vapour pressure
    REAL(dp) :: cpu1, cpu2  ! CPU time in seconds
    REAL(dp) :: Tavg = 0d0, Pavg = 0d0 ! Average Temperature and Pressure

    !speed_up: factor for increasing integration time step for individual prosesses
    !...(1): Photolysis, (2): chemistry; (3):Nucleation; (4):Condensation; (5): Coagulation; (6): Deposition
    INTEGER          :: speed_up(10) = 1
    INTEGER          :: iters
    TYPE(error_type) :: error

    INTEGER, DIMENSION(:), allocatable :: index_cond
    REAL(dp), dimension(:), allocatable:: conc_vapour

    CALL READ_INPUT_DATA ! Declare most variables and read user input and options in input.f90

    ! set up number of ACDC iterations.
    if (ACDC_solve_ss) THEN
      iters = 1
    ELSE
      iters = acdc_iterations
    end if

    IF (Chemistry_flag) CALL CHECK_INPUT_AGAINST_KPP ! Check that the input exists in chemistry, or if not, print warning



    ! Particles are considered -> initialize a particle representation, set initial PSD and determine composition
    IF (Aerosol_flag) THEN

      ! Initialzie the Particle representation
      CALL INITIALIZE_PSD
      ALLOCATE(shitarray(n_bins_particle))
      ! Initialize the saturation vapour concentration and diffusivity with mean temperature and pressure
      DO i = 2, size(CONC_MAT,1)-1
        Tavg = Tavg + (CONC_MAT(i,inm_TempK) .mod. MODS(inm_Tempk))
        Pavg = Pavg + (CONC_MAT(i,inm_PRES) .mod. MODS(inm_PRES))
      END DO

      Tavg = Tavg/(size(CONC_MAT,1)-2)
      Pavg = Pavg/(size(CONC_MAT,1)-2)

      CALL Calculate_SaturationVapourConcentration(VAPOUR_PROP, Tavg)

      if (.not. temp_depend_csat) THEN
        write(buf, '(a,f6.2, a, f8.1,a)') 'Using mean temperature ', Tavg, ' K and pressure ', Pavg, ' [Pa] for condensation'
        print FMT_MSG, buf
      END IF

      ! Send par_data (from input): diameter and first time step for fitting of initial model PSD
      IF (CURRENT_PSD%PSD_style == 1) THEN !only defined procedure for fully stationary
        print*, 'CHECK COMPO INITIALOIZATION'
          ! intialize concentration of condensables in each bin kg/m3
          do  i = 1, CURRENT_PSD%nr_bins
            CURRENT_PSD%composition_fs(i,:) = VAPOUR_PROP%mfractions * CURRENT_PSD%volume_fs(i) * CURRENT_PSD%conc_fs(i)* VAPOUR_PROP%density * 1d6
            conc_pp(i,:) = CURRENT_PSD%composition_fs(i,:) / VAPOUR_PROP%molar_mass*Na *1d6
          end do


        print*,
        print FMT_HDR, 'INITIALIZING PARTICLE STRUCTURES '

        ! Derive composition of the particles form input (XTRAS(I), I...# of noncond (nr_noncond))
        IF (extra_particles /= '') THEN
          PRINT FMT_MSG,'initial particles are composed of:'
          DO I = 1,size(xtras(:))
            write(buf, '(a,3(es12.3))') xtras(I)%name,xtras(I)%options
            PRINT FMT_MSG, TRIM(buf)

            in_channels = size(XTRAS(I)%sections(:))

            CALL GeneratePSDfromInput(xtras(I)%sections(:),xtras(I)%binseries(1,:), xtras(i)%conc_modelbins)

          END DO
        END IF
      END IF
      print FMT_LEND,
    END IF

    ALLOCATE(TSTEP_CONC(N_VARS))
    TSTEP_CONC = 0

    ALLOCATE(CH_GAS(size(SPC_NAMES)))
    CH_GAS = 0

    ! This only called once for KPP in the beginning
    IF (Chemistry_flag .EQV. .true.) THEN
      CALL KPP_SetUp
    ENDIF



    ! All run output goes to this directory. RUN_OUTPUT_DIR is allocated to correct length so we dont need to TRIM every time
    ALLOCATE(CHARACTER(len=LEN(TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(INDEX)//'/'//TRIM(RUN_NAME))) :: RUN_OUTPUT_DIR)
    RUN_OUTPUT_DIR = TRIM(INOUT_DIR)//'/'//TRIM(CASE_NAME)//'_'//TRIM(DATE)//TRIM(INDEX)//'/'//TRIM(RUN_NAME)

    OPEN(100,file=RUN_OUTPUT_DIR//"/diameter.dat",status='replace',action='write')
    OPEN(104,file=RUN_OUTPUT_DIR//"/time.dat",status='replace',action='write')
    OPEN(101,file=RUN_OUTPUT_DIR//"/particle_conc.sum",status='replace',action='write')
    OPEN(105,file=RUN_OUTPUT_DIR//"/avail_species.dat",status='replace',action='write')
    OPEN(106,file=RUN_OUTPUT_DIR//"/index_cond.dat",status='replace',action='write')
    OPEN(107,file=RUN_OUTPUT_DIR//"/index_avail_comp.dat",status='replace',action='write')

    !Open output file
    write(*,*)
    print FMT_HDR, 'INITIALIZING OUTPUT '
    CALL OPEN_FILES( RUN_OUTPUT_DIR, Description, MODS, CH_GAS, VAPOUR_PROP, CURRENT_PSD)

    !Open error output file
    open(unit=333, file=RUN_OUTPUT_DIR//'/error_output.txt')

    ! Save backup from the initfile
    OPEN(UNIT=504, FILE=TRIM(ADJUSTL(Fname_init)), STATUS='OLD', ACTION='READ', iostat=ioi)
    open(unit=334, file=RUN_OUTPUT_DIR//'/InitBackup.txt')
    ioi = 0
    DO while (ioi == 0)
      read(504,'(a)', iostat = ioi) inibk
      write(334,'(a)') TRIM(inibk)
    end do
    close(504)
    close(334)

    ! reading the index of compounds
    ALLOCATE(index_cond(nr_cond))
    ALLOCATE(conc_vapour(nr_species_P))

    index_cond=0
    val1=0
! check how many species we have in common
    DO i = 1,size(SPC_NAMES)
      DO j = 1,nr_cond
        IF (VAPOUR_PROP%vapour_names(j) .eq. SPC_NAMES(i)) THEN ! SPC_NAMES from second-Monitor
          index_cond(j) = i
          val1 = val1+1
          exit
        END IF
      END DO
    END DO

    write(106,'(i0)') index_cond
    write(107,'(i0)') pack(index_cond,index_cond/=0)
    write(105,'(a)') spc_names(pack(index_cond,index_cond/=0))
    CLOSE(105)
    CLOSE(106)
    CLOSE(107)


    CALL PAUSE_FOR_WHILE(wait_for)

    print*,;print FMT_HDR, 'Beginning simulation' ! Information to user
    call cpu_time(cpu1)
    ! =================================================================================================
    DO WHILE (GTIME%SIM_TIME_S - GTIME%sec > -1d-12) ! MAIN LOOP STARTS HERE
        ! =================================================================================================


        ! =================================================================================================
        ! Print time header in a nice way, start with empty row
        if (GTIME%printnow) THEN
          print *, ! Print time
          print FMT_TIME, GTIME%hms
          call cpu_time(cpu2)
        ! if --gui flag was used, print a dot and EOL so the STDOUT-reading in Python GUI will be smoother.
        ELSE
          if (ingui) print'(a)', '.'
        END IF
        ! =================================================================================================
        ! GC_AIR_NOW, TEMPK, PRES and RH are calculated as global variables and are available everywhere.
        ! USE GC_AIR_NOW FOR CURRENT AIR CONCENTRATION IN CM^3
        GTEMPK = interp(timevec, CONC_MAT(:,inm_tempK)) .mod. MODS(inm_tempK)
        GPRES = interp(timevec, CONC_MAT(:,inm_PRES)) .mod. MODS(inm_PRES)
        GRH = interp(timevec, CONC_MAT(:,inm_RH)) .mod. MODS(inm_RH)
        GC_AIR_NOW = C_AIR_cc(GTEMPK, GPRES)
        ! =================================================================================================


        ! =================================================================================================
        ! Assign values to input variables

        DO I = 1, N_VARS ! <-- N_VARS will cycle through all variables that user can provide or tamper, and leave zero if no input or mod was provided
          ! IF ((I==inm_TempK) .or. (MODS(I)%col>0) .or. (MODS(I)%MODE > 0) .or. (ABS(MODS(I)%Shift) > 1d-100)) THEN
              TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I)) .mod. MODS(I)
            ! INDRELAY(I)>0 means that user must have provided a column from an input file; MODS(I)%MODE > 0 means NORMALD is in use
          ! END IF
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
              IF (I == inm_H2SO4 .and. .not. model_H2SO4) THEN
                CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
              ELSE IF (I /= inm_H2SO4) THEN
                CH_GAS(INDRELAY_CH(I)) = TSTEP_CONC(I)
              END IF
            END IF
          END DO

          ! Solar Zenith angle. For this to properly work, lat, lon and Date need to be defined in INIT_FILE
          call BETA(CH_Beta)

            ! write(*,*) 'L 262 before chemcalc sum of ch_gas', sum(CH_GAS)
          Call CHEMCALC(CH_GAS, GTIME%sec, (GTIME%sec + GTIME%dt_chm), GTEMPK, TSTEP_CONC(inm_swr), CH_Beta,  &
                        CH_H2O, GC_AIR_NOW, TSTEP_CONC(inm_CS), TSTEP_CONC(inm_CS_NA), CH_Albedo, CH_RO2)

          if (model_H2SO4) TSTEP_CONC(inm_H2SO4) = CH_GAS(ind_H2SO4)

        END IF ! IF (Chemistry_flag)

        ! =================================================================================================
        ! NUCLEATION
        IF (NUCLEATION .and. (.not. error%error_state) ) THEN
            if (ACDC) THEN
                CALL ACDC_J(TSTEP_CONC, iters)
                J_TOTAL = J_ACDC_DMA + J_ACDC_NH3 ! [particles/s/m^3]
                IF (.not. RESOLVE_BASE) J_TOTAL = J_TOTAL + TSTEP_CONC(inm_JIN) * 1d6 ! [particles/s/cm^3]
            else
              continue
              ! CALL SOME_OTHER_NUCLEATION_TYPE
            END if
        END if
        if (GTIME%savenow .and. RESOLVE_BASE) CALL Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)


        ! =================================================================================================
        ! AEROSOL PART STARTS HERE
        ! =================================================================================================
        IF (Aerosol_flag) THEN


          ! =================================================================================================
          ! Read in background particles
          if (use_dmps .and. GTIME%min/60d0 < DMPS_read_in_time) THEN

            IF (GTIME%printnow) print FMT_SUB, 'Reading in background particles'

            ! Initialization is necessary because GeneratePSDfromInput skips some bins and this would lead to blowup
            CURRENT_PSD%conc_fs = 0d0
            CALL GeneratePSDfromInput( par_data(1,2:),  par_data(min(GTIME%ind_netcdf+1, size(par_data, 1)),2:), CURRENT_PSD%conc_fs )
            where(CURRENT_PSD%conc_fs<1d-28) CURRENT_PSD%conc_fs = 1d-28
            ! NOTE Sumfile is typically in particles /cm^3
            CURRENT_PSD%conc_fs = CURRENT_PSD%conc_fs * 1d6

          END IF

          val1=0
          do i = 1, nr_cond
            if (index_cond(i) /= 0) then
            ! print*, i, index_cond(i)
              conc_vapour(i) =  CH_GAS(index_cond(i))*1D6 ! mol/m3
              val1=val1+1
            else
              conc_vapour(i) = 0d0
            end if
          end do

          conc_vapour(nr_species_P) = CH_GAS(ind_H2SO4)*1d6

          call Nucleation_routine(J_TOTAL,CURRENT_PSD%conc_fs)

          if (Condensation) THEN
            ! CURRENT_PSD%conc_fs = CURRENT_PSD%conc_fs*1d6
            ! CURRENT_PSD%composition_fs = CURRENT_PSD%composition_fs*1d6
            ! =================================================================================================
            if (temp_depend_csat) CALL Calculate_SaturationVapourConcentration(VAPOUR_PROP, GTEMPK)
            CALL Condensation_apc(VAPOUR_PROP,CURRENT_PSD,conc_pp,conc_vapour,dmass)
            ! print*, 'bef mnc conc_pp', sum(conc_pp)
            ! print*, 'bef mnc compo', sum(sum(CURRENT_PSD%composition_fs(:,:),1) * Na / VAPOUR_PROP%molar_mass)
            CALL Mass_Number_Change('condensation')

          end if ! end of condensation

          do i = 1, nr_cond
            if (index_cond(i) /= 0) then
             CH_GAS(index_cond(i)) = conc_vapour(i) *1D-6
            end if
          end do
          CH_GAS(ind_H2SO4) = conc_vapour(nr_species_P)*1d-6

          ! For coagulation
          if (Coagulation) then
            Call Coagulation_routine(GTIME%dt_aer,CURRENT_PSD,dconc_coag)
            Call Mass_Number_Change('coagulation')
          end if

          do i = 1,current_PSD%nr_bins
            conc_pp(i,:) = CURRENT_PSD%composition_fs(i,:) * Na / VAPOUR_PROP%molar_mass
          end do
          ! CURRENT_PSD%conc_fs = CURRENT_PSD%conc_fs*1d-6
          ! CURRENT_PSD%composition_fs = CURRENT_PSD%composition_fs*1d-6

        end if ! if Aerosol_flag


        ! =================================================================================================
        ! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
        ! if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
        ! if (GTIME%savenow) CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, VAPOUR_PROP)
        ! =================================================================================================

        IF (error%error_state) THEN !ERROR handling
          CALL error_handling(error, speed_up,GTIME%dt,GTIME%sec)
          error%error_state = .false.

        ! There was no error during the actual timestep
        ELSE
          ! ============================================================================================================
          ! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
          if (GTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
          if (GTIME%savenow) THEN
            CALL SAVE_GASES(TSTEP_CONC,MODS,CH_GAS,J_ACDC_NH3, J_ACDC_DMA, VAPOUR_PROP, CURRENT_PSD)
            WRITE(100,*) CURRENT_PSD%diameter_fs
            WRITE(101,*) CURRENT_PSD%conc_fs*1d-6
            WRITE(104,*) GTIME%sec
          END IF

          ! ============================================================================================================
          ! Add main timestep to GTIME'
          GTIME = ADD(GTIME)

        END IF ! ERROR handling

    END DO	! Main loop ends
    ! ==================================================================================================================
    ! ==================================================================================================================

    CLOSE(100)
    CLOSE(101)
    CLOSE(104)


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
            ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
      end do
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
    END IF
    ! Speed up program by ignoring nucleation when there is none
    if ((DMA > 1d6 .and. H2SO4>1d9) .or. (.not. skip_acdc)) THEN
      ! Idea of iteration is to not necessarily reach steady state but still get a "better" estimate of instantenous formation rate
      Do i=1,iters
        CALL get_acdc_D(H2SO4,DMA,c_org,C(inm_CS),C(inm_TEMPK),GTIME,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
      end do
    ELSE
        ! This will leave the last value for J stand - small enough to not count but not zero
        if (GTIME%printnow) print FMT_SUB, 'DMA IGNORED'
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
    IF (inm_DMA   /= 0) print FMT10_2CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
    print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
    IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', C(inm_CS) , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
    if ((GTIME%sec)>0) print '("| ",a,i0,a,i0.2,a,i0,a,i0.2,t65,a,f6.2,t100,"|")', 'Elapsed time (m:s) ',int(cpu2 - cpu1)/60,':',modulo(int(cpu2 - cpu1),60) ,' Est. time to finish (m:s) ',&
                            int((cpu2 - cpu1)/GTIME%sec*(GTIME%SIM_TIME_S-GTIME%sec))/60,':', MODULO(int((cpu2 - cpu1)/GTIME%sec*(GTIME%SIM_TIME_S-GTIME%sec)),60),&
                            'Realtime/Modeltime: ', GTIME%sec/(cpu2 - cpu1)

    print FMT_LEND,
  END SUBROUTINE PRINT_KEY_INFORMATION




  ! =================================================================================================
  ! If user had opted a save or print interval which does not concide with last timestep, we print
  ! and save the values at last timesteps
  ! =================================================================================================
  SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY
      implicit none
      GTIME = ADD(GTIME, -GTIME%dt)
      print*
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
                      INDRELAY_CH(I) = J !store the 'key->value map' for input and chemistry
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
    ! character(1) :: buf

    write(*,*)
    ! IF (python) THEN
    !   write(*,'(a,1(" "),a)', advance='no') 'SIMULATION HAS ENDED. Plot general output (requires Python3). y? '
    !   read(*,*) buf
    !   if (UCASE(buf) == 'Y') CALL EXECUTE_COMMAND_LINE('python3 Scripts/PlotNetCDF.py '//'output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'_general.nc')
    ! ELSE
      write(*,'(a,1(" "),a)', advance='no') 'SIMULATION HAS ENDED. '
    ! END IF
    write(*, '(a)') 'SO LONG!'
    write(*,*)
    write(*,*)
    write(333,*) '---------- SIMULATION REACHED END SUCCESFULLY ------------'
    close(333)

  END SUBROUTINE FINISH




END PROGRAM SUPERMODEL
