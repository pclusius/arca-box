PROGRAM Supermodel

USE second_Precision,  ONLY : dp ! KPP Numerical type
USE second_Monitor, ONLY: SPC_NAMES
USE constants
USE AUXILLARIES
use INPUT
use OUTPUT
USE ACDC_NH3
USE ACDC_DMA
USE aerosol_auxillaries

IMPLICIT NONE
! MOST OF THE VARIABLES ARE DEFINED IN INPUT.F90
REAL(dp), ALLOCATABLE :: TSTEP_CONC(:)
REAL(DP), ALLOCATABLE :: CH_GAS_DUMMY(:)
REAL(dp) :: J_ACDC_NH3 = 0d0
REAL(dp) :: J_ACDC_DMA = 0d0
REAL(dp) :: J_NH3_BY_IONS(3) = 0d0
REAL(dp) :: acdc_cluster_diam = 2.17d-9

INTEGER  :: I

!speed_up: factor for increasing integration time step for individual prosesses
!...(1): Photolysis, (2): chemistry; (3):Nucleation; (4):Condensation; (5): Coagulation; (6): Deposition
INTEGER :: speed_up(10) = 1
TYPE(error_type) :: error

type(aerosol_setup)       :: AER_setup
type(particle_properties) :: AER_par_prop

  CALL read_input_data ! Declare most variables and read user input and options in input.f90

  CALL CHECK_MODIFIERS ! Print out which modifiers differ from default values

  CALL CHECK_INPUT_AGAINST_KPP

  ALLOCATE(TSTEP_CONC(N_VARS))
  TSTEP_CONC = 0
  ALLOCATE(CH_GAS_DUMMY(size(SPC_NAMES)))
  CH_GAS_DUMMY = 0

  !Open error output file
  open(unit=333, file='output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'_error_output.txt')

  !Open output file
  CALL OPEN_GASFILE(('output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'.nc'), MODS, Description)

  print*,;print FMT_HDR, 'Beginning simulation' ! Information to user

! =================================================================================================
  DO WHILE (MODELTIME%SIM_TIME_S - MODELTIME%sec > -1d-12) ! MAIN LOOP STARTS HERE
! =================================================================================================

  if (MODELTIME%printnow) print FMT_TIME, MODELTIME%hms ! Print time
  C_AIR_NOW = C_AIR_cc(interp(timevec, CONC_MAT(:,inm_tempK))  .mod. MODS(inm_tempK), interp(timevec, CONC_MAT(:,inm_pres))  .mod. MODS(inm_pres))
! =================================================================================================
! Assign values to input variables
  DO I = 1, N_VARS ! <-- N_VARS will cycle through all variables that user can provide or tamper, and leave zero if no input or mod was provided
    IF ((MODS(I)%col>0) .or. (MODS(I)%MODE > 0) .or. (ABS(MODS(I)%Shift-0d0) > 1d-100)) THEN
      TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I))  .mod. MODS(I)
      ! INDRELAY(I)>0 means that user must have provided a column from an input file; MODS(I)%MODE > 0 means NORMALD is in use
    END IF
  END DO
! =================================================================================================

! Chemistry
  DO I = 1, N_VARS ! <-- N_VARS will cycle through all variables
    IF (INDRELAY_CH(I)>0) THEN ! <-- this will pick those that were paired in CHECK_INPUT_AGAINST_KPP
      CH_GAS_DUMMY(INDRELAY_CH(I)) = TSTEP_CONC(I)
    END IF
  END DO

! =================================================================================================
! NUCLEATION
  if (NUCLEATION .and. (.not. error%error_state)) THEN
    if (ACDC) THEN
      CALL ACDC_J(TSTEP_CONC)
    else
      CALL SOME_OTHER_NUCLEATION_TYPE
    END if
  END if
! =================================================================================================

! Condensation
! Coagulation
! Deposition

! =================================================================================================
! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
  ! if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
  ! if (MODELTIME%savenow) CALL SAVE_GASES(TSTEP_CONC(inm_TempK), TSTEP_CONC(inm_H2SO4), TSTEP_CONC(inm_nh3),&
  !                             TSTEP_CONC(inm_dma),J_ACDC_NH3, J_ACDC_DMA, TSTEP_CONC(inm_cs), TSTEP_CONC(inm_pres), MODS)
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
    if (MODELTIME%savenow) CALL SAVE_GASES(TSTEP_CONC(inm_TempK), TSTEP_CONC(inm_H2SO4), TSTEP_CONC(inm_nh3),&
                                TSTEP_CONC(inm_dma),J_ACDC_NH3, J_ACDC_DMA, TSTEP_CONC(inm_cs), TSTEP_CONC(inm_pres), MODS)
  ! =================================================================================================

  ! =================================================================================================
    ! Add main timestep to modeltime
    MODELTIME = ADD(MODELTIME)
  ! =================================================================================================

END IF !ERROR hanldling
! =================================================================================================
  END DO	! Main loop ends
! =================================================================================================

  !CALL PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

  CALL CLOSE_FILES() !Close output file netcdf

  print *, ACHAR(10)," SIMULATION HAS ENDED. SO LONG!",ACHAR(10)

  if (VAP_logical) then
  call allocate_and_setup(AER_setup, AER_par_prop, vapours)
  end if

CONTAINS

  SUBROUTINE error_handling(error,speed_up,dt_box,sim_time_sec)
  !==================================================================================================
  !Error handling for aerosol dynamics: In case any error appears during the calculation of
  !dynamics, this subroutine is called to handle the error. This means adjusting the timestep (either
  ! via speed_up or dt_box) for the process the error showed up
  !Note that this only works if we calculate "future" changes not what happend in the last time step!
  !INPUTS:  a) errors; b) speed_up; dt_box
  !Outputs: a) updated speed_up or dt_box
  !==================================================================================================
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


SUBROUTINE ACDC_J(C)
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
! =================================================================================================
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

  ! ================================================================================================
  ! END ACDC Nucleation
  ! ================================================================================================
END SUBROUTINE ACDC_J

SUBROUTINE SOME_OTHER_NUCLEATION_TYPE()
  IMPLICIT NONE
  REAL(dp) :: J1, J2
  J1 = 1d0
  J2 = 1d0
END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE

SUBROUTINE PRINT_KEY_INFORMATION(C)
  IMPLICIT NONE
  real(dp), intent(in) :: C(:)
  print FMT10_2CVU,'ACID C: ', C(inm_H2SO4), ' [1/cm3]', 'Temp:', C(inm_TempK), 'Kelvin'
  print FMT10_CVU,'APINE C: ', C(IndexFromName('APINENE')), ' [1/cm3]'
  print FMT10_2CVU,'Pressure: ', C(inm_pres), ' []', 'Air_conc', C_AIR_cc(C(inm_TempK), C(inm_pres)), ' [1/cm3]'
  IF (inm_NH3   /= 0) print FMT10_2CVU, 'NH3 C:', C(inm_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
  IF (inm_DMA   /= 0) print FMT10_2CVU, 'DMA C:', C(inm_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
  print FMT10_3CVU, 'Jion neutral:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion neg:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion pos:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
  IF (inm_IPR   /= 0) print FMT10_2CVU, 'C-sink:', C(inm_CS) , ' [1/s]','IPR:', C(inm_IPR) , ' [1/s/cm3]'
  print FMT_LEND,
END SUBROUTINE PRINT_KEY_INFORMATION

SUBROUTINE CHECK_MODIFIERS()
  IMPLICIT NONE
  type(input_mod) :: test
  integer         :: i,j=0
  character(4)    :: cprf

  print FMT_HDR, 'Check input validity'
  CALL CONVERT_TEMPS_TO_KELVINS
  CALL CONVERT_PRESSURE
  do i=1,size(MODS)
    IF (MODS(i)%MODE > 0) THEN
      print FMT_NOTE0, 'Replacing input for '//TRIM(MODS(i)%name)//' with parametrized function.'
      j=1
    ELSE
      IF (ABS(MODS(i)%multi - test%multi) > 1d-9) THEN
        print FMT_NOTE1, 'Multiplying '//TRIM(MODS(i)%name)//' with: ',MODS(i)%multi
        j=1
      END IF
      if (ABS(MODS(i)%shift - test%shift) > 1d-9) THEN
        if (TRIM(MODS(i)%UNIT) == '#') THEN
          cprf = '/cm3'
        else
          cprf = ''
        end if
        print FMT_NOTE1, 'Adding a constant to '//TRIM(MODS(i)%name)//', [in '//TRIM(MODS(i)%UNIT)//TRIM(cprf)//']: ',MODS(i)%shift
        j=1
      END IF
    END IF
  END DO
  if (j == 1) print FMT_LEND
END SUBROUTINE CHECK_MODIFIERS

SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY
! =================================================================================================
! If user had opted a save or print interval which does not concide with last timestep, we print
! and save the values at last timesteps
  implicit none
  MODELTIME = ADD(MODELTIME, -MODELTIME%dt)
  print*
  print FMT_HDR, 'Main loop ended'

  if (.not. MODELTIME%printnow) THEN
    print FMT_MSG, 'Values from the last timestep:'
    print FMT_TIME, MODELTIME%hms
    CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
  END IF
  if (.not. MODELTIME%savenow) call SAVE_GASES(TSTEP_CONC(inm_tempK), TSTEP_CONC(inm_H2SO4), &
  TSTEP_CONC(inm_nh3), TSTEP_CONC(inm_dma),J_ACDC_NH3, J_ACDC_DMA, TSTEP_CONC(inm_cs), TSTEP_CONC(inm_pres), MODS)
END SUBROUTINE PRINT_FINAL_VALUES_IF_LAST_STEP_DID_NOT_DO_IT_ALREADY

SUBROUTINE CHECK_INPUT_AGAINST_KPP
  implicit none
  integer :: i,j, check
  print FMT_HDR, 'Checking against KPP for chemicals'
  do i=1,N_VARS
    check = 0
    IF (MODS(I)%col > 0 .or. MODS(I)%MODE > 0 .or. ABS(MODS(I)%SHIFT - 0) > 1d-100) THEN
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

SUBROUTINE CONVERT_TEMPS_TO_KELVINS
  use constants, ONLY: UCASE
  IMPLICIT NONE

  if ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == '#') THEN
    print FMT_WARN0, "No unit for temperature. Use either 'K' or 'C'. Now assuming Kelvins."
    TempUnit = 'K'
  elseif ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'K') THEN
    TempUnit = 'K'
  elseif ((TRIM(UCASE(TempUnit)) /= 'K' .and. TRIM(UCASE(TempUnit)) /= 'C') .and. TRIM(UCASE(MODS(inm_TempK)%UNIT)) == 'C') THEN
    TempUnit = 'C'
  END IF

  IF (UCASE(TempUnit) == 'K') THEN
    print FMT_MSG, '- Temperature input in Kelvins.'
  ELSEIF (UCASE(TempUnit) == 'C') THEN
    print FMT_MSG, '- Converting temperature from degrees C -> K.'
    MODS(inm_TempK)%min = MODS(inm_TempK)%min + 273.15d0
    MODS(inm_TempK)%max = MODS(inm_TempK)%max + 273.15d0
    CONC_MAT(:,inm_TempK) = CONC_MAT(:,inm_TempK) + 273.15d0
  ELSE
    print FMT_WARN0, "Could not recognize temperature unit. Use either 'K' or 'C'. Now assuming Kelvins."
    TempUnit = 'K'
  END IF
  IF ((TempUnit == 'C') .and. (  ABS(MODS(inm_TempK)%shift - 273.15)<1d0  )) THEN
    print FMT_WARN1, 'Temperature will be converted to Kelvins, but a suspicious constant is added: ',MODS(inm_TempK)%shift
  END IF
  MODS(inm_TempK)%UNIT = TempUnit
END SUBROUTINE CONVERT_TEMPS_TO_KELVINS


SUBROUTINE CONVERT_PRESSURE
  use constants, ONLY: UCASE
  IMPLICIT NONE
  character(5) :: buf

  buf = UCASE(TRIM(MODS(inm_pres)%UNIT))
  if (TRIM(buf) == 'HPA' .or. TRIM(buf) == 'MBAR') THEN
    CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 100d0
    print FMT_MSG, '- Converting pressure from hPa (mbar) to Pascals.'
  elseif (TRIM(buf) == 'KPA') THEN
    CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 1000d0
    print FMT_MSG, '- Converting pressure from kPa to Pascals.'
  elseif (TRIM(buf) == 'ATM') THEN
    CONC_MAT(:,inm_pres) = CONC_MAT(:,inm_pres) * 1.01325d5
    print FMT_MSG, '- Converting pressure from atm to Pascals.'
  elseif (TRIM(buf) == 'PA') THEN
    print FMT_MSG, '- Pressure in Pascals.'
    continue
  else
    if ((MODS(inm_pres)%MODE > 0) .or. (MODS(inm_pres)%col > 1)  .or. (ABS(MODS(inm_pres)%multi - 1d0)>1d-9) .or. (ABS(MODS(inm_pres)%shift)>1d-16)) THEN
      if (TRIM(buf) == '#') THEN
        print FMT_MSG, '- Assuming Pascals for pressure.'
      else
        print FMT_FAT0, 'Cannot recognize given unit "'//TRIM(MODS(inm_pres)%UNIT)//'" for pressure. Exiting. '
        stop
      end if
    end if
  end if

  do i=1,N_VARS
    buf = UCASE(TRIM(MODS(i)%UNIT))
    if (TRIM(buf) /= '#' .and. i /= inm_pres .and. i /= inm_tempK) THEN
      IF (TRIM(buf) /= 'PPM' .and. TRIM(buf) /= 'PPB' .and. TRIM(buf) /= 'PPT' .and. TRIM(buf) /= 'PPQ') THEN
        print FMT_FAT0, 'Cannot recognize unit "'//TRIM(MODS(i)%UNIT)//'" for '//TRIM(MODS(i)%name)//'. Exiting. '
        stop
      ELSE
        if ((MODS(i)%MODE > 0) .or. (MODS(i)%col > 1)  .or. (ABS(MODS(i)%multi - 1d0)>1d-9) .or. (ABS(MODS(i)%shift)>1d-16)) THEN
          print FMT_MSG, '- Converting '//TRIM(MODS(i)%name)//' from '//TRIM(MODS(i)%UNIT)
        end if
      END IF
    end if
  end do
END SUBROUTINE CONVERT_PRESSURE






END PROGRAM SUPERMODEL
