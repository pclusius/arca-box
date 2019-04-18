PROGRAM Supermodel

USE second_Precision,  ONLY : dp ! KPP Numerical type
USE constants
USE AUXILLARIES
use INPUT
use OUTPUT
#ifdef ISACDC
  USE ACDC_NH3
  USE ACDC_DMA
#endif

IMPLICIT NONE
! MOST OF THE VARIABLES ARE DEFINED IN INPUT.F90
REAL(dp) :: TSTEP_CONC(ind_LAST) = 0d0
REAL(dp) :: J_ACDC_NH3
REAL(dp) :: J_ACDC_DMA
REAL(dp) :: J_NH3_BY_IONS(3)
REAL(dp) :: acdc_cluster_diam = 2.17d-9
INTEGER :: I

! =================================================================================================
  print*, ! PROGRAM STARTS with an empty line
! =================================================================================================

! =================================================================================================
  CALL read_input_data() ! READ USER INPUT AND OPTIONS
! =================================================================================================

! =================================================================================================
  CALL OPEN_GASFILE(('output/'//TRIM(CASE_NAME)//'_'//TRIM(RUN_NAME)//'.nc'), MODS, Description)
! =================================================================================================

! =================================================================================================
  CALL CHECK_MODIFIERS() ! Print out which modifiers differ from default values
! =================================================================================================

! =================================================================================================
  print*,;print FMT_HDR, 'Beginning simulation' ! Information to user
! =================================================================================================
  DO WHILE (MODELTIME%sec < MODELTIME%SIM_TIME_S + 1d-12) ! MAIN LOOP STARTS HERE
! =================================================================================================

! =================================================================================================
  if (MODELTIME%printnow) print FMT_TIME, MODELTIME%hms ! Print time
! =================================================================================================

! =================================================================================================
! Assign values to input variables
  DO I = 1, IND_LAST ! <-- IND_last will cycle through all variables that user can provide or tamper, and leave zero if no input or mod was provided
    IF ((INDRELAY(I,2)>0) .or. (MODS(I)%MODE > 0) .or. (ABS(MODS(I)%Shift-0d0) > 1d-100)) THEN
      TSTEP_CONC(I) = interp(timevec, CONC_MAT(:,I))  .mod. MODS(I)
      ! INDRELAY(I,2)>0 means that user must have provided a column from an input file; MODS(I)%MODE > 0 means NORMALD is in use
    END IF
  END DO
! =================================================================================================

! =================================================================================================
! NUCLEATION
  if (NUCLEATION) THEN
#ifdef ISACDC
    if (ACDC) CALL ACDC_J(TSTEP_CONC) ! SUBROUTINE in CONTAINS of this file
#endif
  END if
! =================================================================================================

! Condensation
! Coagulation
! Deposition

! =================================================================================================
! Write printouts to screen and outputs to netcdf-file, later this will include more optionality
  if (MODELTIME%printnow) CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
  if (MODELTIME%savenow) CALL SAVE_GASES(TSTEP_CONC(ind_temp), TSTEP_CONC(ind_H2SO4), TSTEP_CONC(ind_nh3), TSTEP_CONC(ind_dma),&
                                          J_ACDC_NH3, J_ACDC_DMA, TSTEP_CONC(ind_cs), MODS)
! =================================================================================================

! =================================================================================================
  MODELTIME = ADD(MODELTIME) ! Add main timestep to modeltime
! =================================================================================================

! =================================================================================================
  END DO	! Main loop ends
! =================================================================================================

! =================================================================================================
! If user had opted a save or print interval which does not concide with last timestep, we print
! and save the values at last timesteps
  MODELTIME = ADD(MODELTIME, -MODELTIME%dt)
  print*
  print FMT_HDR, 'Main loop ended'
  print FMT_MSG, 'Some housekeeping...'

  if (.not. MODELTIME%printnow) THEN
    print FMT_MSG, 'Values from the last timestep:'
    print FMT_TIME, MODELTIME%hms
    CALL PRINT_KEY_INFORMATION(TSTEP_CONC)
  END IF
  if (.not. MODELTIME%savenow) call SAVE_GASES(TSTEP_CONC(ind_temp), TSTEP_CONC(ind_H2SO4), &
  TSTEP_CONC(ind_nh3), TSTEP_CONC(ind_dma),J_ACDC_NH3, J_ACDC_DMA, TSTEP_CONC(ind_cs), MODS)
! =================================================================================================

! =================================================================================================
  CALL CLOSE_FILES() !Close output file netcdf
! =================================================================================================

! =================================================================================================
  print *, ACHAR(10)," SIMULATION HAS ENDED. SO LONG!",ACHAR(10)
! =================================================================================================

CONTAINS

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
real(dp) :: c_org = 1d-12
real(dp), intent(in) :: C(:)
real(dp)             :: H2SO4,NH3,DMA,IPR ! these are crated to make the unit conversion

! NUCLEATION BY S-ACID AND NH3 - NOTE: ingoing concentrations are assumed to be in 1/m3!!
  H2SO4 = C(IND_H2SO4)*1d6
  NH3   = C(IND_NH3)*1d6
  DMA   = C(IND_DMA)*1d6
  IPR   = C(IND_IPR)*1d6
  ! Speed up program by ignoring nucleation when there is none
  if (NH3 > 1d12 .or. J_ACDC_NH3 > 1d-6) THEN
    CALL get_acdc_J(H2SO4,NH3,c_org,C(IND_CS),C(IND_TEMP),IPR,MODELTIME,&
          ACDC_solve_ss,J_ACDC_NH3,acdc_cluster_diam, J_NH3_BY_IONS)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (MODELTIME%printnow) print FMT_SUB, 'NH3 IGNORED'
  END IF

  ! Speed up program by ignoring nucleation when there is none
  if (DMA > 1d6 .or. J_ACDC_DMA > 1d-6) THEN
    CALL get_acdc_D(H2SO4,DMA,c_org,C(IND_CS),C(IND_TEMP),MODELTIME,ACDC_solve_ss,J_ACDC_DMA,acdc_cluster_diam)
  ELSE
    ! This will leave the last value for J stand - small enough to not count but not zero
    if (MODELTIME%printnow) print FMT_SUB, 'DMA IGNORED'
  END IF

  ! ================================================================================================
  ! END ACDC Nucleation
  ! ================================================================================================
END SUBROUTINE ACDC_J

SUBROUTINE SOME_OTHER_NUCLEATION_TYPE(J1, J2)
  IMPLICIT NONE
  REAL(dp) :: J1, J2
  J1 = 1d0
  J2 = 1d0
END SUBROUTINE SOME_OTHER_NUCLEATION_TYPE

SUBROUTINE PRINT_KEY_INFORMATION(C)
  IMPLICIT NONE
  real(dp), intent(in) :: C(:)
  print FMT10_2CVU,'ACID C: ', C(IND_H2SO4), ' [1/cm3]', 'Temp:', C(IND_TEMP), 'Kelvin'
  print FMT10_2CVU, 'NH3 C:', C(IND_NH3), ' [1/cm3]','J_NH3:', J_ACDC_NH3*1d-6, ' [1/cm3]'
  print FMT10_2CVU, 'DMA C:', C(IND_DMA) , ' [1/cm3]','J_DMA:', J_ACDC_DMA*1d-6, ' [1/cm3]'
  print FMT10_3CVU, 'Jion1:', J_NH3_BY_IONS(1)*1d-6 , ' [1/s/cm3]','Jion2:', J_NH3_BY_IONS(2)*1d-6 , ' [1/s/cm3]','Jion2:', J_NH3_BY_IONS(3)*1d-6 , ' [1/s/cm3]'
  print FMT10_2CVU, 'C-sink:', C(IND_CS) , ' [1/s]','IPR:', C(IND_IPR) , ' [1/s/cm3]'
  print FMT_LEND,
END SUBROUTINE PRINT_KEY_INFORMATION

SUBROUTINE CHECK_MODIFIERS()
  IMPLICIT NONE
  type(input_mod) :: test
  integer         :: i,j=0
  do i=1,size(MODS)
  IF (MODS(i)%MODE > 0) THEN
    print FMT_WARN0, 'Replacing input for '//TRIM(MODS(i)%name)//' with parametrized function.'
    j=1
  ELSE
    IF (ABS(MODS(i)%multi - test%multi) > 1d-9) THEN
      print FMT_WARN1, 'Multiplying '//TRIM(MODS(i)%name)//' with: ',MODS(i)%multi
      j=1
    END IF
    if (ABS(MODS(i)%shift - test%shift) > 1d-9) THEN
      print FMT_WARN1, 'Adding a constant to '//TRIM(MODS(i)%name)//', value is: ',MODS(i)%shift
      j=1
    END IF
  END IF
  END DO
  if (j == 1) print FMT_LEND
END SUBROUTINE CHECK_MODIFIERS

END PROGRAM SUPERMODEL
