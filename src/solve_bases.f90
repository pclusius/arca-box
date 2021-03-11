MODULE SOLVEBASES
  USE CONSTANTS
  USE AUXILLARIES
  USE INPUT
  USE second_Precision
  USE ACDC_NH3
  USE ACDC_DMA

  IMPLICIT NONE
  REAL(dp) :: conc_limits(2), J_limits(2)
  INTEGER  :: i
  CHARACTER(110) :: buf

  PRIVATE
  PUBLIC :: Get_BASE


contains


  !---------------------------------------------------------------------------------------------------------------------
  ! Subroutine solves the amount of base needed to produce observed formation rates. BASE here is considered to consist
  ! of Ammonia and DMA, so that all BASE is first assigned to ammonia, and BASE * DMA_f is assigned to DMA. The assumption
  ! is that DMA trend follows ammonia (say some unknown [N]) but concentrations are e.g. 5% of ammonia's, then
  ! [D]=[N]*DMA_f = [N]*0.05. The BASE would then be the concentration of ammonia + 5% or BASE = [N]+[D]=[N]+5%*[N]. The
  ! precision of the iteration is defined in resolve_BASE_precision, which is included in NML_MISC and can therefore be
  ! sent in from INITFILE. Precision is the relative difference between target and iterated J: abs(1-Jiter/Jtarget).
  ! Also DMA_f is sent in from NML_MISC.
  !---------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  RESOLVED_BASE ! Resolved concentration of whatever is needed to produce target J
    real(dp), INTENT(INOUT) ::  RESOLVED_J    ! Resolved J with resolved BASE
    real(dp), INTENT(IN)    ::  TSTEP_CONC(:) ! Sulfuric acid in #/cm^3
    real(dp)                ::  SA            ! Sulfuric acid in #/cm^3
    real(dp)                ::  target_J      ! Formation rate in #/s/cm^3
    real(dp)                ::  T             ! Temperature in K
    real(dp)                ::  CS            ! Condensation sink in 1/s
    real(dp)                ::  IPR           ! Ion production rate, ion pairs/s/cm^3
    real(dp)                ::  testbase, dummy1, dummy2, dummy3(3), maxb, minb,iterJ

    ! Fill the work variables and convert to metercubed
    conc_limits = ([1d7, 1d12]) * 1d6
    target_J = TSTEP_CONC(inm_JIN) * 1d6
    SA  = TSTEP_CONC(inm_H2SO4) * 1d6
    T   = TSTEP_CONC(inm_TEMPK)
    CS  = GCS
    IPR = TSTEP_CONC(inm_IPR) * 1d6

    ! If only part of the formation is iterated, here using DMA, the target is the difference between observed and ACDC J
    IF (UCASE(Fill_formation_with) == 'DMA') THEN
      target_J = max(0d0, target_J - J_ACDC_NH3_M3)
      ! limits are shifted down for pure DMA; upper limit is still ample; lower might sometimes be still too much
      conc_limits = conc_limits*0.1d0
    ! If only part of the formation is iterated, here using NH3, the target is the difference between observed and ACDC J
    ELSE IF  (UCASE(Fill_formation_with) == 'NH3') THEN
      target_J = max(0d0, target_J - J_ACDC_DMA_M3)
    END IF

    i=0 ! reset counter for total number of ACDC runs
    iterJ = 0d0
    ! Here the limits for formation is establised. The min/max for dma_f is trying to handle DMA_f = 0 and DMA_f > 1

    J_limits(1) = JACDC(MIN(conc_limits(1), conc_limits(1)/max(1d0,DMA_f)))
    J_limits(2) = JACDC(conc_limits(2))

    IF (J_limits(1) < 1d-100) J_limits(1) = 1d-100

    ! In case someone sends percentages... Not foolproof
    ! IF (resolve_BASE_precision > 1) resolve_BASE_precision = resolve_BASE_precision/100d0

    ! show messages only when other printouts are shown, even if sub is run at filesaves
    if (GTIME%printnow) THEN
      print FMT_HDR
      write(buf, '(a,f5.2,a)') 'SOLVING FOR BASE NEEDED TO PRODUCE OBSERVED FORMATION RATE (within ',100*resolve_BASE_precision,'% of target):'
      print FMT_MSG, TRIM(buf)
    END IF

    ! if target J is smaller than lower limits would produce, don't bother iterating further
    if (target_J < 1d-200) THEN
      if (GTIME%printnow) print FMT_MSG, 'Target J is zero. No need for extra formation.'
      testbase = 0

    ! if target J is smaller than lower limits would produce, don't bother iterating further
    else if (target_J < J_limits(1)) THEN
      if (GTIME%printnow) print FMT_MSG, 'Target J ('//TRIM(f2chr(target_J*1d-6))//') is smaller than what lower limits produce'
      testbase = 0

    ! if target J is larger than upper limits can produce, don't bother iterating further
    else if (target_J > J_limits(2)) THEN
      if (GTIME%printnow) print FMT_MSG, ' Target J '//TRIM(f2chr(target_J*1d-6))//' is too large for current upper limits of bases: '&
                                    //TRIM(f2chr(conc_limits(2)))//' and '//TRIM(f2chr(conc_limits(2)*DMA_f))
      testbase = GC_AIR_NOW*1d6

    ! if target J reasonable, iterate the concentrations
    else
      ! set upper and lower bonds for binary search
      minb = LOG10((MIN(conc_limits(1), conc_limits(1)/max(1d0,DMA_f))))
      maxb = LOG10(conc_limits(2))
      iterJ = J_limits(2)

      ! Binary search to resolve the needed concentration
      do while (.not. ((iterJ/target_J < 1+resolve_BASE_precision) .and. (iterJ/target_J > 1-resolve_BASE_precision)))
        testbase = 10**((maxb + minb)/2)
        iterJ = JACDC(testbase)

        if (iterJ > target_J) THEN
          maxb = LOG10(testbase)
        ELSE
          minb = log10(testbase)
        END IF
        if (i>50) THEN
          PRINT*, 'Situation too instable, not possible to reach the precision, settling for less.'
          EXIT
        END IF
      END DO
    end if
    ! inform the result
    if (GTIME%printnow) THEN
      if (TRIM(UCASE(Fill_formation_with)) == '') THEN
        write(buf,'(4(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. NH3 ',testbase*1d-6, ' DMA:', testbase*DMA_f*1d-6, ' (',i,' iterations)'
      ELSE IF (UCASE(Fill_formation_with) == 'DMA') THEN
        write(buf,'(3(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. DMA:', testbase*1d-6, ' (',i,' iterations)'
      ELSE IF (UCASE(Fill_formation_with) == 'NH3') THEN
        write(buf,'(3(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. NH3:', testbase*1d-6, ' (',i,' iterations)'
      END IF
      print FMT_MSG, TRIM(buf)
      print FMT_LEND
    END IF
    ! save the result
    RESOLVED_J = iterJ*1d-6
    RESOLVED_BASE = testbase*1d-6
    IF (RESOLVED_BASE<0) RESOLVED_BASE = 0d0

  contains

    !-------------------------------------------------------------------------------------------------------------------
    ! Function that handles the ACDC calls, this is basically just for convenience for the sub it is in. This way calling
    ! ACDC is like calling a variable, also different conditions for solver are handled here.
    !-------------------------------------------------------------------------------------------------------------------
    real(dp) function JACDC(test)
      IMPLICIT NONE
      real(DP), intent(in) :: test
      real(DP) :: BASE,testJ ! the cache variable for concentration. BASE is INTENT(INOUT) in ACDC, therefore this expendable var.

      BASE = test
      JACDC = 0d0
      IF (UCASE(Fill_formation_with) /= 'DMA') THEN
        i = i+1
        CALL get_acdc_J(SA,BASE,dummy1,CS,T,IPR,GTIME,.true.,testJ,dummy2,dummy3,GTIME%dt)
        JACDC = testJ
      END IF
      IF (UCASE(Fill_formation_with) /= 'NH3') THEN
        if (DMA_f > 0 .or. UCASE(Fill_formation_with) == 'DMA') THEN
          i = i+1
          BASE = test*DMA_f
          IF (UCASE(Fill_formation_with) == 'DMA') THEN
            BASE = test
            JACDC = 0d0
          END IF
          CALL get_acdc_D(SA,BASE,dummy1,CS,T,GTIME,.true.,testJ,dummy2,GTIME%dt)
          JACDC = JACDC + testJ
          if (JACDC>0d0) RESOLVED_J_FACTR = testJ/JACDC
        END IF
      END IF

    end function JACDC

  END SUBROUTINE Get_BASE

END MODULE SOLVEBASES
