MODULE SOLVEBASES
  USE CONSTANTS
  USE INPUT
  USE second_Precision
  USE ACDC_NH3
  USE ACDC_DMA

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 2 ! Number of elements in the test matrix
  REAL(dp) :: testvalues(2)
  REAL(dp) :: outJ(2)
  INTEGER  :: i
  CHARACTER(110) :: buf
  PRIVATE
  PUBLIC :: Get_BASE


contains


  !---------------------------------------------------------------------------------------------------------------------
  ! Subroutine solves the amount of base needed to produce observed formation rates. BASE here is considered to consist
  ! of Ammonia and DMA, so that BASE is assigned to ammonia, and BASE * DMA_f is assigned to DMA. So if the assumption
  ! is that DMA trend follows ammonia (say some unknown [N]) but concentrations are 5% of ammonia's, then
  ! [D]=[N]*DMA_f = 0.05*[N]. The BASE would then be the concentration of ammonia + 5% or BASE = [N] + 5%*[N]. The
  ! precision of the iteration is defined in resolve_BASE_precision, which is included in NML_MISC and can therefore be
  ! sent in from INITFILE. Also DMA_f is sent in from NML_MISC.
  !---------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Get_BASE(TSTEP_CONC, RESOLVED_BASE, RESOLVED_J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  RESOLVED_BASE ! Resolved concentration of whatever is needed to produce target J
    real(dp), INTENT(INOUT) ::  RESOLVED_J    ! Resolved J with resolved BASE
    real(dp), INTENT(IN)    ::  TSTEP_CONC(:) ! Sulfuric acid in #/cm^3
    real(dp)                ::  SA            ! Sulfuric acid in #/cm^3
    real(dp)                ::  target_J             ! Formation rate in #/s/cm^3
    real(dp)                ::  T             ! Temperature in K
    real(dp)                ::  CS            ! Condensation sink in 1/s
    real(dp)                ::  IPR           ! Ion production rate, ion pairs/s/cm^3
    real(dp)                ::  testJ, testbase, dummy1, dummy2, dummy3(3), maxb, minb,iterJ


    testvalues = ([1d7, 1d12]) * 1d6
    target_J = (TSTEP_CONC(inm_JDMA) + TSTEP_CONC(inm_JIN)) * 1d6
    SA  = TSTEP_CONC(inm_H2SO4) * 1d6
    T   = TSTEP_CONC(inm_TEMPK)
    CS  = TSTEP_CONC(inm_CS)
    IPR = TSTEP_CONC(inm_IPR) * 1d6

    IF (UCASE(Fill_formation_with) == 'DMA') THEN
      target_J = max(0d0, target_J - J_ACDC_NH3*1d6)
      testvalues = testvalues*0.1d0
    ELSE IF  (UCASE(Fill_formation_with) == 'NH3') THEN
      target_J = max(0d0, target_J - J_ACDC_DMA*1d6)
    END IF

    i=0
    ! print'(2(" ",a, es12.3))', 'min NH3', testvalues(1), 'min DMA', DMA_f*(MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f)))
    ! print'(2(" ",a, es12.3))', 'max NH3', testvalues(2), 'max DMA', testvalues(2)*DMA_f
    outJ(1) = JACDC(MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f)))
    outJ(2) = JACDC(testvalues(2))

    ! print'(2(" ",a, es12.3))', 'min J', outJ(1), 'max J', outJ(2)
    if (MODELTIME%printnow) THEN
      print FMT_HDR
      write(buf, '(a,f5.2,a)') 'SOLVING FOR BASE NEEDED TO PRODUCE OBSERVED FORMATION RATE (within ',100*resolve_BASE_precision,'% of target):'
      print FMT_MSG, buf
    END IF

    if (target_J < outJ(1)) THEN
      write(buf,'(a,es9.3,a)') 'Target J (',target_J,') is smaller than what lower limits produce'
      print FMT_MSG, buf
      testbase = 0

    else if (target_J > outJ(2)) THEN
      ! print'(a,es8.2)', 'max target_J :', outJ(2)
      print'(es12.3,a,es12.3,es12.3)', target_J,' Target J is too large for current upper limits of bases', testvalues(2), testvalues(2)*DMA_f
      testbase = C_AIR_NOW*1d6*2

    else
      minb = LOG10((MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f))))
      maxb = LOG10(testvalues(N))
      iterJ = outJ(2)

      ! Binary search to resolve the needed concentration
      do while (.not. ((iterJ/target_J < 1+resolve_BASE_precision) .and. (iterJ/target_J > 1-resolve_BASE_precision)))
        testbase = 10**((maxb + minb)/2)
        iterJ = JACDC(testbase)

        if (iterJ > target_J) THEN
          maxb = LOG10(testbase)
        ELSE
          minb = log10(testbase)
        END IF

      END DO
    end if

    if (MODELTIME%printnow) THEN
      if (UCASE(Fill_formation_with) == '') THEN
        write(buf,'(4(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. NH3 ',testbase*1d-6, ' DMA:', testbase*DMA_f*1d-6, ' (',i,' iterations)'
      ELSE IF (UCASE(Fill_formation_with) == 'DMA') THEN
        write(buf,'(3(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. DMA:', testbase*1d-6, ' (',i,' iterations)'
      ELSE IF (UCASE(Fill_formation_with) == 'NH3') THEN
        write(buf,'(3(a,es8.2),a,i0,a)') 'target J: ',target_J*1d-6, ' resolved J: ', iterJ*1d-6, '. NH3:', testbase*1d-6, ' (',i,' iterations)'
      END IF
      print FMT_MSG, buf
      print FMT_LEND
    END IF

    RESOLVED_J = iterJ*1d-6
    RESOLVED_BASE = testbase*1d-6

  contains
    real(dp) function JACDC(test)
      IMPLICIT NONE
      real(DP), intent(in) :: test
      real(DP) :: BASE
      BASE = test
      IF (UCASE(Fill_formation_with) /= 'DMA') THEN
        i = i+1
        CALL get_acdc_J(SA,BASE,dummy1,CS,T,IPR,MODELTIME,.true.,testJ,dummy2,dummy3)
        JACDC = testJ
      END IF
      IF (UCASE(Fill_formation_with) /= 'NH3') THEN
        if (DMA_f > 0) THEN
          i = i+1
          BASE = test*DMA_f
          IF (UCASE(Fill_formation_with) == 'DMA') THEN
            BASE = test
            JACDC = 0d0
          END IF
          CALL get_acdc_D(SA,BASE,dummy1,CS,T,MODELTIME,.true.,testJ,dummy2)
          JACDC = JACDC + testJ
        END IF
      END IF

    end function JACDC

  END SUBROUTINE Get_BASE

END MODULE SOLVEBASES
