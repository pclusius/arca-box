MODULE SOLVEBASES
  USE CONSTANTS
  USE INPUT
  USE second_Precision
  USE ACDC_NH3
  USE ACDC_DMA

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 2 ! Number of elements in the test matrix
  REAL(dp) :: testvalues(2) = ([1d7, 1d12]) * 1d6
  REAL(dp) :: outJ(2)
  INTEGER  :: i
  CHARACTER(110) :: buf
  PRIVATE
  PUBLIC :: Get_BASE


contains

  SUBROUTINE Get_BASE(TSTEP_CONC)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  TSTEP_CONC(:) ! Sulfuric acid in #/cm^3
    real(dp)                ::  SA            ! Sulfuric acid in #/cm^3
    real(dp)                ::  J             ! Formation rate in #/s/cm^3
    real(dp)                ::  T             ! Temperature in K
    real(dp)                ::  CS            ! Condensation sink in 1/s
    real(dp)                ::  IPR           ! Ion production rate, ion pairs/s/cm^3
    real(dp)                ::  testJ, testbase, dummy1, dummy2, dummy3(3), maxb, minb,iterJ

    SA    = TSTEP_CONC(inm_H2SO4) * 1d6
    J     = (TSTEP_CONC(inm_JDMA) + TSTEP_CONC(inm_JNH3)) * 1d6
    T     = TSTEP_CONC(inm_TEMPK)
    CS    = TSTEP_CONC(inm_CS)
    IPR   = TSTEP_CONC(inm_IPR) * 1d6
    resolve_NH3_precision = 0.005
    i=0
    ! print'(2(" ",a, es12.3))', 'min NH3', testvalues(1), 'min DMA', DMA_f*(MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f)))
    ! print'(2(" ",a, es12.3))', 'max NH3', testvalues(2), 'max DMA', testvalues(2)*DMA_f
    outJ(1) = JACDC(MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f)))
    outJ(2) = JACDC(testvalues(2))
    ! print'(2(" ",a, es12.3))', 'min J', outJ(1), 'max J', outJ(2)
    if (MODELTIME%printnow) THEN
      print FMT_HDR
      write(buf, '(a,f5.2,a)') 'SOLVING FOR BASE NEEDED TO PRODUCE OBSERVED FORMATION RATE (within ',100*resolve_NH3_precision,'% of target):'
      print FMT_MSG, buf
    END IF

    if (J < outJ(1)) THEN
      print'(es12.3,a)', J,' Target J is smaller than what lower limits produce'
      testbase = 0
    else if (J > outJ(2)) THEN
      ! print'(a,es8.2)', 'max J :', outJ(2)
      print'(es12.3,a,es12.3,es12.3)', J,' Target J is too large for current upper limits of bases', testvalues(2), testvalues(2)*DMA_f
      testbase = C_AIR_NOW*1d6*2
    else
      minb = LOG10((MIN(testvalues(1), testvalues(1)/max(1d0,DMA_f))))
      maxb = LOG10(testvalues(N))
      iterJ = outJ(2)
      do while (.not. ((iterJ/J < 1+resolve_NH3_precision) .and. (iterJ/J > 1-resolve_NH3_precision)))
        testbase = 10**((maxb + minb)/2)
        iterJ = JACDC(testbase)
        if (iterJ > J) THEN
          maxb = LOG10(testbase)
        ELSE
          minb = log10(testbase)
        END IF
      END DO
    end if
    if (MODELTIME%printnow) THEN
      write(buf,'(4(a,es8.2),a,i0,a)') 'target J: ',J*1d-6, ' resolved J: ', iterJ*1d-6, '. NH3 ',testbase*1d-6, ' DMA:', testbase*DMA_f*1d-6, ' (',i,' iterations)'
      print FMT_MSG, buf
      print FMT_LEND
    END IF
    RESOLVED_BASE = testbase*1d-6



  contains
    real(dp) function JACDC(test)
      IMPLICIT NONE
      real(DP), intent(in) :: test
      real(DP) :: BASE
      BASE = test
      i = i+1
      CALL get_acdc_J(SA,BASE,dummy1,CS,T,IPR,MODELTIME,.true.,testJ,dummy2,dummy3)
      JACDC = testJ
      if (DMA_f > 0) THEN
        i = i+1
        BASE = test*DMA_f
        CALL get_acdc_D(SA,BASE,dummy1,CS,T,MODELTIME,.true.,testJ,dummy2)
        JACDC = JACDC + testJ
      END IF
    end function JACDC


  END SUBROUTINE Get_BASE
END MODULE SOLVEBASES
