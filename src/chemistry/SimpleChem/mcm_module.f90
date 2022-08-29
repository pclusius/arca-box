MODULE KPP_ROOT_Main
  ! A driver option for calling KPP as a module, from an external program
  ! by Sampo Smolander, 2008

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: KPP_SetUp, KPP_Proceed

  KPP_REAL :: T, DVAL(NSPEC)
  KPP_REAL :: RSTATE(20)
  INTEGER :: ICNTRL(20) ! Added by Sampo Smolander
  INTEGER :: i
  REAL(DP), SAVE :: R_F(NREACT) = 1d0

CONTAINS

!~~~> Initialization

  SUBROUTINE KPP_SetUp(R_F_in)
    REAL(DP), OPTIONAL :: R_F_in(NREACT)
    IF (PRESENT(R_F_in)) R_F = R_F_in

    STEPMIN = 0.0d0
    STEPMAX = 0.0d0
    ICNTRL = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    DO i=1,NVAR
       RTOL(i) = 1.0d-4
       ATOL(i) = 1.0d-3
    END DO
    CALL Initialize()
  END SUBROUTINE KPP_SetUp

!~~~> Time loop

  SUBROUTINE KPP_Proceed(CONS, T1_in, T2_in, TEMP_in, O2_in, N2_in, M_in, H2O_in, RES1_in, RES2_in, J_in)

    IMPLICIT NONE
    REAL(kind=dp), INTENT(INOUT) :: CONS(NSPEC)
    REAL(kind=dp), INTENT(IN) :: T1_in, T2_in
    REAL(kind=dp), INTENT(IN) :: TEMP_in, O2_in, N2_in, M_in, H2O_in
    REAL(kind=dp), INTENT(IN) :: RES1_in, RES2_in
    REAL(kind=dp), INTENT(IN) :: J_in(NPHOT)

    REAL(kind=dp) :: T1,T2

    T1 = T1_in
    T2 = T2_in

    TEMP  = TEMP_in
    O2    = O2_in
    N2    = N2_in
    M     = M_in
    H2O   = H2O_in
    RES1  = RES1_in
    RES2  = RES2_in
    J     = J_in

    ! Update_RCONST and INTEGRATE will operate on C (actually on VAR and FIX)
    C = CONS

    CALL Update_RCONST()

    ! Modify rate constants
    RCONST = RCONST * R_F

    ! In INTEGRATE VAR and FIX are calculated, so we need to put values from C to them
    ! And FIX is meaningful only when NVAR < NSPEC
    VAR(1:NVAR) = C(1:NVAR)
    DO i = NVAR+1, NSPEC
      FIX(i-NVAR) = C(i)
    END DO

    CALL INTEGRATE( TIN = T1, TOUT = T2, RSTATUS_U = RSTATE, ICNTRL_U = ICNTRL )

    ! Set VAR and FIX back to C
    C(1:NVAR) = VAR(1:NVAR)
    DO i = NVAR+1, NSPEC
      C(i) = FIX(i-NVAR)
    END DO

    CONS = C

  END SUBROUTINE KPP_Proceed


END MODULE KPP_ROOT_Main
