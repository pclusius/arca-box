MODULE AUXILLARIES
USE SECOND_PRECISION,  ONLY : dp, ik
USE CONSTANTS
IMPLICIT NONE
PUBLIC

INTERFACE TIME_HMS
  module procedure time_hms4
  module procedure time_hms8
  module procedure time_hmsi
end INTERFACE TIME_HMS

CONTAINS

  ! Function that will convert second in CHARACTER(8) of form hh:mm:ss
  ! Due to overloading accepts all numeric types
  CHARACTER(8) function time_hmsi(s)
    INTEGER :: s
    time_hmsi = time___HMS(real(s))
  end function time_hmsi
  CHARACTER(8) function time_hms4(s)
    real :: s
    time_hms4 = time___HMS(s)
  end function time_hms4
  CHARACTER(8) function time_hms8(s)
    real(dp) :: s
    time_hms8 = time___HMS(real(s))
  end function time_hms8
  CHARACTER(8) function time___HMS(s)
    real :: s
    ! This will return the character string - modify if necessary
    write(time___HMS, '(i2.2, ":" i2.2, ":" i2.2)') int(s+0.5_dp)/3600_ik, &
    int(MODULO(int(s+0.5_dp),3600_ik)/60_ik), MODULO(MODULO(int(s+0.5_dp),3600_ik), 60_ik)
  end function time___HMS


  ! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/cm3
  real(dp) function C_AIR_cc(T,P)
    REAL(dp) :: T, P
    C_AIR_cc = P/(kb*T)*1d-6
  end function C_AIR_cc
  ! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/m3
  real(dp) function C_AIR_m3(P,T)
    REAL(dp) :: T, P
    C_AIR_m3 = P/(kb*T)
  end function C_AIR_m3
  ! Converts hours to seconds
  real(dp) function hrs_to_s(h)
    REAL(dp) :: h
    hrs_to_s = h*3600
  end function hrs_to_s

  real(dp) function sec_to_d(s)
    REAL(dp) :: s
    sec_to_d = s/24d0/3600d0
  end function sec_to_d

! Checks if model time is at some exact minute interval
! t = time in seconds (real8)
! check = minute that is to be checked (integer)
logical function EVENMIN(t,check)
  INTEGER :: check
  REAL(dp):: t
  IF (MODULO(int(t+0.5d0), 60*check) == 0) THEN
    EVENMIN = .true.
  ELSE
    EVENMIN = .false.
  END IF
end function EVENMIN

!==============================================================================
! Function to return y-value at [time] from a normal distribution function with
! standard deviation [sigma], minimum value (DC-level) [mn], maximum value [mx]
! and time of maximum at [mean]. Combine with mean = mean + SIN(time/3600)*k to
! model many realistic diurnal concentrations, or use large deviation with late
! mean to model for example VOC etc. If LG=TRUE, will scale y logaritmically
!..............................................................................
REAL(dp) FUNCTION NORMALD(time,sigma, mn, mx, mean, LG)
  REAL(dp) :: time,sigma, mn, mx, mean
  REAL(dp) :: f
  LOGICAL  :: LG
  NORMALD = 1d0/SQRT(2d0*pi*sigma**2) * EXP(- (time/3600d0-mean)**2/(2d0*sigma**2))
  f = 1d0/SQRT(2d0*pi*sigma**2)
  if (LG) THEN
    f = (LOG10(mx-mn))/f
    NORMALD = 10**(NORMALD*f) + mn
  ELSE
    f = (mx-mn)/f
    NORMALD = NORMALD*f + mn
  END IF
end FUNCTION NORMALD

end MODULE AUXILLARIES
