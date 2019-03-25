MODULE AUXILLARIES
USE SECOND_PRECISION,  ONLY : dp
USE CONSTANTS
IMPLICIT NONE
PUBLIC

CHARACTER(56) :: FMT_Tm       = '("+.",t18,82("."),t3,"Time: ",a,t100, "+")'
CHARACTER(56) :: FMT10_CVU   = '("| ",t3, a,t13, es10.3,a,t100, "|")'
CHARACTER(86) :: FMT10_2CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT_LEND    = '("+",t2, 98(".") t100, "+")'
CHARACTER(86) :: FMT_SUB    = '(t5,":",8(".") a)'

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
  real(dp) function C_AIR_m3(T,P)
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
  IMPLICIT NONE
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
! and time of maximum at [mean]. Combine with mean = mean + SIN(time/3600*c)*k to
! model many realistic diurnal concentrations, or use large deviation with late
! mean to model for example VOC etc. If LG=TRUE, will scale y logaritmically
! If minimum >= maximum, minumum value is used as constant concentration
!..............................................................................
REAL(dp) FUNCTION NORMALD(time,MODS)
  IMPLICIT NONE
  type(parametered_input) :: MODS
  REAL(dp) :: time
  REAL(dp) :: f, D
  D = MODS%peaktime + sin(time/3600*MODS%omega)*MODS%amplitude
  NORMALD = 1d0/SQRT(2d0*pi*MODS%sigma**2) * EXP(- (time/3600d0-D)**2/(2d0*MODS%sigma**2))
  f = 1d0/SQRT(2d0*pi*MODS%sigma**2)
  ! Check if minimum is same or more than maximum and if so, use constant concentration (minumum)
  IF (MODS%max_c-MODS%min_c <1e-9) THEN
    NORMALD = MODS%min_c
  ELSE
    if (MODS%LOGSCALE) THEN
      f = (LOG10(MODS%max_c-MODS%min_c+1))/f
      NORMALD = 10**(NORMALD*f)-1 + MODS%min_c
    ELSE
      f = (MODS%max_c-MODS%min_c)/f
      NORMALD = NORMALD*f + MODS%min_c
    END IF
  END IF
end FUNCTION NORMALD

!==============================================================================
! Function to return y-value at [time] from a sine function. Use to model a
! periodic variable such as temperature or relative humidity. By default uses
! 24 hour period.
!..............................................................................
REAL(dp) FUNCTION SIND(time,MODS)
  IMPLICIT NONE
  type(parametered_input) :: MODS
  REAL(dp) :: time, A
  A = (MODS%max_c-MODS%min_c)/2d0
  SIND = A*SIN(2d0*pi*time/3600d0/24d0 - 3d0*pi/4d0) + MODS%min_c + A
end FUNCTION SIND

end MODULE AUXILLARIES
